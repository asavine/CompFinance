
/*
Written by Antoine Savine in 2018

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications
is freely granted to any person or company who purchased a copy of the book

Modern Computational Finance: AAD and Parallel Simulations
Antoine Savine
Wiley, 2018

As long as this comment is preserved at the top of the file
*/

//  Entry points to the library

#pragma once

#include "mcBase.h"
#include "mcMdl.h"
#include "mcPrd.h"
#include "mrg32k3a.h"
#include "sobol.h"
#include <numeric>
#include <fstream>
using namespace std;

#include "store.h"

struct NumericalParam
{
    bool              parallel;
    bool              useSobol;
    int               numPath;
    int               seed1 = 12345;
    int               seed2 = 1234;
};

//  Price product in model
inline auto value(
    const Model<double>&    model,
    const Product<double>&  product,
    //  numerical parameters
    const NumericalParam&   num)
{
    //  Random Number Generator
    unique_ptr<RNG> rng;
    if (num.useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(num.seed1, num.seed2);

    //  Simulate
    const auto resultMat = num.parallel
        ? mcParallelSimul(product, model, *rng, num.numPath)
        : mcSimul(product, model, *rng, num.numPath);

    //  We return 2 vectors : the payoff identifiers and their values
    struct
    {
        vector<string> identifiers;
        vector<double> values;
    } results;

    const size_t nPayoffs = product.payoffLabels().size();
    results.identifiers = product.payoffLabels();
    results.values.resize(nPayoffs);
    for (size_t i = 0; i < nPayoffs; ++i)
    {
        results.values[i] = accumulate(resultMat.begin(), resultMat.end(), 0.0,
            [i](const double acc, const vector<double>& v) { return acc + v[i]; }
        ) / num.numPath;
    }

    return results;
}

//  Overload that picks product and model by name in the store
inline auto value(
    const string&           modelId,
    const string&           productId,
    //  numerical parameters
    const NumericalParam&   num)
{
    //  Get model and product
    const Model<double>* model = getModel<double>(modelId);
    const Product<double>* product = getProduct<double>(productId);

    if (!model || !product)
    {
        throw runtime_error("value() : Could not retrieve model and product");
    }

    return value(*model, *product, num);
}

//  AAD risk, one payoff
inline auto AADriskOne(
    const string&           modelId,
    const string&           productId,
    const NumericalParam&   num,
    const string&           riskPayoff = "")
{
    //  Get model and product
    const Model<Number>* model = getModel<Number>(modelId);
    const Product<Number>* product = getProduct<Number>(productId);

    if (!model || !product)
    {
        throw runtime_error("AADrisk() : Could not retrieve model and product");
    }

    //  Random Number Generator
    unique_ptr<RNG> rng;
    if (num.useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(num.seed1, num.seed2);

    //  Find the payoff for risk
    size_t riskPayoffIdx = 0;
    if (!riskPayoff.empty())
    {
        const vector<string>& allPayoffs = product->payoffLabels();
        auto it = find(allPayoffs.begin(), allPayoffs.end(), riskPayoff);
        if (it == allPayoffs.end())
        {
            throw runtime_error("AADriskOne() : payoff not found");
        }
        riskPayoffIdx = distance(allPayoffs.begin(), it);
    }

    //  Simulate
    const auto simulResults = num.parallel
        ? mcParallelSimulAAD(*product, *model, *rng, num.numPath,
            [riskPayoffIdx](const vector<Number>& v) {return v[riskPayoffIdx]; })
        : mcSimulAAD(*product, *model, *rng, num.numPath,
            [riskPayoffIdx](const vector<Number>& v) {return v[riskPayoffIdx]; });

    //  We return: a number and 2 vectors : 
    //  -   The payoff identifiers and their values
    //  -   The value of the aggreagte payoff
    //  -   The parameter idenitifiers 
    //  -   The sensititivities of the aggregate to parameters
    struct
    {
        vector<string>  payoffIds;
        vector<double>  payoffValues;
        double          riskPayoffValue;
        vector<string>  paramIds;
        vector<double>  risks;
    } results;

    const size_t nPayoffs = product->payoffLabels().size();
    results.payoffIds = product->payoffLabels();
    results.payoffValues.resize(nPayoffs);
    for (size_t i = 0; i < nPayoffs; ++i)
    {
        results.payoffValues[i] = accumulate(
            simulResults.payoffs.begin(), 
            simulResults.payoffs.end(), 
            0.0,
            [i](const double acc, const vector<double>& v) { return acc + v[i]; }
        ) / num.numPath;
    }
    results.riskPayoffValue = accumulate(
        simulResults.aggregated.begin(),
        simulResults.aggregated.end(),
        0.0) / num.numPath;
    results.paramIds = model->parameterLabels();
    results.risks = move (simulResults.risks);

    return results;
}

//  AAD risk, aggregate portfolio
inline auto AADriskAggregate(
    const string&           modelId,
    const string&           productId,
    const map<string, double>&   notionals,
    const NumericalParam&   num)
{
    //  Get model and product
    const Model<Number>* model = getModel<Number>(modelId);
    const Product<Number>* product = getProduct<Number>(productId);

    if (!model || !product)
    {
        throw runtime_error("AADriskAggregate() : Could not retrieve model and product");
    }

    //  Random Number Generator
    unique_ptr<RNG> rng;
    if (num.useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(num.seed1, num.seed2);

    //  Vector of notionals
    const vector<string>& allPayoffs = product->payoffLabels();
    vector<double> vnots(allPayoffs.size(), 0.0);
    for (const auto& notional : notionals)
    {
        auto it = find(allPayoffs.begin(), allPayoffs.end(), notional.first);
        if (it == allPayoffs.end())
        {
            throw runtime_error("AADriskAggregate() : payoff not found");
        }
        vnots[distance(allPayoffs.begin(), it)] = notional.second;
    }

    //  Aggregator
    auto aggregator = [&vnots](const vector<Number>& payoffs)
    {
        return inner_product(payoffs.begin(), payoffs.end(), vnots.begin(), Number(0.0));
    };

    //  Simulate
    const auto simulResults = num.parallel
        ? mcParallelSimulAAD(*product, *model, *rng, num.numPath, aggregator)
        : mcSimulAAD(*product, *model, *rng, num.numPath, aggregator);

    //  We return: a number and 2 vectors : 
    //  -   The payoff identifiers and their values
    //  -   The value of the aggreagte payoff
    //  -   The parameter idenitifiers 
    //  -   The sensititivities of the aggregate to parameters
    struct
    {
        vector<string>  payoffIds;
        vector<double>  payoffValues;
        double          riskPayoffValue;
        vector<string>  paramIds;
        vector<double>  risks;
    } results;

    const size_t nPayoffs = product->payoffLabels().size();
    results.payoffIds = product->payoffLabels();
    results.payoffValues.resize(nPayoffs);
    for (size_t i = 0; i < nPayoffs; ++i)
    {
        results.payoffValues[i] = accumulate(
            simulResults.payoffs.begin(),
            simulResults.payoffs.end(),
            0.0,
            [i](const double acc, const vector<double>& v) { return acc + v[i]; }
        ) / num.numPath;
    }
    results.riskPayoffValue = accumulate(
        simulResults.aggregated.begin(),
        simulResults.aggregated.end(),
        0.0) / num.numPath;
    results.paramIds = model->parameterLabels();
    results.risks = move(simulResults.risks);

    return results;
}

//  Returns a vector of values and a matrix of risks 
//      with payoffs in columns and parameters in rows
//      along with ids of payoffs and parameters

struct RiskReports
{
    vector<string> payoffs;
    vector<string> params;
    vector<double> values;
    matrix<double> risks;
};

//  Itemized AAD risk, one per payoff
inline RiskReports AADriskMulti(
    const string&           modelId,
    const string&           productId,
    const NumericalParam&   num)
{
    const Model<Number>* model = getModel<Number>(modelId);
    const Product<Number>* product = getProduct<Number>(productId);

    if (!model || !product)
    {
        throw runtime_error("AADrisk() : Could not retrieve model and product");
    }

    RiskReports results;

    //  Random Number Generator
    unique_ptr<RNG> rng;
    if (num.useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(num.seed1, num.seed2);

    //  Simulate
    const auto simulResults = num.parallel
		? mcParallelSimulAADMulti(*product, *model, *rng, num.numPath)
        : mcSimulAADMulti(*product, *model, *rng, num.numPath);

    results.params = model->parameterLabels();
    results.payoffs = product->payoffLabels();
	results.risks = move(simulResults.risks);

	//	Average values across paths
	const size_t nPayoffs = product->payoffLabels().size();
	results.values.resize(nPayoffs);
	for (size_t i = 0; i < nPayoffs; ++i)
	{
		results.values[i] = accumulate(
			simulResults.payoffs.begin(),
			simulResults.payoffs.end(),
			0.0,
			[i](const double acc, const vector<double>& v) { return acc + v[i]; }
		) / num.numPath;
	}

    return results;
}

//  Bump risk, itemized
//  Same result format as AADriskMulti()
inline RiskReports bumpRisk(
    const string&           modelId,
    const string&           productId,
    const NumericalParam&   num)
{
    auto* orig = getModel<double>(modelId);
    const Product<double>* product = getProduct<double>(productId);

    if (!orig || !product)
    {
        throw runtime_error("bumpRisk() : Could not retrieve model and product");
    }

    RiskReports results;

    //  base values
    auto baseRes = value(*orig, *product, num);
    results.payoffs = baseRes.identifiers;
	results.values = baseRes.values;

    //  make copy so we don't modify the model in memory
    auto model = orig->clone();
    
    results.params = model->parameterLabels();
    const vector<double*> parameters = model->parameters();
    const size_t n = parameters.size(), m = results.payoffs.size();
    results.risks.resize(n, m);

    //  bumps
    for (size_t i = 0; i < n; ++i)
    {
        *parameters[i] += 1.e-08;
        auto bumpRes = value(*model, *product, num);
        *parameters[i] -= 1.e-08;

        for (size_t j = 0; j < m; ++j)
        {
            results.risks[i][j] = 1.0e+08 *
                (bumpRes.values[j] - baseRes.values[j]);
        }
    }

    return results;
}

//  Dupire specific

//  Returns a struct with price, delta and vega matrix
inline auto dupireAADRisk(
    //  model id
    const string&           modelId,
    //  product id
    const string&           productId,
    const map<string, double>&   notionals,
    //  numerical parameters
    const NumericalParam&   num)
{
    //  Check that the model is a Dupire
    const Model<Number>* model = getModel<Number>(modelId);
    if (!model)
    {
        throw runtime_error("dupireAADRisk() : Model not found");
    }
    const Dupire<Number>* dupire = dynamic_cast<const Dupire<Number>*>(model);
    if (!dupire)
    {
        throw runtime_error("dupireAADRisk() : Model not a Dupire");
    }

    //  Results
    struct
    {
        double value;
        double delta;
        matrix<double> vega;
    } results;

    //  Go
    auto simulResults = AADriskAggregate(modelId, productId, notionals, num);

    //  Find results

    //  Value
    results.value = simulResults.riskPayoffValue;

    //  Delta

    results.delta = simulResults.risks[0];

    //  Vegas

    results.vega.resize(dupire->spots().size(), dupire->times().size());
    copy(next(simulResults.risks.begin()), simulResults.risks.end(), results.vega.begin());

    return results;
}

//  Returns spots, times and lVols in a struct
inline auto
dupireCalib(
    //  The local vol grid
    //  The spots to include
    const vector<double>& inclSpots,
    //  Maximum space between spots
    const double maxDs,
    //  The times to include, note NOT 0
    const vector<Time>& inclTimes,
    //  Maximum space between times
    const double maxDt,
    //  The IVS we calibrate to
    //  'B'achelier, Black'S'choles or 'M'erton
    const double spot,
    const double vol,
    const double jmpIntens = 0.0,
    const double jmpAverage = 0.0,
    const double jmpStd = 0.0)
{
    //  Create IVS
    MertonIVS ivs(spot, vol, jmpIntens, jmpAverage, jmpStd);

    //  Go
    return dupireCalib(ivs, inclSpots, maxDs, inclTimes, maxDt);
}

//  Superbucket

struct SuperbucketResults
{
    double value;
    double delta;
    vector<double> strikes;
    vector<Time> mats;
    matrix<double> vega;
};

//  Returns value, delta, strikes, maturities 
//      and vega = derivatives to implied vols = superbucket
inline auto
    dupireSuperbucket(
    //  Model parameters that are not calibrated
    const double            spot,
    const double            maxDt,
    //  Product 
    const string&           productId,
    const map<string, double>&   notionals,
    //  The local vol grid
    //  The spots to include
    const vector<double>&   inclSpots,
    //  Maximum space between spots
    const double            maxDs,
    //  The times to include, note NOT 0
    const vector<Time>&     inclTimes,
    //  Maximum space between times
    const double            maxDtVol,
    //  The IVS we calibrate to
    //  Risk view
    const vector<double>&   strikes,
    const vector<Time>&     mats,
    //  Merton params
    const double            vol,
    const double            jmpIntens,
    const double            jmpAverage,
    const double            jmpStd,
    //  Numerical parameters
    const NumericalParam&   num)
{
    //  Results
    SuperbucketResults results;

    //  Start with a clean tape
    auto* tape = Number::tape;
    tape->rewind();

    //  Calibrate the model
    auto params = dupireCalib(
        inclSpots, 
        maxDs, 
        inclTimes, 
        maxDtVol, 
        spot, 
        vol, 
        jmpIntens, 
        jmpAverage, 
        jmpStd);
    const vector<double>& spots = params.spots;
    const vector<Time>& times = params.times;
    const matrix<double>& lvols = params.lVols;

    //  Put in memory
    putDupire(spot, spots, times, lvols, maxDt, "superbucket");

    //  Find delta and microbucket
    auto mdlDerivs = dupireAADRisk(
        "superbucket",
        productId,
        notionals,
        num);
    results.value = mdlDerivs.value;
    results.delta = mdlDerivs.delta;
    const matrix<double>& microbucket = mdlDerivs.vega;

    //  Clear tape
    //  tape->rewind();
    tape->clear();

    //  Convert market inputs to numbers, put on tape
            
    //  Create IVS
    MertonIVS ivs(spot, vol, jmpIntens, jmpAverage, jmpStd);
    
    //  Risk view --> that is the AAD input
    //  Note: that puts the view on tape
    RiskView<Number> riskView(strikes, mats);

    //  Calibrate again, in AAD mode, make tape
    auto nParams = dupireCalib(
        ivs, 
        inclSpots, 
        maxDs, 
        inclTimes, 
        maxDtVol, 
        riskView);
    matrix<Number>& nLvols = nParams.lVols;

    //  Seed local vol adjoints on tape with microbucket results
    for (size_t i = 0; i < microbucket.rows(); ++i)
    {
        for (size_t j = 0; j < microbucket.cols(); ++j)
        {
            nLvols[i][j].adjoint() = microbucket[i][j];
        }
    }
    
    //  Propagate
    Number::propagateAdjoints(prev(tape->end()), tape->begin());

    //  Results: superbucket = risk view

    //  Copy results
    results.strikes = strikes;
    results.mats = mats;
    results.vega.resize(riskView.rows(), riskView.cols());
    transform(riskView.begin(), riskView.end(), results.vega.begin(), 
        [](const Number& n)
    {
        return n.adjoint();
    });

    //  Clear tape
    tape->clear();

    //  Return results
    return results;
}

//  Superbucket with bumps

//  Returns value, delta, strikes, maturities 
//      and vega = derivatives to implied vols = superbucket
inline auto
    dupireSuperbucketBump(
        //  Model parameters that are not calibrated
        const double            spot,
        const double            maxDt,
        //  Product 
        const string&           productId,
        const map<string, double>&   notionals,
        //  The local vol grid
        //  The spots to include
        const vector<double>&   inclSpots,
        //  Maximum space between spots
        const double            maxDs,
        //  The times to include, note NOT 0
        const vector<Time>&     inclTimes,
        //  Maximum space between times
        const double            maxDtVol,
        //  The IVS we calibrate to
        //  Risk view
        const vector<double>&   strikes,
        const vector<Time>&     mats,
        //  Merton params
        const double            vol,
        const double            jmpIntens,
        const double            jmpAverage,
        const double            jmpStd,
        //  Numerical parameters
        const NumericalParam&   num)
{
    //  Results
    SuperbucketResults results;

    //  Calibrate the model
    auto params = dupireCalib(
        inclSpots,
        maxDs,
        inclTimes,
        maxDtVol,
        spot,
        vol,
        jmpIntens,
        jmpAverage,
        jmpStd);
    const vector<double>& spots = params.spots;
    const vector<Time>& times = params.times;
    const matrix<double>& lvols = params.lVols;

    //  Create model
    Dupire<double> model(spot, spots, times, lvols, maxDt);
    
    //  Get product
    const Product<double>* product = getProduct<double>(productId);

    //  Base price
    auto baseVals = value(model, *product, num);

    //  Vector of notionals
    const vector<string>& allPayoffs = baseVals.identifiers;
    vector<double> vnots(allPayoffs.size(), 0.0);
    for (const auto& notional : notionals)
    {
        auto it = find(allPayoffs.begin(), allPayoffs.end(), notional.first);
        if (it == allPayoffs.end())
        {
            throw runtime_error("dupireSuperbucketBump() : payoff not found");
        }
        vnots[distance(allPayoffs.begin(), it)] = notional.second;
    }

    //  Base book value
    results.value = inner_product(vnots.begin(), vnots.end(), baseVals.values.begin(), 0.0);

    //  Create IVS
    MertonIVS ivs(spot, vol, jmpIntens, jmpAverage, jmpStd);

    //  Create risk view 
    RiskView<double> riskView(strikes, mats);

    //  Bumps

    //  bump, recalibrate, reset model, reprice, pick value, unbump

    //  Delta

    //  Recreate model
    Dupire<double> bumpedModel(spot + 1.0e-08, spots, times, lvols, maxDt);
    //  Reprice
    auto bumpedVals = value(bumpedModel, *product, num);
    //  Pick results and differentiate
    results.delta = (
        inner_product(vnots.begin(), vnots.end(), bumpedVals.values.begin(), 0.0)
        - results.value) * 1.0e+08;

    //  Vega

    const size_t n = riskView.rows(), m = riskView.cols();
    results.vega.resize(n, m);
    for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < m; ++j)
    {
        //  Bump
        riskView.bump(i, j, 1.0e-05);
        //  Recalibrate
        auto bumpedCalib = dupireCalib(ivs, inclSpots, maxDs, inclTimes, maxDtVol, riskView);
        //  Recreate model
        Dupire<double> bumpedModel(spot, bumpedCalib.spots, bumpedCalib.times, bumpedCalib.lVols, maxDt);
        //  Reprice
        auto bumpedVals = value(bumpedModel, *product, num);
        //  Pick results and differentiate
        results.vega[i][j] = (
            inner_product(vnots.begin(), vnots.end(), bumpedVals.values.begin(), 0.0)
            - results.value) * 1.0e+05;
        //  Unbump
        riskView.bump(i, j, -1.0e-05);
    }

    //  Copy results and strikes
    results.strikes = strikes;
    results.mats = mats;

    //  Return results
    return results;
}