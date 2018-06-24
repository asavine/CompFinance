
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

#pragma once

#include "mcBase.h"
#include "mcMdl.h"
#include "mcPrd.h"
#include "mrg32k3a.h"
#include "sobol.h"
#include <numeric>
#include <fstream>
using namespace std;

inline double uocDupire(
    //  model parameters
    const double            spot,
    const vector<double>&   spots,
    const vector<Time>&     times,
    //  spot major
    const matrix<double>&   vols,   
    const double            maxDt,
    //  product parameters
    const double            strike,
    //  negative = european
    const double            barrier,
    const Time              maturity,
    const double            monitorFreq,
    //  numerical parameters
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    //  optionals
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Build model, product and rng
    Dupire<double> model(spot, spots, times, vols, maxDt);
    
    //  Product
    unique_ptr<Product<double>> product;
    if (barrier > 0) product = make_unique<UOC<double>>(strike, barrier, maturity, monitorFreq);
    else product = make_unique<European<double>>(strike, maturity);

    //  RNG
    unique_ptr<RNG> rng;
    if (useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1, seed2);

    //  Simulate
    const auto resultMat = parallel
        ? mcParallelSimul(*product, model, *rng, numPath)
        : mcSimul(*product, model, *rng, numPath);

    //  Compute averages among paths
    double result = accumulate(resultMat.begin(), resultMat.end(), 0.0,
        [](const double acc, const vector<double>& v) { return acc + v[0]; }
        ) / numPath;
        
    return result;
}

inline vector<vector<double>> europeansDupire(
    //  model parameters
    const double            spot,
    const vector<double>&   spots,
    const vector<Time>&     times,
    //  spot major
    const matrix<double>&   vols,
    const double            maxDt,
    //  product parameters
    const map<Time, vector<double>>&    options, 
    //  numerical parameters
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Build model, product and rng
    Dupire<double> model(spot, spots, times, vols, maxDt);

    //  Product
    unique_ptr<Product<double>> product = make_unique<Europeans<double>>(options);

    //  RNG
    unique_ptr<RNG> rng;
    if (useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1, seed2);

    //  Simulate
    const auto resultMat = parallel
        ? mcParallelSimul(*product, model, *rng, numPath)
        : mcSimul(*product, model, *rng, numPath);

    //  Compute averages among paths
    const Europeans<double>* prd = static_cast<const Europeans<double>*> (product.get());
    vector<vector<double>> results = prd->strikes();   //  just for the right size

    size_t nPay = 0;
    for (size_t i = 0; i < prd->maturities().size(); ++i)
    {
        for (size_t j = 0; j < prd->strikes()[i].size(); ++j)
        {
            results[i][j] = accumulate(resultMat.begin(), resultMat.end(), 0.0,
                [nPay](const double acc, const vector<double>& payoffs)
                {
                    return acc + payoffs[nPay];
                }
            ) / numPath;

            ++nPay;
        }
    }

    return results;
}

//  Returns a struct with price, delta and vega matrix
inline auto uocDupireBumpRisk(
    //  model parameters
    const double            spot,
    const vector<double>&   spots,
    const vector<Time>&     times,
    const matrix<double>&   vols,   //  spot major
    const double            maxDt,
    //  product parameters
    const double            strike,
    const double            barrier,
    const Time              maturity,
    const double            monitorFreq,
    //  numerical parameters
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    //  optionals
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Results
    struct
    {
        double value;
        double delta;
        matrix<double> vega;
    } results;

    //  base value
    results.value = uocDupire(
        spot,
        spots,
        times,
        vols,
        maxDt,
        strike,
        barrier,
        maturity,
        monitorFreq,
        parallel,
        useSobol,
        numPath,
        seed1,
        seed2);
    
    //  Delta
    double newSpot = spot;
    newSpot += 1.e-08;
    const double v1 = uocDupire(
        newSpot,
        spots,
        times,
        vols,
        maxDt,
        strike,
        barrier,
        maturity,
        monitorFreq,
        parallel,
        useSobol,
        numPath,
        seed1,
        seed2);
    results.delta = 1.e+08 * (v1 - results.value);

    //  Vega
    matrix<double> newVols = vols;
    results.vega.resize(vols.rows(), vols.cols());
    for (auto i = 0; i < newVols.rows(); ++i)
    {
        for (auto j = 0; j < newVols.cols(); ++j)
        {
            newVols[i][j] += 1.e-08;
            const double v1 = uocDupire(
                spot,
                spots,
                times,
                newVols,
                maxDt,
                strike,
                barrier,
                maturity,
                monitorFreq,
                parallel,
                useSobol,
                numPath,
                seed1,
                seed2);
            results.vega[i][j] = 1.e+08 * (v1 - results.value);
            newVols[i][j] -= 1.e-08;
        }
    }

    return results;
}

//  Returns a struct with price, delta and vega matrix
inline auto uocDupireAADRisk(
    //  model parameters
    const double            spot,
    const vector<double>&   spots,
    const vector<Time>&     times,
    //  spot major
    const matrix<double>&   vols,   
    const double            maxDt,
    //  product parameters
    const double            strike,
    //  negative = european
    const double            barrier,
    const Time              maturity,
    const double            monitorFreq,
    //  numerical parameters
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    //  optionals
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Results
    struct
    {
        double value;
        double delta;
        matrix<double> vega;
    } results;

    //  Build model, product and rng
    Dupire<Number> model(spot, spots, times, vols, maxDt);

    //  Product
    unique_ptr<Product<Number>> product;
    if (barrier > 0) product = make_unique<UOC<Number>>(strike, barrier, maturity, monitorFreq);
    else product = make_unique<European<Number>>(strike, maturity);

    //  RNG
    unique_ptr<RNG> rng;
    if (useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1, seed2);

    //  Simulate

    const auto simulResults = parallel
        ? mcParallelSimulAAD(*product, model, *rng, numPath)
        : mcSimulAAD(*product, model, *rng, numPath);

    //  Value

    results.value = accumulate(
        simulResults.payoffs.begin(),
        simulResults.payoffs.end(),
        0.0,
        [](const double acc, const vector<double>& v) { return acc + v[0]; }
            ) / numPath;

    //  Delta

    results.delta = simulResults.risks[0];

    //  Vegas

    results.vega.resize(model.spots().size(), model.times().size());
    copy(++simulResults.risks.begin(), simulResults.risks.end(), results.vega.begin());

    return results;
}

//  Returns a struct with price, delta and vega matrix of portfolio
inline auto europeansDupireAADRisk(
    //  model parameters
    const double            spot,
    const vector<double>&   spots,
    const vector<Time>&     times,
    //  spot major
    const matrix<double>&   vols,
    const double            maxDt,
    //  product parameters
    const vector<Time>&     maturities, 
    const vector<double>&   strikes,
    const vector<double>&   notionals,
    //  numerical parameters
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    //  optionals
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Results
    struct
    {
        double value;
        double delta;
        matrix<double> vega;
    } results;

    //  Build model, product and rng
    Dupire<Number> model(spot, spots, times, vols, maxDt);

    //  RNG
    unique_ptr<RNG> rng;
    if (useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1, seed2);

    //  Put in the right format
    map<Time, vector<double>> options, nots;
    for (size_t i = 0; i < maturities.size(); ++i)
    {
        options[maturities[i]].push_back(strikes[i]);
        nots[maturities[i]].push_back(notionals[i]);
    }

    //  Product
    unique_ptr<Product<Number>> product = make_unique<Europeans<Number>>(options);

    //  Aggregator lambda

    //  Flat notionals
    vector<double> flatNotionals;
    for (const auto& p : nots)
    {
        copy(p.second.begin(), p.second.end(), back_inserter(flatNotionals));
    }

    auto aggregator = [&flatNotionals] (const vector<Number>& payoffs)
    {
        return inner_product(flatNotionals.begin(), flatNotionals.end(), payoffs.begin(), Number(0.0));
    };

    //  Simulate

    const auto simulResults = parallel
        ? mcParallelSimulAAD(*product, model, *rng, numPath, aggregator)
        : mcSimulAAD(*product, model, *rng, numPath, aggregator);

    //  Value
    
    results.value = accumulate(
        simulResults.aggregated.begin(),
        simulResults.aggregated.end(),
        0.0) / numPath;

    //  Delta

    results.delta = simulResults.risks[0];

    //  Vegas

    results.vega.resize(model.spots().size(), model.times().size());
    copy(++simulResults.risks.begin(), simulResults.risks.end(), results.vega.begin());

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
    const char ivsType,
    const double spot,
    const double vol,
    const double jmpIntens = 0.0,
    const double jmpAverage = 0.0,
    const double jmpStd = 0.0)
{
    //  Create IVS
    auto ivs = ivsType == 'B' || ivsType == 'b'
        ? unique_ptr<IVS>(new BachelierIVS(spot, vol))
        : ivsType == 'S' || ivsType == 's'
        ? unique_ptr<IVS>(new BlackScholesIVS(spot, vol))
        : unique_ptr<IVS>(new MertonIVS(spot, vol, jmpIntens, jmpAverage, jmpStd));

    return dupireCalib(*ivs, inclSpots, maxDs, inclTimes, maxDt);
}

//  Returns value, delta, strikes, maturities 
//      and vega = derivatives to implied vols = superbucket
inline auto
    dupireSuperbucket(
    const double            spot,
    //  product parameters
    const double            strike,
    //  negative = european
    const double            barrier,
    const Time              maturity,
    const double            monitorFreq,
    //  numerical parameters
    const double            maxDtSimul,
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
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
    //  'B'achelier, Black'S'choles or 'M'erton
    const char              ivsType,
    const double            vol,
    const double            jmpIntens = 0.0,
    const double            jmpAverage = 0.0,
    const double            jmpStd = 0.0,
    //  optionals
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Results
    struct
    {
        double value;
        double delta;
        vector<double> strikes;
        vector<Time> mats;
        matrix<double> vega;
    } results;

    //  Start with a clean tape
    auto* tape = Number::tape;
    tape->rewind();

    //  Calibrate the model
    auto params = dupireCalib(
        inclSpots, 
        maxDs, 
        inclTimes, 
        maxDtVol, 
        ivsType, 
        spot, 
        vol, 
        jmpIntens, 
        jmpAverage, 
        jmpStd);
    const vector<double>& spots = params.spots;
    const vector<Time>& times = params.times;
    const matrix<double>& lvols = params.lVols;

    //  Find delta and microbucket
    auto mdlDerivs = uocDupireAADRisk(
        spot,
        spots,
        times,
        lvols,
        maxDtSimul,
        strike,
        barrier,
        maturity,
        monitorFreq,
        parallel,
        useSobol,
        numPath,
        seed1,
        seed2);
    results.value = mdlDerivs.value;
    results.delta = mdlDerivs.delta;
    const matrix<double>& microbucket = mdlDerivs.vega;

    //  Clear tape
    //  tape->rewind();
    tape->clear();

    //  Convert market inputs to numbers, put on tape
            
    //  Create IVS
    auto ivs = ivsType == 'B' || ivsType == 'b'
        ? unique_ptr<IVS>(new BachelierIVS(spot, vol))
        : ivsType == 'S' || ivsType == 's'
        ? unique_ptr<IVS>(new BlackScholesIVS(spot, vol))
        : unique_ptr<IVS>(new MertonIVS(spot, vol, jmpIntens, jmpAverage, jmpStd));
    
    //  Risk view --> that is the AAD input
    //  Note: that puts the view on tape
    RiskView<Number> riskView(strikes, mats);

    //  Calibrate again, in AAD mode, make tape
    auto nParams = dupireCalib(*ivs, inclSpots, maxDs, inclTimes, maxDtVol, riskView);
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
    Number::propagateAdjoints(tape->back(), tape->begin());

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

inline double uocBS(
    //  model parameters
    const double            spot,
    const double            vol,
    const bool              normal,
    const double            rate,
    const double            div,
    //  product parameters
    const double            strike,
    //  negative = european
    const double            barrier,
    const Time              maturity,
    const double            monitorFreq,
    //  numerical parameters
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Build model, product and rng
    BlackScholes<double> model(spot, vol, normal, rate, div);

    //  Product
    unique_ptr<Product<double>> product;
    if (barrier > 0) product = make_unique<UOC<double>>(strike, barrier, maturity, monitorFreq);
    else product = make_unique<European<double>>(strike, maturity);

    //  RNG
    unique_ptr<RNG> rng;
    if (useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1, seed2);

    //  Simulate
    const auto resultMat = parallel
        ? mcParallelSimul(*product, model, *rng, numPath)
        : mcSimul(*product, model, *rng, numPath);

    //  Compute averages among paths
    double result = accumulate(resultMat.begin(), resultMat.end(), 0.0,
        [](const double acc, const vector<double>& v) { return acc + v[0]; }
    ) / numPath;

    return result;
}

//  Returns a struct with price, delta and vega matrix
inline auto uocBSAADRisk(
    //  model parameters
    const double            spot,
    const double            vol,   
    const bool              normal,
    const double            rate,
    const double            div,
    //  product parameters
    const double            strike,
    //  negative = european
    const double            barrier,
    const Time              maturity,
    const double            monitorFreq,
    //  numerical parameters
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    //  optionals
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Results
    struct
    {
        double value;
        double delta;
        double vega;
        double rho;
        double ddiv;
    } results;

    //  Build model, product and rng
    BlackScholes<Number> model(spot, vol, normal, rate, div);

    //  Product
    unique_ptr<Product<Number>> product;
    if (barrier > 0) product = make_unique<UOC<Number>>(strike, barrier, maturity, monitorFreq);
    else product = make_unique<European<Number>>(strike, maturity);

    unique_ptr<RNG> rng;
    if (useSobol) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1, seed2);

    //  Simulate

    const auto simulResults = parallel
        ? mcParallelSimulAAD(*product, model, *rng, numPath)
        : mcSimulAAD(*product, model, *rng, numPath);

    //  Value

    results.value = accumulate(
        simulResults.payoffs.begin(),
        simulResults.payoffs.end(),
        0.0,
        [](const double acc, const vector<double>& v) { return acc + v[0]; }
    ) / numPath;

    //  Delta

    results.delta = simulResults.risks[0];
    results.vega = simulResults.risks[1];
    results.rho = simulResults.risks[2];
    results.ddiv = simulResults.risks[3];

    return results;
}
