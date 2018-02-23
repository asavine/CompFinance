#pragma once

#include "mcBase.h"
#include "mcMdl.h"
#include "mcPrd.h"
#include "mrg32k3a.h"
#include "sobol.h"
#include "memory.h"
#include <numeric>
#include <fstream>
using namespace std;

inline double uocDupire(
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
    const bool              antithetic,
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Build model, product and rng
    Dupire<double> model(spot, spots, times, vols, maxDt);
    UOC<double> product(strike, barrier, maturity, monitorFreq);
    unique_ptr<RNG> rng = useSobol ? unique_ptr<RNG>(new Sobol) : unique_ptr<RNG>(new mrg32k3a(seed1, seed2));

    //  Simulate
    auto results = parallel
        ? mcParallelSimul(product, model, *rng, numPath, antithetic)
        : mcSimul(product, model, *rng, numPath, antithetic);

    //  Return average
    return accumulate(results.begin(), results.end(), 0.0) 
        / results.size();
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
    const bool              antithetic,
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
        antithetic,
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
        antithetic,
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
                antithetic,
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
    const bool              antithetic,
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
    UOC<Number> product(strike, barrier, maturity, monitorFreq);
    unique_ptr<RNG> rng = useSobol 
        ? unique_ptr<RNG>(new Sobol) 
        : unique_ptr<RNG>(new mrg32k3a(seed1, seed2));

    //  Simulate

    const auto simulResults = parallel
        ? mcParallelSimulAAD(product, model, *rng, numPath, antithetic)
        : mcSimulAAD(product, model, *rng, numPath, antithetic);

    //  Value

    results.value = accumulate(
        simulResults.payoffs.begin(),
        simulResults.payoffs.end(),
        0.0)
            / numPath;

    //  Risks

    //  Downcast the model, we know it is a Dupire
    Dupire<Number>& resMdl = *static_cast<Dupire<Number>*>(simulResults.model.get());

    results.delta = resMdl.spot().adjoint();
    results.vega.resize(resMdl.vols().rows(), resMdl.vols().cols());
    transform(resMdl.vols().begin(), resMdl.vols().end(), results.vega.begin(),
        [numPath](const Number& vol)
    {
        return vol.adjoint() / numPath;
    });
    //  Normalize
    results.delta /= numPath;

    //  Clear the tape
    Number::tape->clear();

    //  Return value
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
    const double            barrier,
    const Time              maturity,
    const double            monitorFreq,
    //  numerical parameters
    const double            maxDtSimul,
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    const bool              antithetic,
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
        antithetic,
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