#pragma once

#include "mcBase.h"
#include "mcMdl.h"
#include "mcPrd.h"
#include "mrg32k3a.h"
#include "sobol.h"
#include "memory.h"
#include <numeric>

inline double uocDupire(
    //  model parameters
    const double            spot,
    const vector<double>    spots,
    const vector<Time>      times,
    const matrix<double>    vols,   //  spot major
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

//  Returns price, delta and matrix of vegas
inline tuple<double, double, matrix<double>> uocDupireBumpRisk(
    //  model parameters
    const double            spot,
    const vector<double>    spots,
    const vector<Time>      times,
    const matrix<double>    vols,   //  spot major
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
    //  base value
    const double v0 = uocDupire(
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
    const double delta = 1.e+08 * (v1 - v0);

    //  Vega
    matrix<double> newVols = vols;
    matrix<double> vega(vols.rows(), vols.cols());
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
            vega[i][j] = 1.e+08 * (v1 - v0);
            newVols[i][j] -= 1.e-08;
        }
    }

    return make_tuple(v0, delta, vega);
}

//  Returns price, delta and matrix of vegas
inline tuple<double, double, matrix<double>> uocDupireAADRisk(
    //  model parameters
    const double            spot,
    const vector<double>    spots,
    const vector<Time>      times,
    const matrix<double>    vols,   //  spot major
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
    //  Build model, product and rng
    Dupire<Number> model(spot, spots, times, vols, maxDt);
    UOC<Number> product(strike, barrier, maturity, monitorFreq);
    unique_ptr<RNG> rng = useSobol 
        ? unique_ptr<RNG>(new Sobol) 
        : unique_ptr<RNG>(new mrg32k3a(seed1, seed2));

    //  Simulate

    const auto results = parallel
        ? mcParallelSimulAAD(product, model, *rng, numPath, antithetic)
        : mcSimulAAD(product, model, *rng, numPath, antithetic);

    //  Value

    const double value = accumulate(results.first.begin(), 
        results.first.end(), 
        0.0)
            / numPath;

    //  Risks

    //  Downcast the model, we know it is a Dupire
    Dupire<Number>& resMdl = *dynamic_cast<Dupire<Number>*>(results.second.get());

    double delta = resMdl.spot().adjoint();
    matrix<double> vega(resMdl.vols().rows(), resMdl.vols().cols());
    transform(resMdl.vols().begin(), resMdl.vols().end(), vega.begin(),
        [](const Number& vol)
    {
        return vol.adjoint();
    });
    //  Normalize
    delta /= numPath;
    for (auto& risk : vega) risk /= numPath;

    //  Clear the tape
    Number::tape->clear();

    //  Return value
    return make_tuple(value, delta, vega);
}

//  Returns spots, times and local vols
inline tuple<vector<double>, vector<Time>, matrix<double>>
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
//      and derivatives to calls = superbucket
inline tuple<double, double, vector<double>, vector<Time>, matrix<double>>
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
    const vector<double>& spots = get<0>(params);
    const vector<Time>& times = get<1>(params);
    const matrix<double>& lvols = get<2>(params);

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
    const double value = get<0>(mdlDerivs);
    const double delta = get<1>(mdlDerivs);
    const matrix<double>& microbucket = get<2>(mdlDerivs);

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
    matrix<Number>& nLvols = get<2>(nParams);

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
    matrix<double> superbucket(riskView.rows(), riskView.cols());
    transform(riskView.begin(), riskView.end(), superbucket.begin(), 
        [](const Number& n)
    {
        return n.adjoint();
    });

    //  Clear tape
    tape->clear();

    //  Return results
    return make_tuple(value, delta, strikes, mats, superbucket);
}