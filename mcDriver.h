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

//  Returns price, delta and matrix of vegas
inline tuple<double, double, matrix<double>> uocDupireCheckPointedRisk(
    //  market parameters
    const double            spot,
    const vector<double>    strikes,
    const vector<Time>      mats,
    const matrix<double>    calls,   //  strike major
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
    //  Start with a clean tape
    auto* tape = Number::tape;
    tape->rewind();

    //  Calibrate the model
    auto params = dupireCalib(spot, strikes, mats, calls);
    const vector<double>& spots = get<0>(params);
    const vector<Time>& times = get<1>(params);
    const matrix<double>& lvols = get<2>(params);

    //  Find delta and microbucket
    auto mdlDerivs = uocDupireAADRisk(
        spot,
        spots,
        times,
        lvols,
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
    const double value = get<0>(mdlDerivs);
    const double modelDelta = get<1>(mdlDerivs);
    const matrix<double>& microbucket = get<2>(mdlDerivs);

    //  Clear tape
    tape->rewind();

    //  Convert market inputs to numbers, put on tape
    Number nSpot(spot);
    matrix<Number> nCalls(calls.rows(), calls.cols());
    convertCollection(calls.begin(), calls.end(), nCalls.begin());

    //  Calibrate again, in AAD mode, make tape
    auto nParams = dupireCalib(nSpot, strikes, mats, nCalls);
    matrix<Number>& nLvols = get<2>(nParams);

    //  Seed tape
    nSpot.adjoint() = modelDelta;
    for (size_t i = 0; i < microbucket.rows(); ++i)
    {
        for (size_t j = 0; j < microbucket.cols(); ++j)
        {
            nLvols[i][j].adjoint() = microbucket[i][j];
        }
    }
    
    //  Propagate
    Number::propagateAdjoints(tape->back(), tape->begin());

    //  Pack results
    double marketDelta = nSpot.adjoint();
    matrix<double> superbucket(calls.rows(), calls.cols());
    for (size_t i = 0; i < calls.rows(); ++i)
    {
        for (size_t j = 0; j < calls.cols(); ++j)
        {
            superbucket[i][j] = nCalls[i][j].adjoint();
        }
    }

    //  Clear tape
    tape->clear();

    //  Return results
    return make_tuple(value, marketDelta, superbucket);
}    

