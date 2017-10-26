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
    const vector<double>    times,
    const matrix<double>    vols,   //  spot major
    const double            maxDt,
    //  product parameters
    const double            strike,
    const double            barrier,
    const double            maturity,
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

inline double uocDupireBumpRisk(
    //  model parameters
    const double            spot,
    const vector<double>    spots,
    const vector<double>    times,
    const matrix<double>    vols,   //  spot major
    const double            maxDt,
    //  product parameters
    const double            strike,
    const double            barrier,
    const double            maturity,
    const double            monitorFreq,
    //  numerical parameters
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    const bool              antithetic,
    //  risk outputs
    double&                 delta,
    matrix<double>&         vega,
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
    delta = 1.e+08 * (v1 - v0);

    //  Vega
    matrix<double> newVols = vols;
    vega.resize(vols.rows(), vols.cols());
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

    return v0;
}

inline double uocDupireAADRisk(
    //  model parameters
    const double            spot,
    const vector<double>    spots,
    const vector<double>    times,
    const matrix<double>    vols,   //  spot major
    const double            maxDt,
    //  product parameters
    const double            strike,
    const double            barrier,
    const double            maturity,
    const double            monitorFreq,
    //  numerical parameters
    const bool              parallel,
    const bool              useSobol,
    const int               numPath,
    const bool              antithetic,
    //  risk outputs
    double&                 delta,
    matrix<double>&         vega,
    //  optionals
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Build model, product and rng
    Dupire<Number> model(spot, spots, times, vols, maxDt);
    UOC<Number> product(strike, barrier, maturity, monitorFreq);
    unique_ptr<RNG> rng = useSobol ? unique_ptr<RNG>(new Sobol) : unique_ptr<RNG>(new mrg32k3a(seed1, seed2));

    //  Simulate

    const auto results = parallel
        ? mcParallelSimulAAD(product, model, *rng, numPath, antithetic)
        : mcSimulAAD(product, model, *rng, numPath, antithetic);

    //  Value

    const double value = accumulate(results.first.begin(), results.first.end(), 0.0)
        / numPath;

    //  Risks

    //  Downcast the model, we know it is a Dupire
    Dupire<Number>& resMdl = *dynamic_cast<Dupire<Number>*>(results.second.get());

    delta = resMdl.spot().adjoint();
    vega.resize(resMdl.vols().rows(), resMdl.vols().cols());
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
    return value;
}