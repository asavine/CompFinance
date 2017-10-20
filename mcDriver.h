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
    const bool              useSobol,
    const int               numPath,
    const int               seed1 = 12345,
    const int               seed2 = 12346)
{
    //  Build model, product and rng
    Dupire<double> model(spot, spots, times, vols, maxDt);
    UOC<double> product(strike, barrier, maturity, monitorFreq);
    unique_ptr<RNG> rng = useSobol ? unique_ptr<RNG>(new Sobol) : unique_ptr<RNG>(new mrg32k3a(seed1, seed2));

    //  Simulate
    auto results = mcSimul(product, model, *rng, numPath);

    //  Return average
    return accumulate(results.begin(), results.end(), 0.0) / results.size();
}