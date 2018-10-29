#pragma once

#include "mcBase.h"
#include "interp.h"

inline double toyDupireBarrierMc(
    //  Spot
    const double            S0,
    //  Local volatility
    const vector<double>    spots,
    const vector<Time>      times,
    const matrix<double>    vols,
    //  Product parameters
    const double            maturity,
    const double            strike,
    const double            barrier,
    //  Number of paths and time steps
    const int               Np,
    const int               Nt,
    //  Random number generator
    RNG&                    random)
{
    //  Initialize
    double result = 0;
    random.init(Nt);
    vector<double> gaussianIncrements(Nt);
    const double dt = maturity / Nt, sdt = sqrt(dt);

    //  Loop over paths
    for (size_t i = 0; i < Np; ++i)
    {
        //  Generate Nt Gaussian Numbers
        random.nextG(gaussianIncrements);

        //  Step by step
        double spot = S0, time = 0;
        bool alive = true;
        for (size_t j = 0; j < Nt; ++j)
        {
            //  Interpolate volatility
            const double vol = interp2D(spots, times, vols, spot, time);
            time += dt;

            //  Simulate return
            spot *= exp(-0.5 * vol * vol * dt + vol * sdt * gaussianIncrements[j]);

            //  Monitor barrier
            if (spot > barrier)
            {
                alive = false;
                break;
            }
        }

        //  Payoff
        if (alive && spot > strike) result += spot - strike;

    }   //  paths

    return result / Np;
}