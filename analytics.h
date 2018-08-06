#pragma once

#include "gaussians.h"

//  Black and Scholes formula, templated

//  Price
template<class T, class U, class V, class W>
inline T blackScholes(
    const U spot,
    const V strike,
    const T vol,
    const W mat)
{
    const auto std = vol * sqrt(mat);
    if (std <= EPS) return max<T>(T(0.0), T(spot - strike));
    const auto d2 = log(spot / strike) / std - 0.5 * std;
    const auto d1 = d2 + std;
    return spot * normalCdf(d1) - strike * normalCdf(d2);
}

//  Implied vol, untemplated
inline double blackScholesIvol(
    const double spot,
    const double strike,
    const double prem,
    const double mat)
{
    if (prem <= max(0.0, spot - strike) + EPS) return 0.0;

    double p, pu, pl;
    double u = 0.5;
    while (blackScholes(spot, strike, u, mat) < prem) u *= 2;
    double l = 0.05;
    while (blackScholes(spot, strike, l, mat) > prem) l /= 2;
    pu = blackScholes(spot, strike, u, mat);
    blackScholes(spot, strike, l, mat);

    while (u - l > 1.e-12)
    {
        const double m = 0.5 * (u + l);
        p = blackScholes(spot, strike, m, mat);
        if (p > prem)
        {
            u = m;
            pu = p;
        }
        else
        {
            l = m;
            pl = p;
        }
    }

    return l + (prem - pl) / (pu - pl) * (u - l);
}

//  Merton, templated

template<class T, class U, class V, class W, class X>
inline T merton(
    const U spot,
    const V strike,
    const T vol,
    const W mat,
    const X intens,
    const X meanJmp,
    const X stdJmp)
{
    const auto varJmp = stdJmp * stdJmp;
    const auto mv2 = meanJmp + 0.5 * varJmp;
    const auto comp = intens * (exp(mv2) - 1);
    const auto var = vol * vol;
    const auto intensT = intens * mat;

    unsigned fact = 1;
    X iT = 1.0;
    const size_t cut = 10;
    T result = 0.0;
    for (size_t n = 0; n < cut; ++n)
    {
        const auto s = spot*exp(n*mv2 - comp*mat);
        const auto v = sqrt(var + n * varJmp / mat);
        const auto prob = exp(-intensT) * iT / fact;
        result += prob * blackScholes(s, strike, v, mat);
        fact *= n + 1;
        iT *= intensT;
    }

    return result;
}

//	Up and out call in Black-Scholes, untemplated

inline double BlackScholesKO(
    const double spot, 
    const double rate,
    const double div, 
    const double strike, 
    const double barrier, 
    const double mat, 
    const double vol) 
{
    const double std = vol * sqrt(mat);
    const double fwdFact = exp((rate - div)*mat);
    const double fwd = spot * fwdFact;
    const double disc = exp(-rate * mat);
    const double v = rate - div - 0.5 * vol * vol;
    const double d2 = (log(spot / barrier) + v * mat) / std;
    const double d2prime = (log(barrier / spot) + v * mat) / std;

    const double bar = blackScholes(fwd, strike, vol, mat)
        - blackScholes(fwd, barrier, vol, mat)
        - (barrier - strike) * normalCdf(d2)
        - pow(barrier / spot, 2 * v / vol / vol) *
        (
            blackScholes(fwdFact * barrier* barrier / spot, strike, vol, mat)
            - blackScholes(fwdFact * barrier* barrier / spot, barrier, vol, mat)
            - (barrier - strike) * normalCdf(d2prime) 
            );

    return disc * bar;
}
