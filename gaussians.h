#pragma once

#include <math.h>
#include <vector>
#include <algorithm>
using namespace std;

//  Normal density
inline double normalDens(const double x)
{
    return x<-10.0 || 10.0<x ? 0.0 : exp(-0.5*x*x) / 2.506628274631;
}

//	Normal CDF (N in Black-Scholes)
inline double normalCdf( const double x)
{
	//	checks
	if (x<-10.0) return 0.0;
	if (x>10.0) return 1.0;
	if (x<0.0) return 1.0 - normalCdf(-x);

	//  calc pol 

	//	constants
	static const double p = 0.2316419;
	static const double b1 = 0.319381530;
	static const double b2 = -0.356563782;
	static const double b3 = 1.781477937;
	static const double b4 = -1.821255978;
	static const double b5 = 1.330274429;

	//	transform
	double t = 1.0 / (1.0 + p*x);

	//	finally pol
	double pol = t*(b1 + t*(b2 + t*(b3 + t*(b4 + t*b5))));

	//	calc pdf
	double pdf = x<-10.0 || 10.0<x ? 0.0 : exp(-0.5*x*x) / 2.506628274631; // sqrt (2 * pi())

	//	return cdf
	return 1.0 - pdf * pol;
}

//	Inverse CDF (for generation of Gaussians out of Uniforms)
inline double invNormalCdf( const double p)
{
	//	to ensure symmetry
    const bool sup = p > 0.5;
    const double up = sup ? 1.0 - p : p;

	//	constants
	static const double a0 = 2.50662823884;
	static const double a1 = -18.61500062529;
	static const double a2 = 41.39119773534;
	static const double a3 = -25.44106049637;

	static const double b0 = -8.47351093090;
	static const double b1 = 23.08336743743;
	static const double b2 = -21.06224101826;
	static const double b3 = 3.13082909833;

	static const double c0 = 0.3374754822726147;
	static const double c1 = 0.9761690190917186;
	static const double c2 = 0.1607979714918209;
	static const double c3 = 0.0276438810333863;
	static const double c4 = 0.0038405729373609;
	static const double c5 = 0.0003951896511919;
	static const double c6 = 0.0000321767881768;
	static const double c7 = 0.0000002888167364;
	static const double c8 = 0.0000003960315187;

	//	send x negative in all cases
	double x = up - 0.5;
	double r;

	//	polymonomial approx
	if (fabs(x)<0.42)
	{
		r = x*x;
		r = x*(((a3*r + a2)*r + a1)*r + a0) / ((((b3*r + b2)*r + b1)*r + b0)*r + 1.0);
		return sup ? -r: r;
	}

	//	log log approx
	r = up;
	r = log(-log(r));
	r = c0 + r*(c1 + r*(c2 + r*(c3 + r*(c4 + r*(c5 + r*(c6 + r*(c7 + r*c8)))))));

	//	flip sign 

	//	done
	return sup? r: -r;
}

//  turn a uniform vector into a gaussian vector
inline void u2g(const vector<double>& u, vector<double>& g)
{
    transform(u.begin(), u.end(), g.begin(), invNormalCdf);
}

//  Bachelier's formula and implied volatility

inline double bachelier(
    const double spot,
    const double strike,
    const double vol,
    const double mat)
{
    const double std = vol * sqrt(mat);
    if (std <= 0) return spot - strike;
    const double d = (spot - strike) / std;
    return (spot - strike) * normalCdf(d) + std * normalDens(d);
}

inline double bachelierIvol(
    const double spot,
    const double strike,
    const double prem,
    const double mat)
    {
    if (prem <= spot - strike) return 0.0;
    double u = 0.5 * spot;
    while (bachelier(spot, strike, u, mat) < prem) u *= 2;
    double l = 0.05 * spot;
    while (bachelier(spot, strike, l, mat) > prem) u /= 2;

    while (u - l > 1.e-06)
    {
        const double m = 0.5 * (u + l);
        if (bachelier(spot, strike, m, mat) > prem)
        {
            u = m;
        }
        else
        {
            l = m;
        }
    }

    return 0.5 * (u + l);
}

inline double blackScholes(
    const double spot,
    const double strike,
    const double vol,
    const double mat)
{
    const double std = vol * sqrt(mat);
    if (std <= 0) return spot - strike;
    const double d2 = log(spot/strike) / std - 0.5 * std;
    const double d1 = d2 + std;
    return spot * normalCdf(d1) - strike * normalCdf(d2);
}

inline double merton(
    const double spot,
    const double strike,
    const double vol,
    const double mat,
    const double intens,
    const double meanJmp,
    const double stdJmp)
{
    const double varJmp = stdJmp * stdJmp;
    const double mv2 = meanJmp + 0.5 * varJmp;
    const double comp = intens * (exp(mv2) - 1);
    const double var = vol * vol;
    const double intensT = intens * mat;

    unsigned fact = 1;
    double iT = 1;
    const size_t cut = 10;
    double result = 0.0;
    for (size_t n = 0; n < cut; ++n)
    {
        const double s = spot*exp(n*mv2 - comp*mat);
        const double v = sqrt(var + n * varJmp / mat);
        const double prob = exp(-intensT) * iT / fact;
        result += prob * blackScholes(s, strike, v, mat);
        fact *= n + 1;
        iT *= intensT;
    }

    return result;
}


