#pragma once

#include <math.h>
#include <vector>
#include <algorithm>
using namespace std;

#define EPS 1.0e-08

//  Gaussian functions, including Bachelier, BS and Merton
//  Templated as required

//  Normal density
template<class T>
inline T normalDens(const T x)
{
    return x<-10.0 || 10.0<x ? T(0.0) : exp(-0.5*x*x) / 2.506628274631;
}

//	Normal CDF (N in Black-Scholes)
template<class T>
inline T normalCdf( const T x)
{
	//	checks
	if (x < -10.0) return T(0.0);
	if (x > 10.0) return T(1.0);
    if (x < 0.0) return 1.0 - normalCdf(-x);

	//  calc pol 

	//	constants
	static const double p = 0.2316419;
	static const double b1 = 0.319381530;
	static const double b2 = -0.356563782;
	static const double b3 = 1.781477937;
	static const double b4 = -1.821255978;
	static const double b5 = 1.330274429;

	//	transform
	const T t = 1.0 / (1.0 + p*x);

	//	finally pol
    const T pol = t*(b1 + t*(b2 + t*(b3 + t*(b4 + t*b5))));

	//	calc pdf
    const T pdf = x<-10.0 || 10.0<x ? T(0.0) : exp(-0.5*x*x) / 2.506628274631; // sqrt (2 * pi())

	//	return cdf
	return 1.0 - pdf * pol;
}

//	Inverse CDF (for generation of Gaussians out of Uniforms)
//  Untemplated
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

//  Bachelier's formula and implied volatility
template<class T, class U, class V, class W>
inline T bachelier(
    const U spot,
    const V strike,
    const T vol,
    const W mat)
{
    const auto std = vol * sqrt(mat);
    if (std < EPS) return max( T(0.0), T(spot - strike));
    const auto d = (spot - strike) / std;
    return (spot - strike) * normalCdf(d) + std * normalDens(d);
}

//  Vega
inline double bachelierVega(
    const double spot,
    const double strike,
    const double vol,
    const double mat)
{
    const double std = vol * sqrt(mat);
    if (std < EPS) return 0.0;
    const double d = (spot - strike) / std;
    return sqrt(mat) * normalDens(d);
}

//  BS
template<class T, class U, class V, class W>
inline T blackScholes(
    const U spot,
    const V strike,
    const T vol,
    const W mat)
{
    const auto std = vol * sqrt(mat);
    if (std <= EPS) return max(T(0.0), T(spot - strike));
    const auto d2 = log(spot/strike) / std - 0.5 * std;
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

//  Vega
inline double blackScholesVega(
    const double spot,
    const double strike,
    const double vol,
    const double mat)
{
    const double smat = sqrt(mat), std = vol * smat;
    if (std < EPS) return 0.0;
    const double d2 = log(spot / strike) / std - 0.5 * std;
    return strike * smat * normalDens(d2);
}

//  Merton
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