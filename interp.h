#pragma once

#include <algorithm>
using namespace std;

//  Utility for interpolation
//  Interpolates the vector y against knots x in value x0 
//  Interpolation is linear, extrapolation is flat
template <class ITX, class ITY, class T>
inline auto interp(
    //	sorted on xs 
    ITX                         xBegin,
    ITX                         xEnd,
    //	corresponding ys
    ITY                         yBegin,
    ITY                         yEnd,
    //	interpolate for point x0
    const T&					x0,
    //  smooth?
    const bool                  smoothStep = false)
{
    //	STL binary search, returns iterator on 1st no less than x0
    //  upper_boung guarantees logarithmic search
    auto it = upper_bound(xBegin, xEnd, x0);

    //  Extrapolation?
    if (it == xEnd) return *(yEnd - 1);
    if (it == xBegin) return *yBegin;

    //  Interpolation
    size_t n = distance(xBegin, it) - 1;
    auto x1 = xBegin[n];
    auto y1 = yBegin[n];
    if (x0 == x1) return y1;
    auto x2 = xBegin[n + 1];
    auto y2 = yBegin[n + 1];

    auto t = (x0 - x1) / (x2 - x1);

    return smoothStep
        ? y1 + (y2 - y1) * t * t * (3.0 - 2 * t)
        : y1 + (y2 - y1) * t;
}

//  2D
template <class T, class U, class V, class W, class X>
inline auto interp2D(
    //	sorted on xs 
    const vector<T>&            x,
    //	sorted on ys
    const vector<U>&            y,
    //  zs in a matrix
    const matrix<V>&            z,
    //	interpolate for point (x0,y0)
    const W&					x0,
    const X&					y0,
    //  smooth?
    const bool                  smoothStep = false)
{
    const size_t n = x.size();
    const size_t m = y.size();

    //	STL binary search, returns iterator on 1st no less than x0
    //  upper_boung guarantees logarithmic search
    auto it = upper_bound(x.begin(), x.end(), x0);
    const size_t n2 = distance(x.begin(), it);

    //  Extrapolation in x?
    if (n2 == n)
        return interp(y.begin(), y.end(), z[n2 - 1], z[n2 - 1] + m, y0, smoothStep);
    if (n2 == 0)
        return interp(y.begin(), y.end(), z[0], z[0] + m, y0, smoothStep);

    //  Interpolation in x
    const size_t n1 = n2 - 1;
    auto x1 = x[n1];
    auto x2 = x[n2];
    auto z1 = interp(y.begin(), y.end(), z[n1], z[n1] + m, y0, smoothStep);
    auto z2 = interp(y.begin(), y.end(), z[n2], z[n2] + m, y0, smoothStep);

    //  Smooth step
    auto t = (x0 - x1) / (x2 - x1);
    return smoothStep
        ? z1 + (z2 - z1) * t * t * (3.0 - 2 * t)
        : z1 + (z2 - z1) * t;
}