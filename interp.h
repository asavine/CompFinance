#pragma once

#include <algorithm>
using namespace std;

//  Utility for interpolation
//  Interpolates the vector y against knots x in value x0 
//  Interpolation is linear, extrapolation is flat
template <class ITX, class ITY, class T>
inline auto linterp(
    //	sorted non xs 
    ITX                         xBegin,
    ITX                         xEnd,
    //	corresponding ys
    ITY                         yBegin,
    ITY                         yEnd,
    //	interpolate for point x0
    const T&					x0)
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
    return y1 + (y2 - y1) / (x2 - x1) * (x0 - x1);
}