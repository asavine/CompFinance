
/*
Written by Antoine Savine in 2018

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications
is freely granted to any person or company who purchased a copy of the book

Modern Computational Finance: AAD and Parallel Simulations
Antoine Savine
Wiley, 2018

As long as this comment is preserved at the top of the file
*/

#pragma once

#include <algorithm>
using namespace std;

//  Interpolation utility 

//  Interpolates the vector y against knots x in value x0 
//  Interpolation is linear or smooth, extrapolation is flat
template <bool smoothStep=false, class ITX, class ITY, class T>
inline auto interp(
    //	sorted on xs 
    ITX                         xBegin,
    ITX                         xEnd,
    //	corresponding ys
    ITY                         yBegin,
    ITY                         yEnd,
    //	interpolate for point x0
    const T&					x0)
    ->remove_reference_t<decltype(*yBegin)>
{
    //	STL binary search, returns iterator on 1st no less than x0
    //  upper_bound guarantees logarithmic search
    auto it = upper_bound(xBegin, xEnd, x0);

    //  Extrapolation?
    if (it == xEnd) return *(yEnd - 1);
    if (it == xBegin) return *yBegin;

    //  Interpolation
    size_t n = distance(xBegin, it) - 1;
    auto x1 = xBegin[n];
    auto y1 = yBegin[n];
    auto x2 = xBegin[n + 1];
    auto y2 = yBegin[n + 1];

    auto t = (x0 - x1) / (x2 - x1);

    if constexpr (smoothStep)
    {
        return y1 + (y2 - y1) * t * t * (3.0 - 2 * t);
    }
    else
    {
        return y1 + (y2 - y1) * t;
    }
}

//  2D
template <bool smoothStep=false, class T, class U, class V, class W, class X>
inline V interp2D(
    //	sorted on xs 
    const vector<T>&            x,
    //	sorted on ys
    const vector<U>&            y,
    //  zs in a matrix
    const matrix<V>&            z,
    //	interpolate for point (x0,y0)
    const W&					x0,
    const X&					y0)
{
    const size_t n = x.size();
    const size_t m = y.size();

    //	STL binary search, returns iterator on 1st no less than x0
    //  upper_boung guarantees logarithmic search
    auto it = upper_bound(x.begin(), x.end(), x0);
    const size_t n2 = distance(x.begin(), it);

    //  Extrapolation in x?
    if (n2 == n)
        return interp<smoothStep>(y.begin(), y.end(), z[n2 - 1], z[n2 - 1] + m, y0);
    if (n2 == 0)
        return interp<smoothStep>(y.begin(), y.end(), z[0], z[0] + m, y0);

    //  Interpolation in x
    const size_t n1 = n2 - 1;
    auto x1 = x[n1];
    auto x2 = x[n2];
    auto z1 = interp<smoothStep>(y.begin(), y.end(), z[n1], z[n1] + m, y0);
    auto z2 = interp<smoothStep>(y.begin(), y.end(), z[n2], z[n2] + m, y0);

    //  Smooth step
    auto t = (x0 - x1) / (x2 - x1);
    if constexpr (smoothStep)
    {
        return z1 + (z2 - z1) * t * t * (3.0 - 2 * t);
    }
    else
    {
        return z1 + (z2 - z1) * t;;
    }
}
