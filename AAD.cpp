
#include "AADnumber.h"
#include <algorithm>

const struct Number::LeafType Number::leaf;
const struct Number::UnaryType Number::unary;
const struct Number::BinaryType Number::binary;

Tape globalTape;
thread_local Tape* Number::tape = &globalTape;

//  Check-pointing

template<class T>
inline T H(const vector<T>& Y)
{
    return T();
}

template<class T>
inline vector<T> G(const vector<T>& X)
{
    return vector<T>();
}

//  Implements check-pointing
//  Takes input X
//  Computes and returns F(X) = H[G(X)] and its derivatives
inline pair<double, vector<double>> checkPoint(const vector<double>& X)
{
    //  Start with a clean tape
    auto* tape = Number::tape;
    tape->clear();

    //  1
    //  Compute Y
    vector<double> Y = G(X);
    //  Convert to numbers
    vector<Number> iY(Y.size());
    convertCollection(Y.begin(), Y.end(), iY.begin());
    //  Note that also puts iY on tape

    //  2
    Number z = H(iY);

    //  3
    //  Propagate
    z.propagateToStart();
    //  Store derivatives
    vector<double> dhdy(iY.size());
    transform(iY.begin(), iY.end(), dhdy.begin(),
        [](const Number& y) {return y.adjoint(); });

    //  4
    tape->clear();

    //  5
    vector<Number> iX(X.size());
    convertCollection(X.begin(), X.end(), iX.begin());

    //  6
    vector<Number> oY = G(iX);

    //  7
    for (size_t i = 0; i < oY.size(); ++i)
    {
        oY[i].adjoint() = dhdy[i];
    }

    //  8
    Number::propagateAdjoints(tape->back(), tape->begin());

    //  9
    pair<double, vector<double>> results;
    results.first = z.value();
    results.second.resize(iX.size());
    transform(iX.begin(), iX.end(), results.second.begin(),
        [](const Number& x) {return x.adjoint(); });

    //  10
    tape->clear();

    return results;
}