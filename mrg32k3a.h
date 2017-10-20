#pragma once

#include "mcBase.h"
#include "gaussians.h"

//  Positive modulus required by mrg32ka
inline long long mod(const long long a, const long long b)
{
    const long long c = a % b;
    return c > 0 ? c : c + b;
}

class mrg32k3a : public RNG
{
    //  Dimension
    size_t          myDim;

    //  Pre-allocated placeholder for random numbers
    vector<double>  myNumbers;

    //  State
    long long       myXn, myXn1, myXn2, myYn, myYn1, myYn2;

    //  Constants
    static const long long m1 = 4294967087;
    static const long long m2 = 4294944443;
    static const long long a11 = 0;
    static const long long a12 = 1403580;
    static const long long a13 = -810728;
    static const long long a21 = 527612;
    static const long long a22 = 0;
    static const long long a23 = -1370589;

    //  Produce next number and update state
    double nextNumber()
    {
        //  Update states
        const long long x = mod(a11 * myXn + a12 * myXn1 + a13 * myXn2, m1);
        const long long y = mod(a21 * myYn + a22 * myYn1 + a23 * myYn2, m2);
        myXn2 = myXn1;
        myXn1 = myXn;
        myXn = x;
        myYn2 = myYn1;
        myYn1 = myYn;
        myYn = x;

        //  Compute uniform
        const double u = double(mod(x - y, m1)) / m1;
        //  Return Gaussian
        return invNormalCdf(u);
    }

public:

    //  Constructor with seed
    mrg32k3a(const unsigned a = 12345, const unsigned b = 12346) :
        myXn(a), myXn1(a), myXn2(a), myYn(b), myYn1(b), myYn2(b)
    {}

    //  Virtual copy constructor
    unique_ptr<RNG> clone() const override
    {
        return unique_ptr<RNG>(new mrg32k3a(*this));
    }

    //  Initializer 
    void init(const size_t simDim) override      
    {
        myDim = simDim;
        myNumbers.resize(simDim);
    }

    //  We cache the vector of numbers and return it by const reference
    const vector<double>& nextG() override
    {
        for (size_t i = 0; i < myDim; ++i) myNumbers[i] = nextNumber();
        return myNumbers;
    }
};