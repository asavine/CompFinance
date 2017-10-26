#pragma once

#include "mcBase.h"
#include "gaussians.h"

class mrg32k3a : public RNG
{
    //  Dimension
    size_t                      myDim;

    //  State
    unsigned long long          myXn, myXn1, myXn2, myYn, myYn1, myYn2;

    //  Constants
    static const unsigned long long m1 = 4294967087;
    static const unsigned long long m2 = 4294944443;
    static const unsigned long long a12 = 1403580;
    static const unsigned long long a13 = 810728;
    static const unsigned long long a21 = 527612;
    static const unsigned long long a23 = 1370589;
    static const unsigned long long m1Minusa13 = 4294156359;
    static const unsigned long long m2Minusa23 = 4293573854;

    //  Produce next number and update state
    double nextNumber()
    {
        //  Update states
        const unsigned long long x = (a12 * myXn1) % m1 + (m1Minusa13 * myXn2) % m1;
        const unsigned long long y = (a21 * myYn) % m2 + (m2Minusa23 * myYn2) % m2;
        myXn2 = myXn1;
        myXn1 = myXn;
        myXn = x > m1 ? x - m1 : x;
        myYn2 = myYn1;
        myYn1 = myYn;
        myYn = y > m2 ? y - m2 : y;

        //  Compute uniform
        const double u = x > y ? double((x - y) % m1) / m1 : double((x + m1 - y) % m1) / m1;
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
    }

    //  We cache the vector of numbers and return it by const reference
    void nextG(vector<double>& gaussVec) override
    {
        for (size_t i = 0; i < myDim; ++i) gaussVec[i] = nextNumber();
    }

    //  Access dimension
    size_t simDim() const override
    {
        return myDim;
    }

    //  Matrix product with modulus
    static void mPrd(const unsigned long long lhs[3][3],
        const unsigned long long rhs[3][3],
        const unsigned long long mod,
        unsigned long long result[3][3])
    {
        unsigned long long temp[3][3];

        for (size_t j = 0; j<3; j++)
        {
            for (size_t k = 0; k<3; k++)
            {
                unsigned long long s = 0;
                for (size_t l = 0; l<3; l++)
                {
                    s = (s + lhs[j][l] * rhs[l][k]) % mod;
                }
                //  Temp in case result == lhs or rhs
                temp[j][k] = s;
            }
        }
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                result[j][k] = temp[j][k];
            }
        }
    }

    //  Matrix by vector
    static void vPrd(const unsigned long long lhs[3][3],
        const unsigned long long rhs[3],
        const unsigned long long mod,
        unsigned long long result[3])
    {
        for (size_t j = 0; j<3; j++)
        {
            unsigned long long s = 0;
            for (size_t l = 0; l<3; l++)
            {
                s = (s + lhs[j][l] * rhs[l]) % mod;
            }
            result[j] = s;
        }
    }

    //  Skip ahead
    void skipTo(const long b) override
    {
        if ( b <= 0) return;
        long skip = b * myDim;

        unsigned long long Ab[3][3] = {
            { 1, 0 ,0 },        
            { 0, 1, 0 },
            { 0, 0, 1 }
        },
            Bb[3][3] = {
                { 1, 0 ,0 },    
                { 0, 1, 0 },
                { 0, 0, 1 }
        },
            Ai[3][3] = {        //  A0 = A
                { 0, a12 , m1Minusa13 },
                { 1, 0, 0 },
                { 0, 1, 0 }
        },
            Bi[3][3] = {        //  B0 = B
                { a21, 0 , m2Minusa23 },
                { 1, 0, 0 },
                { 0, 1, 0 }
        };

        while (skip>0) 
        {
            if (skip & 1)   //  i.e. ai == 1
            {
                //  cumulate Ab and Bb
                mPrd(Ab, Ai, m1, Ab);
                mPrd(Bb, Bi, m2, Bb);
            }

            //  Recusion on Ai and Bi 
            mPrd(Ai, Ai, m1, Ai);
            mPrd(Bi, Bi, m2, Bi);

            skip >>= 1;
        }

        //  Final result
        unsigned long long X0[3] =
        {
            myXn,
            myXn1,
            myXn2
        },
            Y0[3] =
        {
            myYn,
            myYn1,
            myYn2
        },
            temp[3];
        
        vPrd(Ab, X0, m1, temp);
        myXn = temp[0];
        myXn1 = temp[1];
        myXn2 = temp[2];

        vPrd(Bb, Y0, m2, temp);
        myYn = temp[0];
        myYn1 = temp[1];
        myYn2 = temp[2];
    }

};