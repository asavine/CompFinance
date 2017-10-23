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

    //  Skip ahead
    void skipTo(const long b) override
    {
        if ( b <= 0) return;
        long skip = b * myDim;

        unsigned long long mat1[3][3] = {
            { 1, 0 ,0 },
            { 0, 1, 0 },
            { 0, 0, 1 }
        },
            mat2[3][3] = {
                { 1, 0 ,0 },
                { 0, 1, 0 },
                { 0, 0, 1 }
        },
            t1[3][3] = {
                { 0, 0 ,0 },
                { 0, 0, 0 },
                { 0, 0, 0 }
        },
            t2[3][3] = {
                { 0, 0 ,0 },
                { 0, 0, 0 },
                { 0, 0, 0 }
        },
            temp1[3][3], temp2[3][3];

        //	Transition matrices
        t1[0][1] = a12;
        t1[0][2] = m1Minusa13;
        t1[1][0] = 1;
        t1[2][1] = 1;

        t2[0][0] = a21;
        t2[0][2] = m2Minusa23;
        t2[1][0] = 1;
        t2[2][1] = 1;

        while (skip>0) 
        {
            //  Update matrix product
            if (skip & 1) 
            {
                for (int j = 0; j<3; j++) 
                {
                    for (int k = 0; k<3; k++) 
                    {
                        unsigned long long s1 = 0, s2 = 0;
                        for (int l = 0; l<3; l++) 
                        {
                            s1 = (s1 + t1[j][l] * mat1[l][k]) % m1;
                            s2 = (s2 + t2[j][l] * mat2[l][k]) % m2;
                        }
                        temp1[j][k] = s1;
                        temp2[j][k] = s2;
                    }
                }
                for (int j = 0; j<3; j++) for (int k = 0; k<3; k++) 
                {
                    mat1[j][k] = temp1[j][k];
                    mat2[j][k] = temp2[j][k];
                }
            }

            //  Compute the next square
            for (int j = 0; j<3; j++) 
            {
                for (int k = 0; k<3; k++) 
                {
                    unsigned long long s1 = 0, s2 = 0;
                    for (int l = 0; l<3; l++) 
                    {
                        s1 = (s1 + t1[j][l] * t1[l][k]) % m1;
                        s2 = (s2 + t2[j][l] * t2[l][k]) % m2;
                    }
                    temp1[j][k] = s1;
                    temp2[j][k] = s2;
                }
            }
            for (int j = 0; j<3; j++) for (int k = 0; k<3; k++) 
            {
                t1[j][k] = temp1[j][k];
                t2[j][k] = temp2[j][k];
            }

            skip >>= 1;
        }

        const unsigned long long tempS10 = 
            (mat1[0][0] * myXn) % m1 + 
            (mat1[0][1] * myXn1) % m1 + 
            (mat1[0][2] * myXn2) % m1;
        const unsigned long long tempS11 = 
            (mat1[1][0] * myXn) % m1 + 
            (mat1[1][1] * myXn1) % m1 + 
            (mat1[1][2] * myXn2) % m1;
        const unsigned long long tempS12 = 
            (mat1[2][0] * myXn) % m1 + 
            (mat1[2][1] * myXn1) % m1 + 
            (mat1[2][2] * myXn2) % m1;
        myXn = tempS10 % m1;
        myXn1 = tempS11 % m1;
        myXn2 = tempS12 % m1;

        const unsigned long long tempS20 = 
            (mat2[0][0] * myYn) % m2 + 
            (mat2[0][1] * myYn1) % m2 + 
            (mat2[0][2] * myYn2) % m2;
        const unsigned long long tempS21 = 
            (mat2[1][0] * myYn) % m2 + 
            (mat2[1][1] * myYn1) % m2 + 
            (mat2[1][2] * myYn2) % m2;
        const unsigned long long tempS22 = 
            (mat2[2][0] * myYn) % m2 + 
            (mat2[2][1] * myYn1) % m2 + 
            (mat2[2][2] * myYn2) % m2;
        myYn = tempS20 % m2;
        myYn1 = tempS21 % m2;
        myYn2 = tempS22 % m2;
    }

};