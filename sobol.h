// this file is an edited version of
// the file provided with
// © 2002 "Monte Carlo Methods in Finance"
// on the CD accompanying the book
// "Monte Carlo Methods in Finance" by Peter Jaeckel.
//
// the original copyright is:
// ===========================================================================
// Copyright (C) 2002 "Monte Carlo Methods in Finance". All rights reserved.
//
// Permission to use, copy, modify, and distribute this software is freely
// granted, provided that this notice is preserved.
// ===========================================================================

//  Edited by Antoine Savine in 2018: essentially added a skip ahead (skipTo)

#pragma once
#pragma warning(disable : 4018)

#include "mcBase.h"
#include "gaussians.h"
#include "matrix.h"
#include "primitivepolynomials.h"
#include "initializers.h"

//	content:	generates sobol sequences up to dimension 1110
//	
//				Joe and Kou directional numbers guarantee Sobol property A to be satisfied
//
//				Code is based on QuantLib 0.3.6 and routines from Jaeckel
//	
//				
//
//				References:
//
//				Jaeckel (2002): Monte-Carlo Methods in Finance, Wiley
//
//				Joe & Kuo (2003): Remark on Algorithm 659: Implementing Sobol's 
//				Quasirandom Sequence Generator. ACM Transactions on Mathematical Software 29, 1, 43-57
//

#define ONEOVER2POW32 2.32830643653870E-10

class Sobol : public RNG
{
    //  Dimension
    size_t      myDim;

    //	All initializers and polynomials
    const unsigned long* const*	myInitializers;
    const long* const*			myPrimitivePolynomials;

    //  Used polynomials and degrees
    vector<long>		        myPpmt;
    vector<unsigned>            myDegree;

    //	Direction numbers
    matrix<unsigned long>	    myDirectionIntegers;

    //	Cached Uniforms
    vector<unsigned long>	    myIntegerSequence;  //  The current numbers, as integers

    //  The current point number
    unsigned long				mySimCount;

public:

    //  Virtual copy constructor
    unique_ptr<RNG> clone() const override
    {
        return unique_ptr<RNG>(new Sobol(*this));
    }

    //  Initializer 
    void init(const size_t simDim) override
    {
        myDim = simDim;

        //	Allocate direction numbers
        myDirectionIntegers.resize(myDim, 32);

        //	Set pointers
        myInitializers = get_jk_initializers();

        //  The array jk_PrimitivePolynomials is identical to Jaeckel's original in all 
        //  except the ordering of the polynomials of degree 4 to 8. These appear as a 
        //  renamed copy of the original towards the end of primitivepolynomials.cpp
        myPrimitivePolynomials = jk_PrimitivePolynomials;

        //	Initialize coeffcients and degree of the k'th primitive polynomial
        myPpmt.resize(myDim);
        myDegree.resize(myDim);
        myPpmt[0] = 0;
        myDegree[0] = 0;
        int k, index, currentDegree;
        for (k = 1, index = 0, currentDegree = 1; k<myDim; ++k, ++index)
        {
            myPpmt[k] = myPrimitivePolynomials[currentDegree - 1][index];
            if (myPpmt[k] == -1)
            {
                ++currentDegree;
                index = 0;
                myPpmt[k] = myPrimitivePolynomials[currentDegree - 1][index];
            }
            myDegree[k] = currentDegree;
        }

        //	Initialize recurrence for direction numbers for each dimension 

        //	Dimension 0 is degenerate
        int j;
        for (j = 0; j<32; ++j)
        {
            myDirectionIntegers[0][j] = (1UL << (32 - j - 1));
        }

        //	Dimension 1 to myDim jk numbers are used for direction integers 
        for (k = 1; k<(int)myDim; ++k)
        {
            j = 0;

            //	end of line marked with 0UL
            while (myInitializers[k - 1][j] != 0UL)
            {
                myDirectionIntegers[k][j] = myInitializers[k - 1][j];
                myDirectionIntegers[k][j] <<= (32 - j - 1);
                ++j;
            }
        }

        //	Run recurrence for direction numbers

        for (k = 1; k<(int)myDim; ++k)
        {
            unsigned int gk = myDegree[k];
            for (int l = gk; l<32; ++l)
            {
                //	see eq. 8.19 and comments in Jackel 
                unsigned long n = (myDirectionIntegers[k][l - gk] >> gk);

                //	a(k,j) := ppmt(k)>>(gk-j-1) are the coefficients 
                //      of the monomials in ppmt(k)
                for (long j = 1; j<(long)gk; ++j)
                {
                    //	use the xor operator
                    if ((myPpmt[k] >> (gk - j - 1)) & 1UL)
                    {
                        n ^= myDirectionIntegers[k][l - j];
                    }
                }

                //	The last xor
                n ^= myDirectionIntegers[k][l - gk];
                myDirectionIntegers[k][l] = n;
            }
        }

        //	Resize
        myIntegerSequence.resize(myDim);

        //  Generate the first point
        unsigned i;
        for (i = 0; i<myDim; ++i)
        {
            myIntegerSequence[i] = myDirectionIntegers[i][0];
        }

        //	Set the count
        mySimCount = 1;
    }
	
	//	next point
    void nextG(vector<double>& gaussVec) override
    {
        //	For the n'th draw use the gray code
        unsigned long n = mySimCount;
        unsigned int  j = 0;
        while (n & 1)
        {
            n >>= 1;
            ++j;
        }

        //	XOR the appropriate direction number into each component of the integer sequence
        for (int i = 0; i<myDim; ++i)
        {
            myIntegerSequence[i] ^= myDirectionIntegers[i][j];
        }
        if (!gaussVec.empty())
        {
            for (int i = 0; i < myDim; ++i)
            {
                gaussVec[i] = invNormalCdf(ONEOVER2POW32 * myIntegerSequence[i]);
            }
        }

        //	Update count
        ++mySimCount;
    }

    //  Access dimension
    size_t simDim() const override
    {
        return myDim;
    }

    //  Skip ahead (from 0 to b)
    void skipTo(const long b) override
    {
        //	Check skip
        if (b <= 0) return;

        //	Reset Sobol to entry 0 
        //  (not 1, hence must reset even though reset has already been called in init)
        for (size_t i = 0; i<myDim; ++i) myIntegerSequence[i] = 0;

        //	The actual Sobol skipping algo
        unsigned long long im = static_cast<unsigned long long>(b);
        unsigned int	   two_i = 1, two_i_plus_one = 2;

        unsigned i = 0;
        while (two_i <= im)
        {
            if (((im + two_i) / two_i_plus_one) & 1)
            {
                for (unsigned int k = 0; k<myDim; ++k)
                {
                    myIntegerSequence[k] ^= myDirectionIntegers[k][i];
                }
            }

            two_i <<= 1;
            two_i_plus_one <<= 1;
            ++i;
        }

        //	End of skipping algo

        //	Update next entry
        mySimCount = unsigned long(b);
        nextG(vector<double>());
    }
};



