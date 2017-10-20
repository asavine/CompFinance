#pragma once

#include "mcBase.h"
#include "gaussians.h"
#include "matrix.h"
#include "primitivepolynomials.h"
#include "initializers.h"

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
    vector<double>			    myGaussians;        //  The Gaussian vector, cached

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
        myGaussians.resize(myDim);

        //  Generate the first point
        unsigned i;
        for (i = 0; i<myDim; ++i)
        {
            myIntegerSequence[i] = myDirectionIntegers[i][0];
            myGaussians[i] = invNormalCdf(myIntegerSequence[i]*ONEOVER2POW32);
        }

        //	Set the count
        mySimCount = 1;
    }
	
	//	next point
    const vector<double>& nextG() override
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
        unsigned int i;
        for (i = 0; i<myDim; ++i)
        {
            myIntegerSequence[i] ^= myDirectionIntegers[i][j];
            myGaussians[i] = invNormalCdf(myIntegerSequence[i]*ONEOVER2POW32);
        }

        //	Update count
        ++mySimCount;

        //  Return
        return myGaussians;
    }
};



