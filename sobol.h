
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
#pragma warning(disable : 4018)

//  Implementation of Sobol's sequence,
//  See chapters 5 and 6

#include "mcBase.h"
#include "gaussians.h"

#define ONEOVER2POW32 2.3283064365387E-10

const unsigned * const * getjkDir();

class Sobol : public RNG
{
    //  Dimension
    size_t                      myDim;

    //  State Y
    vector<unsigned>	        myState;  

    //  Current index in the sequence
    unsigned                    myIndex;

    //  The direction numbers listed in sobol.cpp
    //  Note jkDir[i][dim] gives the i-th (0 to 31) 
    //      direction number of dimension dim
    const unsigned * const *    jkDir;

public:

    //  Virtual copy constructor
    unique_ptr<RNG> clone() const override
    {
        return make_unique<Sobol>(*this);
    }

    //  Initializer 
    void init(const size_t simDim) override
    {
        //  Set pointer on direction numbers 
        jkDir = getjkDir();

        //  Dimension
        myDim = simDim;
        myState.resize(myDim);

        //  Reset to 0
        reset();
    }

    void reset()
    {
        //  Set state to 0
        memset(myState.data(), 0, myDim * sizeof(unsigned));
        //  Set index to 0
        myIndex = 0;
    }
	
	//	Next point
	void next() 
	{
		//	Gray code, find position j 
		//		of rightmost zero bit of current index n
		unsigned n = myIndex, j = 0;
		while (n & 1)
		{
			n >>= 1;
			++j;
		}

        //  Direction numbers
        const unsigned* dirNums = jkDir[j];

		//	XOR the appropriate direction number 
		//		into each component of the integer sequence
        for (int i = 0; i<myDim; ++i)
		{
			myState[i] ^= dirNums[i];
		}

		//	Update count
		++myIndex;
	}

	void nextU(vector<double>& uVec) override
	{
		next();
		transform(myState.begin(), myState.end(), uVec.begin(),
			[](const unsigned long i) 
				{return ONEOVER2POW32 * i; });
	}

	void nextG(vector<double>& gaussVec) override
    {
		next();
		transform(myState.begin(), myState.end(), gaussVec.begin(),
			[](const unsigned long i) 
				{return invNormalCdf(ONEOVER2POW32 * i); });
    }

    //  Skip ahead (from 0 to b)
    void skipTo(const unsigned b) override
    {
        //	Check skip
        if (!b) return;

        //	Reset Sobol to 0 
        reset();

        //	The actual Sobol skipping algo
        unsigned im = b;
        unsigned two_i = 1, two_i_plus_one = 2;

        unsigned i = 0;
        while (two_i <= im)
        {
            if (((im + two_i) / two_i_plus_one) & 1)
            {
                for (unsigned k = 0; k<myDim; ++k)
                {
                    myState[k] ^= jkDir[i][k];
                }
            }

            two_i <<= 1;
            two_i_plus_one <<= 1;
            ++i;
        }

        //	End of skipping algo

        //	Update next entry
        myIndex = unsigned(b);
    }
};
