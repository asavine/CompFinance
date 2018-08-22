
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

//  Implementation of the mrg32k3a RNG,
//  See chapters 5 and 6

#include "mcBase.h"
#include "gaussians.h"

class mrg32k3a : public RNG
{

	//	Seed
	const double	myA, myB;
	
	//  Dimension
    size_t			myDim;

	//  State
    double			myXn, myXn1, myXn2, myYn, myYn1, myYn2;

	//	Antithetic
	bool			myAnti;
	//	false: generate new, true: negate cached
	vector<double>	myCachedUniforms;
	vector<double>	myCachedGaussians;

    //  Constants
    static constexpr  double	m1 = 4294967087;
	static constexpr  double	m2 = 4294944443;
	static constexpr  double	a12 = 1403580;
	static constexpr  double	a13 = 810728;
	static constexpr  double	a21 = 527612;
	static constexpr  double	a23 = 1370589;
	//	We divide the final uniform 
	//		by m1 + 1 so we never hit 1
	static constexpr  double	m1p1 = 4294967088;

    //  Produce next number and update state
    double nextNumber()
    {
        //  Update X
		//	Recursion
		double x = a12 * myXn1 - a13 * myXn2;
		//	Modulus
		x -= long(x / m1) * m1;
		if (x < 0) x += m1;
		//	Update
		myXn2 = myXn1;
		myXn1 = myXn;
		myXn = x;

		//	Same for Y
		double y = a21 * myYn - a23 * myYn2;
		y -= long(y / m2) * m2;
		if (y < 0) y += m2;
		myYn2 = myYn1;
		myYn1 = myYn;
		myYn = y;

        //  Uniform 
        const double u = x > y 
			? (x - y) / m1p1 
			: (x - y + m1) / m1p1;
        return u;
    }

public:

    //  Constructor with seed
    mrg32k3a(const unsigned a = 12345, const unsigned b = 12346) :
        myA(a), myB(b)
    {
        reset();
    }

	//	Reset state to 0 (seed)
    void reset()
    {
		//	Reset state
        myXn = myXn1 = myXn2 = myA;
        myYn = myYn1 = myYn2 = myB;
		
		//	Anti = false: generate next
		myAnti = false;
    }

    //  Virtual copy constructor
    unique_ptr<RNG> clone() const override
    {
        return make_unique<mrg32k3a>(*this);
    }

    //  Initializer 
    void init(const size_t simDim) override      
    {
        myDim = simDim;
		myCachedUniforms.resize(myDim);
		myCachedGaussians.resize(myDim);
    }

	void nextU(vector<double>& uVec) override
	{
		if (myAnti)
		{
			//	Do not generate, negate cached
			transform(
				myCachedUniforms.begin(),
				myCachedUniforms.end(),
				uVec.begin(),
				[](const double d) { return 1.0 - d; });
			
			//	Generate next
			myAnti = false;
		}
		else
		{
			//	Generate and cache
			generate(
				myCachedUniforms.begin(), 
				myCachedUniforms.end(), 
				[this]() { return nextNumber(); });
			
			//	Copy
			copy(
				myCachedUniforms.begin(),
				myCachedUniforms.end(),
				uVec.begin());
			
			//	Do not generate next
			myAnti = true;
		}
	}

    void nextG(vector<double>& gaussVec) override
    {
		if (myAnti)
		{
			//	Do not generate, negate cached
			//	Note: we reuse the Gaussian numbers,
			//		we save not only generation, but also
			//		Gaussian transaformation
			transform(
				myCachedGaussians.begin(),
				myCachedGaussians.end(),
				gaussVec.begin(),
				[](const double n) { return -n; });

			//	Generate next
			myAnti = false;
		}
		else
		{
			//	Generate and cache
			generate(
				myCachedGaussians.begin(),
				myCachedGaussians.end(),
				[this]() { return invNormalCdf(nextNumber()); });

			//	Copy
			copy(
				myCachedGaussians.begin(),
				myCachedGaussians.end(),
				gaussVec.begin());

			//	Do not generate next
			myAnti = true;
		}
	}

	//	Skip ahead logic
	//	See chapter 7
	//	To avoid overflow, we nest mods in innermost results
	//		and use 64bit unsigned long long for storage

	//  Skip ahead
	void skipTo(const unsigned b) override
	{
		//	First reset to 0
		reset();

		//	How many numbers to skip
		unsigned skipnums = b * myDim;
		bool odd = false;

		//	Antithetic: skip only half
		if (skipnums & 1)
		{
			//	Odd
			odd = true;
			skipnums = (skipnums - 1) / 2;
		}
		else
		{
			//	Even
			skipnums /= 2;
		}

		//	Skip state
		skipNumbers(skipnums);

		//	If odd, pre-generate for antithetic
		if (odd)
		{
			myAnti = true;

			//	Uniforms
			generate(
				myCachedUniforms.begin(),
				myCachedUniforms.end(),
				[this]() { return nextNumber(); });

			//	Gaussians
			generate(
				myCachedGaussians.begin(),
				myCachedGaussians.end(),
				[this]() { return invNormalCdf(nextNumber()); });
		}
		else
		{
			myAnti = false;
		}
	}

private:

	//  Matrix product with modulus
	static void mPrd(
		const unsigned long long	lhs[3][3],
		const unsigned long long	rhs[3][3],
		const unsigned long long	mod,
		unsigned long long			result[3][3])
	{
		//	Result go to temp, in case result points to lhs or rhs
		unsigned long long temp[3][3];

		for (size_t j = 0; j<3; j++)
		{
			for (size_t k = 0; k<3; k++)
			{
				unsigned long long s = 0;
				for (size_t l = 0; l<3; l++)
				{
					//	Apply modulus to innermost product
					unsigned long long tmpNum = lhs[j][l] * rhs[l][k];
					//	Apply mod
					tmpNum %= mod;
					//	Result
					s += tmpNum;
					//	Reapply mod
					s %= mod;
				}
				//  Store result in temp 
				temp[j][k] = s;
			}
		}

		//	Now product is done, copy temp to result
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				result[j][k] = temp[j][k];
			}
		}
	}

	//  Matrix by vector, exact same logic
	//	Except we don't implement temp,
	//		we never point result to lhs or rhs
	static void vPrd(
		const unsigned long long	lhs[3][3],
		const unsigned long long	rhs[3],
		const unsigned long long	mod,
		unsigned long long			result[3])
	{
		for (size_t j = 0; j<3; j++)
		{
			unsigned long long s = 0;
			for (size_t l = 0; l<3; l++)
			{
				unsigned long long tmpNum = lhs[j][l] * rhs[l];
				tmpNum %= mod;
				s += tmpNum;
				s %= mod;
			}
			result[j] = s;
		}
	}

	void skipNumbers(const unsigned b) 
    {
        if ( b <= 0) return;
        unsigned skip = b;

		static constexpr unsigned long long
			m1l = unsigned long long(m1);
		static constexpr unsigned long long
			m2l = unsigned long long(m2);

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
                { 
					0, 
					unsigned long long (a12) , 
					unsigned long long (m1 - a13) 
					//	m1 - a13 instead of -a13
					//	so results are always positive
					//	and we can use unsigned long longs
					//	after modulus, we get the same results
				},
                { 1, 0, 0 },
                { 0, 1, 0 }
        },
            Bi[3][3] = {        //  B0 = B
                { 
					unsigned long long (a21), 
					0 , 
					unsigned long long (m2 - a23) 
					//	same logic: m2 - a32
				},
                { 1, 0, 0 },
                { 0, 1, 0 }
        };

        while (skip>0) 
        {
            if (skip & 1)   //  i.e. ai == 1
            {
                //  accumulate Ab and Bb
                mPrd(Ab, Ai, m1l, Ab);
                mPrd(Bb, Bi, m2l, Bb);
            }

            //  Recursion on Ai and Bi 
            mPrd(Ai, Ai, m1l, Ai);
            mPrd(Bi, Bi, m2l, Bi);

            skip >>= 1;
        }

        //  Final result
		unsigned long long X0[3] =
        {
			unsigned long long (myXn),
			unsigned long long (myXn1),
			unsigned long long (myXn2)
        },
            Y0[3] =
        {
			unsigned long long (myYn),
			unsigned long long (myYn1),
			unsigned long long (myYn2)
        },
            temp[3];
        
		//	From initial to final state
        vPrd(Ab, X0, m1l, temp);
        
		//	Convert back to doubles
		myXn = double(temp[0]);
        myXn1 = double(temp[1]);
        myXn2 = double(temp[2]);

		//	Same for Y
        vPrd(Bb, Y0, m2l, temp);
        myYn = double(temp[0]);
        myYn1 = double(temp[1]);
        myYn2 = double(temp[2]);
    }
};