//  primitive polynomials
#include "primitivepolynomials.h"


//	--------------- kSobol class ------------

//	bits on machine
const int kSobol::myBits    = 8*sizeof(unsigned long);

//	2^-myBits, written as 0.5/2^(myBits-1)
const double kSobol::myNorm = 0.5/(1UL<<(kSobol::myBits-1));

//	max table dimension
unsigned
kSobol::maxTableDim() const
{
	return sizeof(jk_initializers)/sizeof(unsigned long *) + 1;
}

//	max sobol dimension
unsigned 
kSobol::maxSobolDim() const
{
	return PPMT_MAX_DIM;
}

//	default constructor
kSobol::kSobol()
:	myTableDim(0),
	mySobolDim(0),
	mySimDim(0),
	myInitializers(0),
	myPrimitivePolynomials(0),
	mySimCount(0)
{}

//	added AS, init and skip ahead all in one
kSobol::kSobol( 
	unsigned			totalDim,
	unsigned			sobolDim,
	unsigned			jkDim,
	unsigned			randomize,
	unsigned long long	skip,
	long				seed11,
	long				seed12,
	long				seed21,
	long				seed22)
:	myTableDim(0),
	mySobolDim(0),
	mySimDim(0),
	myInitializers(0),
	myPrimitivePolynomials(0),
	mySimCount(0)
{
	kValArray<long> seeds(4);
	seeds(0) = seed11;
	seeds(1) = seed12;
	seeds(2) = seed21;
	seeds(3) = seed22;

	init(jkDim,sobolDim,totalDim,randomize,seeds);

	skipAhead(skip);
}

//	init
void
kSobol::init(
	unsigned				tableDim, 
	unsigned				sobolDim, 
	unsigned				simDim, 
	unsigned				randomize, 
	const kValArray<long>&	seeds)
{
	//	set sobol dim
	sobolDim = kInlines<unsigned>::min(sobolDim,maxSobolDim());
	sobolDim = kInlines<unsigned>::min(sobolDim,simDim);

	//	set table dim
	tableDim = kInlines<unsigned>::min(tableDim,maxTableDim());
	tableDim = kInlines<unsigned>::min(tableDim,sobolDim);

	//	set dimensions
	myTableDim = tableDim;
	mySobolDim = sobolDim;
	mySimDim   = simDim;

	LOG( MoreInfo) << "Sobol" << endl;
	LOG( MoreInfo) << "Total dim = " << mySimDim << endl;
	LOG( MoreInfo) << "Sobol dim = " << mySobolDim << endl;
	LOG( MoreInfo) << "JK dim = " << myTableDim << endl;
	LOG( MoreInfo) << "Randomize = " << randomize << endl;

	//	resize
	myDirectionIntegers.resize(mySobolDim,myBits);

	//	init random generator for initialization of sequences above tableDim, we flip seeds
//	myRandomGenerator.init(seeds(1),seeds(0));
	myRandomGenerator.init(seeds(0),seeds(1));

	//	set to jk data
	myInitializers		   = jk_initializers;
	myPrimitivePolynomials = jk_PrimitivePolynomials;

	//	init the sobols
	if(mySobolDim)
	{
		//	initialize coeffcients and degree of the k'th primitive polynomial
        myPpmt.resize(mySobolDim);
        myDegree.resize(mySobolDim);
		myPpmt(0)   = 0;
		myDegree(0) = 0;
		int k, index, currentDegree;
		for(k=1,index=0,currentDegree=1;k<(int)mySobolDim;++k,++index)
		{
			myPpmt(k) = myPrimitivePolynomials[currentDegree-1][index];
			if(myPpmt(k)==-1)
			{
				++currentDegree;
				index   = 0;
				myPpmt(k) = myPrimitivePolynomials[currentDegree-1][index];
			}
			myDegree(k) = currentDegree;
		}

		//	initialize myBits direction for each dimension and store in myDirectionIntegers

		//	dimension 0 is degenerate
		int j;
		for(j=0;j<myBits;++j)
		{
			myDirectionIntegers(0,j) = (1UL<<(myBits-j-1));
		}

		//	dimension 1 to myTableDim jk numbers are used for direction integers 
		for(k=1;k<(int)myTableDim;++k)
		{
			j = 0;

			//	end of line marked with 0UL
			while(myInitializers[k-1][j]!=0UL)
			{
				myDirectionIntegers(k,j) = myInitializers[k-1][j];
				myDirectionIntegers(k,j) <<= (myBits - j - 1);
				++j;
			}
		}

		//	random numbers for higher dimensions
		double u;
		for(k=myTableDim;k<(long)mySobolDim;++k)
		{
			for(long l=1;l<=(long)myDegree(k);++l)
			{
				//	do until we hit an odd direction
				do
				{
					//	draw uniform
					u = myRandomGenerator.uniform();

					//	the direction integer has at most the l rightmost bits non-zero
					myDirectionIntegers(k,l-1) = (unsigned long)(u*(1UL<<l));
				}
				while(!(myDirectionIntegers(k,l-1) & 1UL));

				//	shifting myBits-1 to the left, we guarantee that the l'th leftmost but us set and 
				//	only the the first l leftmost bits can be non-zero
				myDirectionIntegers(k,l-1) <<= (myBits-l);
			}
		}

		//	computation of myDirectionIntegers(k,l) for l>=degree(k) by recurrence relation
		for(k=1;k<(int)mySobolDim;++k)
		{
			unsigned int gk = myDegree(k);
			for(int l=gk;l<myBits;++l)
			{
				//	see eq. 8.19 and comments in Jackel 
				unsigned long n = (myDirectionIntegers(k,l-gk)>>gk);

				//	a(k,j) := ppmt(k)>>(gk-j-1) are the coefficients of the monomials in ppmt(k)
                for(long j=1;j<(long)gk;++j) 
                {
                    //	use the xor operator
                    if((myPpmt(k)>>(gk-j-1)) & 1UL)
                    {
                        n ^= myDirectionIntegers(k,l-j);
                    }
                }

                //	the last xor
                n ^= myDirectionIntegers(k,l-gk);
                myDirectionIntegers(k,l) = n;
			}
		}
	}

	//	resize
	myIntegerSequence.resize(mySobolDim);
	myRandomizers.resize(mySimDim);
	myUniforms.resize(mySimDim);

	//	re-init random generator for the randomizers, we flip seeds back
//	myRandomGenerator.init(seeds(0),seeds(1));
	if(mySimDim)
	{
		myRandomizers = 0.0;
		unsigned h, i;
		for(h=0;h<randomize;++h)
		{
			for(i=0;i<mySimDim;++i)
			{
				myRandomizers(i) = myRandomGenerator.uniform();		
			}
		}
	}

	//	re-init random generator again for filling above sobolDim, we use seeds(2..3)
	myRandomGenerator.init(seeds(2),seeds(3));

	//	reset
	reset();
}

//	reset
void
kSobol::reset()
{
	unsigned i;
	for(i=0;i<mySobolDim;++i)
	{
		myIntegerSequence(i) = myDirectionIntegers(i,0);
		myUniforms(i)        = myIntegerSequence(i)*myNorm + myRandomizers(i);
		if(myUniforms(i)>1.0) myUniforms(i) -= 1.0;
	}

	myRandomGenerator.reset();
	for(i=mySobolDim;i<mySimDim;++i)
	{
		myUniforms(i) = myRandomGenerator.uniform() + myRandomizers(i);
		if(myUniforms(i)>1.0) myUniforms(i) -= 1.0;
	}

	//	set sim count
	mySimCount = 1;
}

//	skip ahead
void 
kSobol::skipAhead(
	unsigned long long	skip) //	Argument = number of entries to skip
{	
	
	//	Check skip
	if (!skip) return;

	unsigned int i;

	//	Reset Sobol to entry 0 (not 1, hence must reset even though reset has already been called in init)
	for(i=0;i<mySobolDim;++i) myIntegerSequence(i) = 0;
	myRandomGenerator.reset();	//	Also reset the random generator
		
	//	The actual Sobol skipping algo
	unsigned long long im = skip;
	unsigned int	   two_i = 1, two_i_plus_one = 2;

	i = 0;
	while (two_i <= im) 
	{
		if ( ( (im + two_i) / two_i_plus_one ) & 1 ) 
		{
			for (unsigned int k=0; k<mySobolDim; ++k)
			{
				myIntegerSequence( k) ^= myDirectionIntegers( k, i);
			}
		}

		two_i <<= 1;
		two_i_plus_one <<= 1;
		++i;
	}

	//	End of skipping algo

	//	The random generator must also skip
	if(mySimDim>mySobolDim) myRandomGenerator.skipAhead( skip * (mySimDim - mySobolDim));

	//	Update next entry
	mySimCount = unsigned long(skip);
	next();
}

//	next 
void
kSobol::next()
{
	//	for the n'th draw use the gray code
	unsigned long n = mySimCount;
	unsigned int  j = 0;
	while(n&1)
	{
		n >>= 1;
		++j;
	}

	//	xor the appropriate direction number into each component of the integer sequence
	unsigned int i;
	for(i=0;i<mySobolDim;++i)
	{
		myIntegerSequence(i) ^= myDirectionIntegers(i,j);
		myUniforms(i)        =  myIntegerSequence(i)*myNorm + myRandomizers(i);
		if(myUniforms(i)>1.0) myUniforms(i) -= 1.0;
	}

	//	fill with randoms
	for(i=mySobolDim;i<mySimDim;++i)
	{
		myUniforms(i) = myRandomGenerator.uniform() + myRandomizers(i);
		if(myUniforms(i)>1.0) myUniforms(i) -= 1.0;
	}

	//	next
	++mySimCount;
}







