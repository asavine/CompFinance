
/*

- finalize and recode a, b <-> ATM, skew, rem equations in comments, will ref text
check

- check initialized model, a, b and choldc
check

- upgrade engine to multiple, named, assets, update all models and products, check adequation of model and product in engine
check

- reg tests

- Autocalls

- Test, test, test

- Risks

- Rem public data members, comments, test functions etc.

*/



/*
Written by Antoine Savine in 2019

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications
is granted by explicit application, antoine@asavine.com

This comment must be preserved at all times
*/

//  Multiple correlated assets, all following a displaced Black-Scholes (DLM) dynamics, subject to discrete proportional dividends
//	DLMs are calibrated to short term ATM and skew, see https://www.researchgate.net/publication/335146739_Displaced_Lognormal_Models

//	On a dividend date: S+ = (1-div) * S-		--	by conv, dividend is applied end of day, that is, after everything else
//	Prop divs don't affect BS implied surface

//	Correl matrix is increased with a factor lambda
//	used_correl = (1-lambda) * correl + lambda * 1

#pragma once

#include <set>

#include "matrix.h"
#include "choldc.h"

template <class T>
class MultiDisplaced : public Model<T>
{

	//	For debugging only, goes private after that
public:
	//	

    enum Dynamics
    {
        Normal = 0,
        Surnormal = 1,
        Subnormal = 2
    };

    //  Model parameters

	//	Number of assets and reference prices
	size_t				myNumAssets;

    //  Today's market
    
	//  Constant rate 
    T                   myRate;

	//	Spots
	vector<T>			mySpots;

	//	Divs
	
	//	Dates where at least one asset pays a dividend
	vector<Time>		myDivDates;

	//	Matrix (numDivDates, numAsets) of dividends by date and asset
	matrix<T>			myDivs;
    
	//	Volatility, by asset

    //  ATM vol
    vector<T>			myAtms;
	//	skew
	vector<T>			mySkews;

	//	Approx translation into displacement parameters
	//	Derived
	vector<T>			myAlphas;
	vector<T>			myBetas;
    vector<Dynamics>    myDynamics;

	//	Finally, correl
	//	Lower diagonal matrix
	matrix<T>			myCorrel;
	T					myLambda;

	//	Derived
	matrix<T>			myUsedCorrel;
	//	Cholesky coefs, also lower diagonal
	matrix<T>			myChol;

    //  Similuation timeline = today + div dates + event dates
    vector<Time>		myTimeline;
    //  true (1) if the time step is an event date
    //  false (0) if it is an additional simulation step
    vector<bool>		myCommonSteps;
	//	true (1) if and only if this is a dividend step
	vector<bool>		myDivSteps;

    //  The pruduct's defline byref
    const vector<SampleDef>*    myDefline;

    //  Pre-calculated on initialization

    //  pre-calculated stds
    matrix<T>           myStds;		//	[time, asset]
    //  pre-calculated drifts 
    matrix<T>           myDrifts;	//	[time, asset]
    
    //  pre-caluclated numeraires exp(-r * t)
    vector<T>           myNumeraires;	//	[time]
    //  pre-calculated discounts exp(r * (T - t))
    vector<vector<T>>   myDiscounts;	//	[time]
    //  and rates = (exp(r * (T2 - T1)) - 1) / (T2 - T1)
    vector<vector<T>>   myLibors;	//	[time]

    //  forward factors 
    matrix<vector<T>>   myForwardFactors;	//	[time][asset][mats]

    //  Exported parameters
    vector<T*>          myParameters;
    vector<string>      myParameterLabels;

public:

    //  Constructor: store data

    template <class U>
    MultiDisplaced(
		const size_t		numAssets,
        const U             rate,
        const vector<U>&    spots,
		const vector<Time>&	divDates,
		const matrix<U>&	divs,
        const vector<U>&    atms,
        const vector<U>&    skews,
		const matrix<U>&	correl,
		const U&			lambda) : 
            myNumAssets(numAssets),
            myRate(rate),
			myDivDates(divDates),
			myDivs(divs),
			myLambda(lambda)
    {
		//	Copy data
		mySpots.resize(numAssets);
		copy(spots.begin(), spots.end(), mySpots.begin());
		myAtms.resize(numAssets);
		copy(atms.begin(), atms.end(), myAtms.begin());
		mySkews.resize(numAssets);
		copy(skews.begin(), skews.end(), mySkews.begin());
		myCorrel.resize(numAssets, numAssets);
		copy(correl.begin(), correl.end(), myCorrel.begin());
				
        //  Set parameter labels once 

		size_t numParams = 1 + numAssets + numAssets * divDates.size() + numAssets + numAssets + numAssets * (numAssets - 1) / 2 + 1;
		myParameters.resize(numParams);
		myParameterLabels.resize(numParams);

        myParameterLabels[0] = "rate";

		size_t paramNum = 1;

		for (size_t i = 0; i < numAssets; ++i)
		{
			myParameterLabels[paramNum] = string("spot ") + to_string(i + 1);
			++paramNum;
		}
		
		for (size_t i = 0; i < myDivDates.size(); ++i)
		{
			for (size_t j = 0; j < numAssets; ++j)
			{
                ostringstream ost;
                ost << setprecision(2) << fixed;
                ost << "div " << j << " " << myDivDates[i]; 
                myParameterLabels[paramNum] = ost.str();
				++paramNum;
			}		
		}

		for (size_t i = 0; i < numAssets; ++i)
		{
			myParameterLabels[paramNum] = string("ATM ") + to_string(i + 1);
			++paramNum;
		}

		for (size_t i = 0; i < numAssets; ++i)
		{
			myParameterLabels[paramNum] = string("skew ") + to_string(i + 1);
			++paramNum;
		}

		for (size_t i = 1; i < numAssets; ++i)
		{
			for (size_t j = 0; j < i; ++j)
			{
                myParameterLabels[paramNum] = string("correl ") + to_string(i + 1) + to_string(j + 1);
				++paramNum;
			}		
		}

		myParameterLabels[paramNum] = "lambda";

        setParamPointers();
    }

private:

    //  Must reset on copy
    void setParamPointers()
    {
        myParameters[0] = &myRate;

		size_t paramNum = 1;

		for (size_t i = 0; i < myNumAssets; ++i)
		{
			myParameters[paramNum] = & mySpots[i];
			++paramNum;
		}
		
		for (size_t i = 0; i < myDivDates.size(); ++i)
		{
			for (size_t j = 0; j < myNumAssets; ++j)
			{
                myParameters[paramNum] = & myDivs[i][j];
				++paramNum;
			}		
		}

		for (size_t i = 0; i < myNumAssets; ++i)
		{
			myParameters[paramNum] = & myAtms[i];
			++paramNum;
		}

		for (size_t i = 0; i < myNumAssets; ++i)
		{
			myParameters[paramNum] = & mySkews[i];
			++paramNum;
		}

		for (size_t i = 1; i < myNumAssets; ++i)
		{
			for (size_t j = 0; j < i; ++j)
			{
                myParameters[paramNum] = & myCorrel[i][j];
				++paramNum;
			}		
		}

		myParameters[paramNum] = & myLambda;
    }

public:

    //  Read access to parameters

	size_t numAssets() const override
	{
		return myNumAssets;
	}

    const T rate() const
    {
        return myRate;
    }

	const vector<T>& spots() const
    {
        return mySpots;
    }

	const vector<Time>& divDates() const
	{
		return myDivDates;
	}

	const matrix<T>& divs() const
	{
		return myDivs;
	}

    const vector<T>& atms() const
	{
		return myAtms;
	}

    const vector<T>& skews() const
	{
		return mySkews;
	}

    const matrix<T>& correl() const
	{
		return myCorrel;
	}

	const T lambda() const
	{
		return myLambda;
	}

    //  Access to all the model parameters
    const vector<T*>& parameters() override
    {
        return myParameters;
    }
    const vector<string>& parameterLabels() const override
    {
        return myParameterLabels;
    }

    //  Virtual copy constructor
    unique_ptr<Model<T>> clone() const override
    {
        auto clone = make_unique<MultiDisplaced<T>>(*this);
        clone->setParamPointers();
        return clone;
    }

    //  Initialize timeline
    void allocate(
        const vector<Time>&         productTimeline, 
        const vector<SampleDef>&    defline) 
            override
    {
		myAlphas.resize(myNumAssets);
		myBetas.resize(myNumAssets);
        myDynamics.resize(myNumAssets);

		myUsedCorrel.resize(myNumAssets, myNumAssets);		
		myChol.resize(myNumAssets, myNumAssets);

        //  Simulation timeline = today + div dates + product timeline

		set<Time> times;
		times.insert(systemTime);
        for (const auto& time : myDivDates)
        {
            times.insert(time);
        }
		for (const auto& time : productTimeline)
        {
            times.insert(time);
        }

        myTimeline.resize(times.size());
		copy(times.begin(), times.end(), myTimeline.begin());

        //  Mark steps on timeline that are on the product timeline
        myCommonSteps.resize(myTimeline.size());
        transform(myTimeline.begin(), myTimeline.end(), myCommonSteps.begin(), 
            [&](const Time t)
        {
            return binary_search(productTimeline.begin(), productTimeline.end(), t);
        });

        //  Mark div dates
        myDivSteps.resize(myTimeline.size());
        transform(myTimeline.begin(), myTimeline.end(), myDivSteps.begin(), 
            [&](const Time t)
        {
            return binary_search(myDivDates.begin(), myDivDates.end(), t);
        });

        //  Take a reference on the product's defline
        myDefline = &defline;

        //  Allocate the standard devs and drifts 
        //      over simulation timeline
        myStds.resize(myTimeline.size() - 1, myNumAssets);
        myDrifts.resize(myTimeline.size() - 1, myNumAssets);

        //  Allocate the numeraires, discount and forward factors 
        //      over product timeline
        const size_t n = productTimeline.size();
        myNumeraires.resize(n);
        
        myDiscounts.resize(n);
        for (size_t j = 0; j < n; ++j)
        {
            myDiscounts[j].resize(defline[j].discountMats.size());
        }

        myLibors.resize(n);
        for (size_t j = 0; j < n; ++j)
        {
            myLibors[j].resize(defline[j].liborDefs.size());
        }

        myForwardFactors.resize(n, myNumAssets);
        for (size_t j = 0; j < n; ++j) for (size_t k = 0; k < myNumAssets; ++k)
        {
            myForwardFactors[j][k].resize(defline[j].forwardMats[k].size());
        }
    }

    void init(
        const vector<Time>&         productTimeline, 
        const vector<SampleDef>&    defline) 
            override
    {
		//	Find alphas and betas out of ATMs and skews, see https://www.researchgate.net/publication/335146739_Displaced_Lognormal_Models
 
		for (size_t asset=0; asset<myNumAssets; ++asset)
		{
			myBetas[asset] = myAtms[asset] + 2 * mySkews[asset];
			if (fabs(myBetas[asset]) < 1.0e-05)
			{
                myDynamics[asset] = Normal;
			    myAlphas[asset] = - 2 * mySpots[asset] * mySkews[asset];	
				myBetas[asset] = 0.0;
			}
            else if (myBetas[asset] > 0)
            {
                myDynamics[asset] = Surnormal;
			    myAlphas[asset] = - 2 * mySpots[asset] / myBetas[asset] * mySkews[asset];	
            }
            else    //  neg beta
            {
                myDynamics[asset] = Subnormal;
                myBetas[asset] *= -1.0;   //  pos now
			    myAlphas[asset] = - 2 * mySpots[asset] / myBetas[asset] * mySkews[asset];	            
            }
		}

		//	Apply lambda factor to correl
		transform(myCorrel.begin(), myCorrel.end(), myUsedCorrel.begin(),
			[&](const T& rawCorr) { return myLambda * (1.0 - rawCorr) + rawCorr; });

		//	Apply Choleski decomposition (see e.g. Numerical Recipes) to correlation matrix
		//	If Wi are independent standard Gaussians, then:
		//	Bi = sum(CHOLij * Wj)
		//	are correlated standard Gaussians with the correct correl
		choldc(myUsedCorrel, myChol);

        //  Pre-compute the standard devs and drifts over simulation timeline        
        const T mu = myRate;

		const size_t n1 = myTimeline.size() - 1;
		const size_t n2 = myNumAssets;
        for (size_t i = 0; i < n1; ++i)
        {
            const double dt = myTimeline[i + 1] - myTimeline[i];

			for (size_t j = 0; j < n2; ++j)
			{
				myStds[i][j] = myBetas[j] * sqrt(dt);
				myDrifts[i][j] = (mu - 0.5*myBetas[j]*myBetas[j])*dt;
			}
        }

        //  Pre-compute the numeraires, discount and forward factors 
        //      on event dates
        const size_t m = productTimeline.size();

		for (size_t i = 0; i < m; ++i)
		{
			//  Numeraire
			if (defline[i].numeraire)
			{
                //  Under the risk neutral measure, 
                //      numeraire is deterministic in Black-Scholes = exp(rate * t)
                myNumeraires[i] = exp(myRate * productTimeline[i]);
			}

			//  Discount factors
			const size_t pDF = defline[i].discountMats.size();
			for (size_t j = 0; j < pDF; ++j)
			{
				myDiscounts[i][j] =
					exp(-myRate * (defline[i].discountMats[j] - productTimeline[i]));
			}

			//  Libors
			const size_t pL = defline[i].liborDefs.size();
			for (size_t j = 0; j < pL; ++j)
			{
				const double dt
					= defline[i].liborDefs[j].end - defline[i].liborDefs[j].start;
				myLibors[i][j] = (exp(myRate*dt) - 1.0) / dt;
			}

			//  Forward factors, beware the dividends
			
            /*
            const size_t pFF = defline[i].forwardMats.size();
			for (size_t j = 0; j < pFF; ++j)
			{
				myForwardFactors[i][j] =
					exp(mu * (defline[i].forwardMats[j] - productTimeline[i]));
			}
            */

		}   //  loop on event dates
	}

    //  MC Dimension
    size_t simDim() const override
    {
        return myNumAssets * (myTimeline.size() - 1);
    }

private:

    //  Helper function, fills a Sample given the spot
    inline void fillScen(
        const size_t        idx,    //  index on product timeline
        const T&            spot,   //  spot
        Sample<T>&          scen,   //  Sample to fill
        const SampleDef&    def)    //  and its definition
            const
    {
        if (def.numeraire)
        {
            scen.numeraire = myNumeraires[idx];
            if (mySpotMeasure) scen.numeraire *= spot;
        }
        
        for (size_t a = 0; a < myNumAssets; ++a)
        {
            transform(myForwardFactors[idx][a].begin(), myForwardFactors[idx][a].end(), 
                scen.forwards[a].begin(), 
                [&spot](const T& ff)
                {
                    return spot * ff;
                });
        }

        copy(myDiscounts[idx].begin(), myDiscounts[idx].end(), 
            scen.discounts.begin());

        copy(myLibors[idx].begin(), myLibors[idx].end(),
            scen.libors.begin());
    }

public:

    //  Generate one path, consume Gaussian vector
    //  path must be pre-allocated 
    //  with the same size as the product timeline
    void generatePath(
        const vector<double>&   gaussVec, 
        Scenario<T>&            path) 
            const override
    {
		/*

        //  The starting spot
        //  We know that today is on the timeline
        T spot = mySpot;
        //  Next index to fill on the product timeline
        size_t idx = 0;
        //  Is today on the product timeline?
        if (myTodayOnTimeline)
        {
            fillScen(idx, spot, path[idx], (*myDefline)[idx]);
            ++idx;
        }

        //  Iterate through timeline, apply sampling scheme
        const size_t n = myTimeline.size() - 1;
        for (size_t i = 0; i < n; ++i)
        {
            //  Apply known conditional distributions 
                //  Black-Scholes
            spot = spot * exp(myDrifts[i] 
                + myStds[i] * gaussVec[i]);
            //  Store on the path
            fillScen(idx, spot, path[idx], (*myDefline)[idx]);
            ++idx;
        }
	
		*/
    }
};