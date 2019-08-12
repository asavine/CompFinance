
/*
Written by Antoine Savine in 2019

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications
is granted by explicit application, antoine@asavine.com

This comment must be preserved at all times
*/

//  Multiple correlated assets, all following a displaced Black-Scholes dynamics, subjects to discrete proportional dividends

//	Dynamics for each asset: dS / S = (alpha + S) * beta / S dW		--	i.e. local vol = (alpha + S) * beta
//	On a dividend date: S+ = (1-div) * S-							--	by conv, dividend is applied end of day, that is, after everything else

#pragma once

#include "matrix.h"
#include "choldc.h"

template <class T>
class MultiDisplaced : public Model<T>
{
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
    vector<double>		myAtms;
	//	skew
	vector<double>		mySkews;
	//	Approx translation into displacement parameters
	vector<T>			myAlphas;
	vector<T>			myBetas;

	//	Finally, correl
	//	Lower diagonal matrix
	matrix<double>		myCorrel;
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
    vector<T>           myStds;
    //  pre-calculated drifts 
    vector<T>           myDrifts;
    
    //  pre-caluclated numeraires exp(-r * t)
    vector<T>           myNumeraires;
    //  pre-calculated discounts exp(r * (T - t))
    vector<vector<T>>   myDiscounts;
    //  forward factors exp((r - d) * (T - t))
    vector<vector<T>>   myForwardFactors;
    //  and rates = (exp(r * (T2 - T1)) - 1) / (T2 - T1)
    vector<vector<T>>   myLibors;

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
		const matrix<U>&	correl) : 
            myNumAssets(numAssets),
            myRate(rate),
			myDivDates(divDates),
			myDivs(divs),
            mySpotMeasure(spotMeasure)
    {
		//	Copy data
		mySpots.resize(numAssets);
		copy(spots.begin(), spots.end(), mySpots.begin());
		myAtms.resize(numAssets);
		copy(atms.begin(), atms.end(), myAtms.begin());
		mySkews.resize(numAssets);
		copy(skews.begin(), skews.end(), mySkews.begin());
		myCorrel.resize(numAssets, numAssets);
		copy(correl.begin(), correl.end(), myCorrel.begin())

		//	Apply approx 
		//	atm = "local vol at todays spot" = alpha * beta / S0 + beta
		//	skew = "vol points for 1% increase in strike from S0" = 100 * S0 * d_sigma /dK = - 50 * alpha * beta / S0  
		//	beta = atm + skew / 50  , alpha = - skew * S0 / (50 * beta)
		//	See Savine, Notes from Volatility Lectures at Copenhagen Uni, part II

		myAlpha.resize(numAssets);
		myBeta.resize(numAssets);
		for (size_t asset=0; asset<numAssets; ++asset)
		{
			myBeta[asset] = myAtms[asset] + mySkews[asset] / 50;
			myAlpha[asset] = -mySkews[asset] * mySpots[asset] / 50 / myBeta[asset];
		}

		//	Apply Choleski decomposition (see e.g. Numerical Recipes) to correlation matrix
		//	If W0, W1, ..., Wn-1 are independent standard Gaussians,
		//	B0 = W0, B1 = chol[1][0] * W0 + chol[1][1] * W1, B2 = chol[2][0] * W0 + chol[2][1] * W1 + chol[2][2] * W2, ... are correlated standard Gaussians with the correct correl

		myChol.resize(numAssets, numAssets);
		choldc(correl, myChol);
				
        //  Set parameter labels once 

		size_t numParams = 1 + numAssets + numAssets * divdDates.size() + numAssets + numAssets + numAssets * (numAssets + 1) / 2;
		myParameters.resize(numParamns);
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
			myParameterLabels[paramNum] = string("alpha ") + to_string(i + 1);
			++paramNum;
		}

		for (size_t i = 0; i < numAssets; ++i)
		{
			myParameterLabels[paramNum] = string("beta ") + to_string(i + 1);
			++paramNum;
		}

		for (size_t i = 0; i < numAssets; ++i)
		{
			for (size_t j = 0; j <= i; ++j)
			{
                myParameterLabels[paramNum] = string("chol ") + to_string(i + 1) + to_string(j + 1);
				++paramNum;
			}		
		}

        setParamPointers();
    }

private:

    //  Must reset on copy
    void setParamPointers()
    {
        myParameters[0] = &myRate;

		size_t paramNum = 1;

		for (size_t i = 0; i < numAssets; ++i)
		{
			myParameters[paramNum] = & mySpots[i];
			++paramNum;
		}
		
		for (size_t i = 0; i < myDivDates.size(); ++i)
		{
			for (size_t j = 0; j < numAssets; ++j)
			{
                myParameters[paramNum] = & myDivs[i][j];
				++paramNum;
			}		
		}

		for (size_t i = 0; i < numAssets; ++i)
		{
			myParameters[paramNum] = & myAlphas[i];
			++paramNum;
		}

		for (size_t i = 0; i < numAssets; ++i)
		{
			myParameters[paramNum] = & myBetas[i];
			++paramNum;
		}

		for (size_t i = 0; i < numAssets; ++i)
		{
			for (size_t j = 0; j <= i; ++j)
			{
                myParameters[paramNum] = myChol[i][j];
				++paramNum;
			}		
		}
    }

public:

    //  Read access to parameters

	const size_t numAssets() const 
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

    const vector<double>& atms() const
	{
		return myAtms;
	}

    const vector<double>& skews() const
	{
		return mySkews;
	}

    const vector<T>& alphas() const
	{
		return myAlphas;
	}

	const vector<T>& betas() const
	{
		return myBetas;
	}

    const matrix<double>& correl() const
	{
		return myCorrel;
	}

    const matrix<T>& chol() const
	{
		return myChol;
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
		/*

        //  Simulation timeline = today + product timeline
        myTimeline.clear();
        myTimeline.push_back(systemTime);
        for (const auto& time : productTimeline)
        {
            if (time > systemTime) myTimeline.push_back(time);
        }

        //  Is today on the timeline?
        myTodayOnTimeline = (productTimeline[0] == systemTime);

        //  Take a reference on the product's defline
        myDefline = &defline;

        //  Allocate the standard devs and drifts 
        //      over simulation timeline
        myStds.resize(myTimeline.size() - 1);
        myDrifts.resize(myTimeline.size() - 1);

        //  Allocate the numeraires, discount and forward factors 
        //      over product timeline
        const size_t n = productTimeline.size();
        myNumeraires.resize(n);
        
        myDiscounts.resize(n);
        for (size_t j = 0; j < n; ++j)
        {
            myDiscounts[j].resize(defline[j].discountMats.size());
        }

        myForwardFactors.resize(n);
        for (size_t j = 0; j < n; ++j)
        {
            myForwardFactors[j].resize(defline[j].forwardMats.size());
        }

        myLibors.resize(n);
        for (size_t j = 0; j < n; ++j)
        {
            myLibors[j].resize(defline[j].liborDefs.size());
        }

		*/
    }

    void init(
        const vector<Time>&         productTimeline, 
        const vector<SampleDef>&    defline) 
            override
    {
		/*

        //  Pre-compute the standard devs and drifts over simulation timeline        
        const T mu = myRate - myDiv;
        const size_t n = myTimeline.size() - 1;

        for (size_t i = 0; i < n; ++i)
        {
            const double dt = myTimeline[i + 1] - myTimeline[i];

            //  Var[logST2 / ST1] = vol^2 * dt
            myStds[i] = myVol * sqrt(dt);
            
            if (mySpotMeasure)
            {
                //  under spot measure 
                //      E[logST2 / ST1] = logST1 + ( (r - d) + 0.5 * vol ^ 2 ) * dt
                myDrifts[i] = (mu + 0.5*myVol*myVol)*dt;
            }
            else
            {
                //  under risk neutral measure 
                //  E[logST2 / ST1] = logST1 + ( (r - d) - 0.5 * vol ^ 2 ) * dt
                myDrifts[i] = (mu - 0.5*myVol*myVol)*dt;
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
				if (mySpotMeasure)
				{
                    //  Under the spot measure, 
                    //      the numeraire is the spot with reinvested dividend
                    //      num(t) = spot(t) / spot(0) * exp(div * t)
                    //      we precalculate exp(div * t) / spot(0)
                    myNumeraires[i] = exp(myDiv * productTimeline[i]) / mySpot;
				}
				else
				{
                    //  Under the risk neutral measure, 
                    //      numeraire is deterministic in Black-Scholes = exp(rate * t)
                    myNumeraires[i] = exp(myRate * productTimeline[i]);
				}
			}

			//  Discount factors
			const size_t pDF = defline[i].discountMats.size();
			for (size_t j = 0; j < pDF; ++j)
			{
				myDiscounts[i][j] =
					exp(-myRate * (defline[i].discountMats[j] - productTimeline[i]));
			}

			//  Forward factors
			const size_t pFF = defline[i].forwardMats.size();
			for (size_t j = 0; j < pFF; ++j)
			{
				myForwardFactors[i][j] =
					exp(mu * (defline[i].forwardMats[j] - productTimeline[i]));
			}

			//  Libors
			const size_t pL = defline[i].liborDefs.size();
			for (size_t j = 0; j < pL; ++j)
			{
				const double dt
					= defline[i].liborDefs[j].end - defline[i].liborDefs[j].start;
				myLibors[i][j] = (exp(myRate*dt) - 1.0) / dt;
			}
		}   //  loop on event dates

		*/
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
        
        transform(myForwardFactors[idx].begin(), myForwardFactors[idx].end(), 
            scen.forwards.begin(), 
            [&spot](const T& ff)
            {
                return spot * ff;
            }
        );

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