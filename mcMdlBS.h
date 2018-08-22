
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

//  Black and Scholes model, chapter 6

#pragma once

template <class T>
class BlackScholes : public Model<T>
{
    //  Model parameters

    //  Today's spot
    //  That would be today's linear market in a production system
    T                   mySpot;
    //  Constant rate and dividend yield
    T                   myRate;
    T                   myDiv;
    
    //  Constant vol
    T                   myVol;

    //  false = risk neutral measure
    //  true = spot measure
    const bool          mySpotMeasure;

    //  Similuation timeline = today + event dates
    vector<Time>        myTimeline;
    //  Is today on product timeline?
    bool                myTodayOnTimeline;  
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
    BlackScholes(
        const U             spot,
        const U             vol,
        const bool          spotMeasure = false,
        const U             rate = U(0.0),
        const U             div = U(0.0)) : 
            mySpot(spot),
            myVol(vol),
            myRate(rate),
            myDiv(div),
            mySpotMeasure(spotMeasure),
            myParameters(4),
            myParameterLabels(4) 
    {
        //  Set parameter labels once 
        myParameterLabels[0] = "spot";
        myParameterLabels[1] = "vol";
        myParameterLabels[2] = "rate";
        myParameterLabels[3] = "div";

        setParamPointers();
    }

private:

    //  Must reset on copy
    void setParamPointers()
    {
        myParameters[0] = &mySpot;
        myParameters[1] = &myVol;
        myParameters[2] = &myRate;
		myParameters[3] = &myDiv;
    }

public:

    //  Read access to parameters

    T spot() const
    {
        return mySpot;
    }

    const T vol() const
    {
        return myVol;
    }

    const T rate() const
    {
        return myRate;
    }

    const T div() const
    {
        return myDiv;
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
        auto clone = make_unique<BlackScholes<T>>(*this);
        clone->setParamPointers();
        return clone;
    }

    //  Initialize timeline
    void allocate(
        const vector<Time>&         productTimeline, 
        const vector<SampleDef>&    defline) 
            override
    {
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
    }

    void init(
        const vector<Time>&         productTimeline, 
        const vector<SampleDef>&    defline) 
            override
    {
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
	}

    //  MC Dimension
    size_t simDim() const override
    {
        return myTimeline.size() - 1;
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
    }
};