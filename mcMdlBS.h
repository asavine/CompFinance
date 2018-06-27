#pragma once

template <class T>
class BlackScholes : public Model<T>
{
    //  Model parameters

    //  Today's spot
    //  That would be today's linear market in a production system
    T                   mySpot;
    //  Normal : in units of spot, Lognormal : in %
    T                   myVol;

    //  Constant rate and dividend yield
    T                   myRate;
    T                   myDiv;

    //  Normal specification = Bachelier
    //      or Lognormal = Black-Scholes
    const bool          myNormal;

    //  Similuation timeline = today + product timeline
    vector<Time>        myTimeline;
    bool                myTodayOnTimeline;  //  Is today on product timeline?
    //  The pruduct's dataline byref
    const vector<simulData>*    myDataline;

    //  Pre-calculated on initialization

    //  pre-calculated stds
    vector<T>           myStds;
    //  pre-calculated drifts 
    vector<T>           myDrifts;
    
    //  pre-caluclated numeraires exp(-r * t)
    vector<T>           myNumeraires;
    //  pre-calculated discounts exp(r * (T - t))
    vector<vector<T>>   myDiscounts;
    //  and forward factors exp((r - d) * (T - t))
    vector<vector<T>>   myForwardFactors;

    //  Exported parameters
    vector<T*>          myParameters;
    vector<string>      myParameterLabels;

public:

    //  Constructor: store data

    template <class U>
    BlackScholes(
        const U             spot,
        const U             vol,
        const bool          normal,
        const U             rate = U(0.0),
        const U             div = U(0.0)) : 
            mySpot(spot),
            myVol(vol),
            myRate(rate),
            myDiv(div),
            myNormal(normal),
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
        auto clone = new BlackScholes<T>(*this);
        clone->setParamPointers();
        return unique_ptr<Model<T>>(clone);
    }

    //  Initialize timeline
    void allocate(const vector<Time>& productTimeline, const vector<simulData>& dataline) override
    {
        //  Simulation timeline = today + product timeline
        myTimeline.clear();
        myTimeline.push_back(systemTime);
        for (const auto& time : productTimeline)
        {
            if (time > systemTime) myTimeline.push_back(time);
        }

        //  Mark steps on timeline that are on the product timeline
        myTodayOnTimeline = (productTimeline[0] == systemTime);

        //  Take a reference on the product's dataline
        myDataline = &dataline;

        //  Allocate the standard devs and drifts over simulation timeline
        myStds.resize(myTimeline.size() - 1);
        myDrifts.resize(myTimeline.size() - 1);

        //  Allocate the numeraires, discount and forward factors over product timeline
        const size_t n = productTimeline.size();
        myNumeraires.resize(n);
        
        myDiscounts.resize(n);
        for (size_t j = 0; j < n; ++j)
        {
            myDiscounts[j].resize(dataline[j].discountMats.size());
        }

        myForwardFactors.resize(productTimeline.size());
        for (size_t j = 0; j < n; ++j)
        {
            myForwardFactors[j].resize(dataline[j].forwardMats.size());
        }
    }

    void init(const vector<Time>& productTimeline, const vector<simulData>& dataline) override
    {
        //  Pre-compute the standard devs and drifts over simulation timeline        
        const T mu = myRate - myDiv;
        const size_t n = myTimeline.size() - 1;
        if (myNormal)
        {
            for (size_t i = 0; i < n; ++i)
            {
                const double dt = myTimeline[i + 1] - myTimeline[i];

                //  normal model
                //  Var[ST2 / ST1] = vol^2 * ( exp ( 2 * (r - d) * dt ) - 1) / ( 2 * (r - d) ) 
                //  E[ST2 / ST1] = ST1 * exp ( (r - d) * dt )
                myStds[i] = fabs(mu) > EPS
                    ? T(myVol * sqrt((exp(2 * mu*dt) - 1) / (2 * mu)))
                    : T(myVol * sqrt(dt));
                myDrifts[i] = exp((myRate - myDiv)*dt);
            }
        }
        else
        {
            for (size_t i = 0; i < n; ++i)
            {
                const double dt = myTimeline[i + 1] - myTimeline[i];

                //  lognormal model
                //  Var[logST2 / ST1] = vol^2 * dt
                //  E[logST2 / ST1] = logST1 + ( (r - d) - 0.5 * vol ^ 2 ) * dt
                myStds[i] = myVol * sqrt(dt);
                myDrifts[i] = (mu - 0.5*myVol*myVol)*dt;
            }
        }

        //  Pre-compute the numeraires, discount and forward factors over product timeline
        const size_t m = productTimeline.size();

        for (size_t i = 0; i < m; ++i)
        {
            //  Numeraire
            myNumeraires[i] = exp(myRate * productTimeline[i]);
        }

        for (size_t i = 0; i < m; ++i)
        {
            //  Discount factors
            const size_t p = dataline[i].discountMats.size();
            for (size_t j = 0; j < p; ++j)
            {
                myDiscounts[i][j] = exp(-myRate * (dataline[i].discountMats[j] - productTimeline[i]));
            }
        }

        for (size_t i = 0; i < m; ++i)
        {
            //  Discount factors
            const size_t p = dataline[i].forwardMats.size();
            //  Forward factors
            for (size_t j = 0; j < p; ++j)
            {
                myForwardFactors[i][j] = exp(mu * (dataline[i].forwardMats[j] - productTimeline[i]));
            }
        }
    }

    //  MC Dimension
    size_t simDim() const override
    {
        return myTimeline.size() - 1;
    }

private:

    //  Helper function, fills a scenario given the spot
    inline void fillScen(
        const size_t idx,   //  index on product timeline
        const T& spot,      //  spot
        scenario<T>& scen   //  scenario to fill
    ) const
    {
        scen.numeraire = myNumeraires[idx];
        
        transform(myForwardFactors[idx].begin(), myForwardFactors[idx].end(), 
            scen.forwards.begin(), 
            [&spot](const T& ff)
            {
                return spot * ff;
            }
        );

        copy(myDiscounts[idx].begin(), myDiscounts[idx].end(), 
            scen.discounts.begin());
    }

public:

    //  Generate one path, consume Gaussian vector
    //  path must be pre-allocated 
    //  with the same size as the product timeline
    void generatePath(const vector<double>& gaussVec, vector<scenario<T>>& path) const override
    {
        //  The starting spot
        //  We know that today is on the timeline
        T spot = mySpot;
        //  Next index to fill on the product timeline
        size_t idx = 0;
        //  Is today on the product timeline?
        if (myTodayOnTimeline)
        {
            fillScen(idx, spot, path[idx]);
            ++idx;
        }

        const size_t n = myTimeline.size() - 1;
        if (myNormal)
        {
            //  Iterate through timeline
            for (size_t i = 0; i < n; ++i)
            {
                //  Apply known conditional distributions 
                    //  Bachelier
                spot = spot * myDrifts[i]
                    + myStds[i] * gaussVec[i];
                //  Store on the path
                fillScen(idx, spot, path[idx]);
                ++idx;
            }
        }
        else
        {
            //  Iterate through timeline
            for (size_t i = 0; i < n; ++i)
            {
                //  Apply known conditional distributions 
                    //  Black-Scholes
                spot = spot * exp(myDrifts[i] 
                    + myStds[i] * gaussVec[i]);
                //  Store on the path
                fillScen(idx, spot, path[idx]);
                ++idx;
            }
        }
    }
};