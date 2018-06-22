#pragma once

template <class T>
class SimpleBlackScholes : public Model<T>
{
    //  Model parameters

    //  Today's spot
    //  That would be today's linear market in a production system
    T                   mySpot;
    //  Local vols function of time
    T                   myVol;

    //  Normal specification = Bachelier
    //      or Lognormal = Black-Scholes
    const bool          myNormal;

    //  Similuation timeline = today + product timeline
    vector<Time>        myTimeline;
    bool                myTodayOnTimeline;

    //  Pre-calculated on initialization

    //  pre-calculated stds = vol * sqrt (dt)
    vector<T>           myStds;

    //  pre-calculated ito terms 0.5 * vol ^ 2 * dt
    vector<T>           myItoTerms;

    //  Exported parameters
    vector<T*>          myParameters;
    vector<string>      myParameterLabels;

public:

    //  Constructor: store data

    template <class U>
    SimpleBlackScholes(
        const U             spot,
        const U             vol,
        const bool          normal) : 
            mySpot(spot),
            myVol(vol),
            myNormal(normal),
            myParameters(2),
            myParameterLabels(2)
    {
        //  Set parameter labels once 
        myParameterLabels[0] = "spot";
        myParameterLabels[1] = "vol";

        setParamPointers();
    }

    //  Must reset on copy
    void setParamPointers()
    {
        myParameters[0] = &mySpot;
        myParameters[1] = &myVol;
    }

    //  Read access to parameters

    T spot() const
    {
        return mySpot;
    }

    const T vol() const
    {
        return myVol;
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
        auto clone = new SimpleBlackScholes<T>(*this);
        clone->setParamPointers();
        return unique_ptr<Model<T>>(clone);
    }

    //  Initialize timeline
    void init(const vector<Time>& productTimeline) override
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

        //  Allocate and compute the standard devs and ito terms

        myStds.resize(myTimeline.size() - 1);
        if (!myNormal) myItoTerms.resize(myTimeline.size() - 1);
        for (size_t i = 0; i < myTimeline.size() - 1; ++i)
        {
            const double dt = myTimeline[i + 1] - myTimeline[i];
            const double sqrtdt = sqrt(dt);
            myStds[i] = sqrtdt * myVol;
            if (!myNormal)
            {
                myItoTerms[i] = -0.5 * myStds[i] * myStds[i];
            }
        }
    }

    //  MC Dimension
    size_t simDim() const override
    {
        return myTimeline.size() - 1;
    }

    //  Generate one path, consume Gaussian vector
    //  path must be pre-allocated 
    //  with the same size as the product timeline
    void generatePath(const vector<double>& gaussVec, vector<scenario<T>>& path) const override
    {
        //  The starting spot
        //  We know that today is on the timeline
        T spot = mySpot;
        Time current = systemTime;
        //  Next index to fill on the product timeline
        size_t idx = 0;
        //  Is today on the product timeline?
        if (myTodayOnTimeline) path[idx++].spot = spot;

        //  Iterate through timeline
        for (size_t i = 1; i<myTimeline.size(); ++i)
        {
            //  Interpolate volatility in spot
            const T std = myStds[i - 1];
            //  vol comes out * sqrt(dt)

            //  Apply Euler's scheme
            if (myNormal)
            {
                //  Bachelier
                spot += std * gaussVec[i - 1];
            }
            else
            {
                //  Black-Scholes
                spot *= exp(myItoTerms[i - 1] + std * gaussVec[i - 1]);
            }

            //  Store on the path
            path[idx++].spot = spot;
        }
    }
};