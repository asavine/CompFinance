
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

#include "matrix.h"
#include "mcBase.h"
#include "interp.h"
#include "utility.h"

#include <string>
#include <iterator>
using namespace std;

template <class T>
class BlackScholes : public Model<T>
{
    //  Model parameters

    //  Today's spot
    //  That would be today's linear market in a production system
    T                   mySpot;
    //  Local volatility structure
    vector<Time>        myTimes;
    //  Local vols function of time
    vector<T>           myVols;

    //  Normal specification = Bachelier
    //      or Lognormal = Black-Scholes
    bool                myNormal;

    //  Similuation timeline
    vector<Time>        myTimeline;
    //  true (1) if the time step is on the product timeline
    //  false (0) if it is an additional simulation step
    //  note we use vector<int> because vector<bool> is broken in C++
    vector<int>         myCommonSteps;

    //  Pre-calculated on initialization

    //  volatilities pre-interpolated in time for each time step
    vector<T>           myInterpVols;
    //  volatilities as stored are multiplied by sqrt(dt) 
    //  so there is no need to do that during paths generation

    //  pre-calculated ito terms 0.5 * vol ^ 2 * dt
    vector<T>           myItoTerms;

    //  Exported parameters
    vector<T*>          myParameters;
    vector<string>      myParameterLabels;

public:

    //  Constructor: store data

    template <class U>
    BlackScholes(const U spot,
        const vector<Time> times,
        const matrix<U> vols,
        const bool normal)
        : mySpot(spot),
        myTimes(times),
        myVols(vols),
        myNormal(normal)
    {
        //  All parameters in a vector
        myParameters.reserve(myVols.size() + 1);
        myParameterLabels.reserve(myVols.size() + 1);

        myParameters.push_back(& mySpot);
        myParameterLabels.push_back("spot");

        for (size_t i = 0; i < myVols.size(); ++i)
        {
            myParameters.push_back(& myVols[i]);
            myParameterLabels.push_back("vol " + to_string(myTimes[j]));
        }
    }

    //  Read access to parameters
    
    T spot() const
    {
        return mySpot;
    }

    const vector<Time>& times() const
    {
        return myTimes;
    }

    const vector<T>& vols() const
    {
        return myVols;
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
        return unique_ptr<Model<T>>(new BlackScholes<T>(*this));
    }

    //  Initialize timeline
    void init(const vector<Time>& productTimeline) override
    {
        //  Simulation timeline = today + product timeline
        myTimeline.clear();
        myTimeline.push_back(systemTime);
        for (const auto time& : productTimeline)
        {
            if (time > systemTime) myTimeLine.push_back(time);
        }

        //  Mark steps on timeline that are on the product timeline
        myCommonSteps.resize(myTimeline.size());
        fill(myTimeline.begin(), myTimeline.end(), true);
        if (productTimeline[0] > systemTime) myCommonSteps[0] = false;

        //  Allocate and compute the local volatilities
        //      pre-interpolated in time and multiplied by sqrt(dt)
        myInterpVols.resize(myTimeline.size() - 1);
        if (!myNormal) myItoTerms.resize(myTimeline.size() - 1);
        for (size_t i = 0; i < myTimeline.size() - 1; ++i)
        {
            const double dt = myTimeline[i + 1] - myTimeline[i];
            const double sqrtdt = sqrt(dt);
            onst double vol = sqrtdt * interp(
                myTimes.begin(),
                myTimes.end(),
                myVols.begin(),
                myVols.end(),
                myTimeline[i]);
            myInterpVols[i] = vol;
            if (!myNormal)
            {
                myItoTerms[i] = -0.5 * vol * vol * dt;
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
        if (myCommonSteps[idx]) path[idx++].spot = spot;

        //  Iterate through timeline
        for (size_t i = 1; i<myTimeline.size(); ++i)
        {
            //  Interpolate volatility in spot
            const T vol = myInterpVols[i - 1];
            //  vol comes out * sqrt(dt)

            //  Apply Euler's scheme
            if (myNormal)
            {
                //  Bachelier
                spot += vol * gaussVec[i - 1];
            }
            else
            {
                //  Black-Scholes
                spot *= exp(myItoTerms[i] + vol * gaussVec[i - 1]);
            }

            //  Store on the path?
            if (myCommonSteps[i]) path[idx++].spot = spot;
        }
    }
};

template <class T>
class Dupire : public Model<T>
{
    //  Model parameters

    //  Today's spot
    //  That would be today's linear market in a production system
    T                   mySpot;
    //  Local volatility structure
    vector<double>      mySpots;
    vector<Time>        myTimes;
    //  Local vols
    //  Spot major: sigma(spot i, time j) = myVols[i][j]
    matrix<T>           myVols;

    //  Numerical parameters

    //  Maximum space between time steps
    Time                myMaxDt;

    //  Similuation timeline
    vector<Time>        myTimeline;
    //  true (1) if the time step is on the product timeline
    //  false (0) if it is an additional simulation step
    //  note we use vector<int> because vector<bool> is broken in C++
    vector<int>         myCommonSteps;

    //  Pre-calculated on initialization

    //  volatilities pre-interpolated in time for each time step
    //  here time major: iv(time i, spot j) = myInterpVols[i][j]
    matrix<T>           myInterpVols;
    //  volatilities as stored are multiplied by sqrt(dt) 
    //  so there is no need to do that during paths generation

    //  Exported parameters
    vector<T*>          myParameters;
    vector<string>      myParameterLabels;

public:

    //  Constructor: store data
    
    template <class U>
    Dupire(const U spot,
        const vector<double> spots, 
        const vector<Time> times, 
        const matrix<U> vols,
        const Time maxDt = 0.25)
        : mySpot(spot), 
        mySpots(spots), 
        myTimes(times), 
        myVols(vols), 
        myMaxDt(maxDt),
        myParameters(myVols.rows() * myVols.cols() + 1),
        myParameterLabels(myVols.rows() * myVols.cols() + 1)
    {
        //  Set parameter labels once 
        myParameterLabels[0] = "spot";

        size_t p = 0;
        for (size_t i = 0; i < myVols.rows(); ++i)
        {
            for (size_t j = 0; j < myVols.cols(); ++j)
            {
                myParameterLabels[++p] = 
                    "vol " + to_string(mySpots[i]) + " " + to_string(myTimes[j]);
            }
        }

        setParamPointers();
    }

    //  Must reset on copy
    void setParamPointers()
    {
        myParameters[0] = & mySpot;

        size_t p = 0;
        for (size_t i = 0; i < myVols.rows(); ++i)
        {
            for (size_t j = 0; j < myVols.cols(); ++j)
            {
                myParameters[++p] = & myVols[i][j];
            }
        }
    }

    //  Read access to parameters
    T spot() const
    {
        return mySpot;
    }

    const vector<double>& spots() const
    {
        return mySpots;
    }

    const vector<Time>& times() const
    {
        return myTimes;
    }

    const vector<T>& vols() const
    {
        return myVols;
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
        auto clone = new Dupire<T>(*this);
        clone->setParamPointers();
        return unique_ptr<Model<T>>(clone);
    }

    //  Initialize timeline
    void init(const vector<Time>& productTimeline) override
    {
        //  Fill from product timeline
        
        //  Do the fill
        myTimeline = fillData(
            productTimeline, // Original (product) timeline
            myMaxDt, // Maximum space allowed
            &systemTime, // Include system time
            &systemTime + 1,
            0.00136986301369863);  // Minimum distance = half day
        
        //  Mark steps on timeline that are on the product timeline
        myCommonSteps.resize(myTimeline.size());
        transform(myTimeline.begin(), myTimeline.end(), myCommonSteps.begin(), 
            [&](const Time t)
        {
            return binary_search(productTimeline.begin(), productTimeline.end(), t);
        });

        //

        //  Allocate and compute the local volatilities
        //      pre-interpolated in time and multiplied by sqrt(dt)
        myInterpVols.resize(myTimeline.size() - 1, mySpots.size());
        for (size_t i = 0; i < myTimeline.size() - 1; ++i)
        {
            const double sqrtdt = sqrt(myTimeline[i+1] - myTimeline[i]);
            for (size_t j = 0; j < mySpots.size(); ++j)
            {
                myInterpVols[i][j] = sqrtdt * interp(
                    myTimes.begin(), 
                    myTimes.end(), 
                    myVols[j], 
                    myVols[j] + myTimes.size(), 
                    myTimeline[i]);
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
        if (myCommonSteps[idx]) path[idx++].spot = spot;

        //  Iterate through timeline
        for(size_t i=1; i<myTimeline.size(); ++i)
        {
            //  Interpolate volatility in spot
            const T vol = interp(
                mySpots.begin(), 
                mySpots.end(), 
                myInterpVols[i-1], 
                myInterpVols[i-1] + mySpots.size(), 
                spot);
            //  vol comes out * sqrt(dt)

            //  Apply Euler's scheme
            spot += vol * gaussVec[i - 1];

            //  Store on the path?
            if (myCommonSteps[i]) path[idx++].spot = spot;
        }
    }
};

//  Calibration

#include "ivs.h"

//  Calibrates one maturity
//  Main calibration function below
template <class IT, class OT, class T = double>
inline void dupireCalibMaturity(
    //  IVS we calibrate to
    const IVS& ivs,
    //  Maturity to calibrate
    const Time maturity,
    //  Spots for local vol
    IT spotsBegin,
    IT spotsEnd,
    //  Results, by spot
    //  With (random access) iterator, STL style
    OT  lVolsBegin,
    //  Risk view
    const RiskView<T>& riskView = RiskView<double>())
{
    //  Number of spots
    IT spots = spotsBegin;
    const size_t nSpots = distance(spotsBegin, spotsEnd);

    //  Estimate ATM so we cut the grid 2 stdevs away to avoid instabilities
    const double atmCall = convert<double>(ivs.call(ivs.spot(), maturity));
    //  Standard deviation, approx. atm call * sqrt(2pi)
    const double std = atmCall * 2.506628274631;

    //  Skip spots below and above 2.5 std
    int il = 0;
    while (spots[il] < ivs.spot() - 2.5 * std) ++il;
    int ih = nSpots - 1;
    while (spots[ih] > ivs.spot() + 2.5 * std) --ih;

    //  Loop on spots
    for (int i = il; i <= ih; ++i)
    {
        //  Dupire's formula
        lVolsBegin[i] = ivs.localVol(spots[i], maturity, &riskView);
    }

    //  Extrapolate flat outside std
    for (int i = 0; i <= il - 1; ++i)
        lVolsBegin[i] = lVolsBegin[il];
    for (int i = ih + 1; i < nSpots; ++i)
        lVolsBegin[i] = lVolsBegin[ih - 1];
}

//  Returns a struct with spots, times and lVols
template<class T = double>
inline auto dupireCalib(
        //  The IVS we calibrate to
        const IVS& ivs,
        //  The local vol grid
        //  The spots to include
        const vector<double>& inclSpots,
        //  Maximum space between spots
        const double maxDs,
        //  The times to include, note NOT 0
        const vector<Time>& inclTimes,
        //  Maximum space between times
        const double maxDt,
        //  Risk view if required
        //  omitted: T = double , no risk view
        const RiskView<T>& riskView = RiskView<double>())
{
    //  Results
    struct
    {
        vector<double> spots;
        vector<Time> times;
        matrix<T> lVols;
    } results;

    //  Spots and times
    results.spots = fillData(inclSpots, maxDs);
    results.times = fillData(inclTimes, maxDt, &maxDt, &maxDt + 1);

    //  Allocate local vols, transposed maturity first
    matrix<T> lVolsT(results.times.size(), results.spots.size());

    //  Maturity by maturity
    for (size_t j = 0; j < results.times.size(); ++j)
    {
        dupireCalibMaturity(
            ivs, 
            results.times[j], 
            results.spots.begin(),
            results.spots.end(),
            lVolsT[j], 
            riskView);
    }

    //  transpose is defined in matrix.h
    results.lVols = transpose(lVolsT);

    return results;
}
