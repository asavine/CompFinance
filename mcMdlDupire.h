
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

//  Dupire's model, chapter 6

#pragma once
#pragma warning(disable : 4018)

#include "matrix.h"
#include "interp.h"
#include "utility.h"

#define HALF_DAY 0.00136986301369863

template <class T>
class Dupire : public Model<T>
{
    //  Model parameters

    //  Today's spot
    //  That would be today's linear market in a production system
    T                       mySpot;
    //  Local volatility structure
    const vector<double>    mySpots;
    //  We keep log spots to interpolate in log space
    vector<double>          myLogSpots;
    const vector<Time>      myTimes;
    //  Local vols
    //  Spot major: sigma(spot i, time j) = myVols[i][j]
    matrix<T>               myVols;

    //  Numerical parameters

    //  Maximum space between time steps
    const Time              myMaxDt;

    //  Similuation timeline
    vector<Time>            myTimeline;
    //  true (1) if the time step is an event date
    //  false (0) if it is an additional simulation step
    vector<bool>            myCommonSteps;

    //  The pruduct's defline byref
    const vector<SampleDef>*    myDefline;

    //  Pre-calculated on initialization

    //  volatilities pre-interpolated in time for each time step
    //  here time major: iv(time i, spot j) = myInterpVols[i][j]
    matrix<T>               myInterpVols;
    //  volatilities as stored are multiplied by sqrt(dt) 
    //  so there is no need to do that during paths generation

    //  Exported parameters
    vector<T*>              myParameters;
    vector<string>          myParameterLabels;

public:

    //  Constructor: store data
    
    template <class U>
    Dupire(const U              spot,
        const vector<double>    spots,
        const vector<Time>      times,
        const matrix<U>         vols,
        const Time maxDt =      0.25)
        : mySpot(spot),
        mySpots(spots),
        myLogSpots(mySpots.size()),
        myTimes(times),
        myVols(vols),
        myMaxDt(maxDt),
        myParameters(myVols.rows() * myVols.cols() + 1),
        myParameterLabels(myVols.rows() * myVols.cols() + 1)
    {
        //  Compute log spots
		transform(mySpots.begin(), mySpots.end(), myLogSpots.begin(), 
            [](const double s) {return log(s); });

        //  Set parameter labels once 
        myParameterLabels[0] = "spot";

        size_t p = 0;
        for (size_t i = 0; i < myVols.rows(); ++i)
        {
            for (size_t j = 0; j < myVols.cols(); ++j)
            {
                ostringstream ost;
                ost << setprecision(2) << fixed;
                ost << "lvol " << mySpots[i] << " " << myTimes[j];
 
                myParameterLabels[++p] = ost.str();
            }
        }

        setParamPointers();
    }

private:

    //  Must reset on copy
    void setParamPointers()
    {
        myParameters[0] = &mySpot;
        transform(myVols.begin(), myVols.end(), next(myParameters.begin()), 
            [](auto& vol) {return &vol; });
    }

public:

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
        auto clone = make_unique<Dupire<T>>(*this);
        clone->setParamPointers();
        return clone;
    }

    //  Initialize timeline
    void allocate(
        const vector<Time>&         productTimeline, 
        const vector<SampleDef>&    defline) 
            override
    {
        //  Fill from product timeline
        
        //  Do the fill
        myTimeline = fillData(
            productTimeline, // Original (product) timeline
            myMaxDt, // Maximum space allowed
            HALF_DAY, // Minimum distance = half day
            &systemTime, &systemTime + 1);  //  Hack to include system time
        
        //  Mark steps on timeline that are on the product timeline
        myCommonSteps.resize(myTimeline.size());
        transform(myTimeline.begin(), myTimeline.end(), myCommonSteps.begin(), 
            [&](const Time t)
        {
            return binary_search(productTimeline.begin(), productTimeline.end(), t);
        });

        //  Take a reference on the product's defline
        myDefline = &defline;

        //  Allocate the local volatilities
        //      pre-interpolated in time over simulation timeline
        myInterpVols.resize(myTimeline.size() - 1, mySpots.size());
    }

    void init(
        const vector<Time>&         productTimeline, 
        const vector<SampleDef>&    defline) 
            override
    {
        //  Compute the local volatilities
        //      pre-interpolated in time and multiplied by sqrt(dt)
        const size_t n = myTimeline.size() - 1;
        for (size_t i = 0; i < n; ++i)
        {
            const double sqrtdt = sqrt(myTimeline[i + 1] - myTimeline[i]);
            const size_t m = myLogSpots.size();
            for (size_t j = 0; j < m; ++j)
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

private:

    //  Helper function, fills a sample given the spot
    inline static void fillScen(const T& spot, Sample<T>& scen)
    {
        fill(scen.forwards.begin(), scen.forwards.end(), spot);
    }

public:

    //  Generate one path, consume Gaussian vector
    //  path must be pre-allocated 
    //  with the same size as the product timeline
    void generatePath(
        const vector<double>& gaussVec, 
        Scenario<T>& path) 
            const override
    {
        //  The starting spot
        //  We know that today is on the timeline
        T logspot = log(mySpot);
        Time current = systemTime;
        //  Next index to fill on the product timeline
        size_t idx = 0;
        //  Is today on the product timeline?
        if (myCommonSteps[idx])
        {
            fillScen(exp(logspot), path[idx]);
            ++idx;
        }

        //  Iterate through timeline
        const size_t n = myTimeline.size() - 1;
        const size_t m = myLogSpots.size();
        for (size_t i = 0; i < n; ++i)
        {
            //  Interpolate volatility in spot
            T vol = interp(
                myLogSpots.begin(),
                myLogSpots.end(),
                myInterpVols[i],
                myInterpVols[i] + m,
                logspot);
            //  vol comes out * sqrt(dt)

            //  Apply Euler's scheme
            logspot += vol * (- 0.5 * vol + gaussVec[i]);

            //  Store on the path?
            if (myCommonSteps[i + 1])
            {
                fillScen(exp(logspot), path[idx]);
                ++idx;
            }
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
    const double atmCall = double(ivs.call(ivs.spot(), maturity));
    //  Standard deviation, approx. atm call * sqrt(2pi)
    const double std = atmCall * 2.506628274631;

    //  Skip spots below and above 2.5 std
    int il = 0;
    while (il < nSpots && spots[il] < ivs.spot() - 2.5 * std) ++il;
    int ih = nSpots - 1;
    while (ih >= 0 && spots[ih] > ivs.spot() + 2.5 * std) --ih;

    //  Loop on spots
    for (int i = il; i <= ih; ++i)
    {
        //  Dupire's formula
        lVolsBegin[i] = ivs.localVol(spots[i], maturity, &riskView);
    }

    //  Extrapolate flat outside std
    for (int i = 0; i < il; ++i)
        lVolsBegin[i] = lVolsBegin[il];
    for (int i = ih + 1; i < nSpots; ++i)
        lVolsBegin[i] = lVolsBegin[ih];
}

#define ONE_HOUR 0.000114469

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
    results.spots = fillData(inclSpots, maxDs, 0.01); //  min space = 0.01
    results.times = fillData(inclTimes, maxDt,
        ONE_HOUR,               //  min space = 1 hour
        &maxDt, &maxDt + 1      //  dirty trick to include maxDt
    );

    //  Allocate local vols, transposed maturity first
    matrix<T> lVolsT(results.times.size(), results.spots.size());

    //  Maturity by maturity
    const size_t n = results.times.size();
    for (size_t j = 0; j < n; ++j)
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
