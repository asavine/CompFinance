#pragma once

#include "matrix.h"
#include "mcBase.h"
#include "interp.h"
#include "utility.h"

#include <iterator>

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

    //  Simulation timeline

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
        myMaxDt(maxDt)
    { }

    //  Read access to parameters
    T spot() const
    {
        return mySpot;
    }

    const matrix<T>&vols()
    {
        return myVols;
    }

    //  Virtual copy constructor
    unique_ptr<Model<T>> clone() const override
    {
        return unique_ptr<Model<T>>(new Dupire<T>(*this));
    }

    //  Initialize timeline
    void init(const vector<Time>& productTimeline) override
    {
        //  Fill from product timeline
        
        //  Do the fill
        myTimeline = fillData<Time>(
            productTimeline, // Original (product) timeline
            myMaxDt, // Maximum space allowed
            &vector<Time>(1, systemTime), // Include system time
            0.002739726);  // Minimum distance = 1 day
        
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

    //  Access to all parameters by copy
    vector<T> parameters() const override
    {
        vector<T> params;
        params.reserve(myVols.rows() * myVols.cols() + 1);
        params.push_back(mySpot);
        for (auto& vol : myVols) params.push_back(vol);

        return params;
    }

    //  AAD enabled
private:
    //  Implementation, for Number instances only

    //  Put parameters on tape 
    template <class U> void putOnTapeI() {}
    template<> void putOnTapeI<Number>()
    {
        mySpot.putOnTape();
        //  Free function putOnTape in AADNumber.h
        ::putOnTape(myVols.begin(), myVols.end());
    }

public:
    //  Interface    

    //  Put parameters on tape 
    void putOnTape() override
    {
        putOnTapeI<T>();
    }
};

//  Calibration

#include "ivs.h"

//  Calibrates one maturity
//  Main calibration function below
template <class IT, class T = double>
inline void dupireCalibMaturity(
    //  IVS we calibrate to
    const IVS& ivs,
    //  Maturity to calibrate
    const Time maturity,
    //  Spots for local vol
    const vector<double> spots,
    //  Results, by spot
    //  With (random access) iterator, STL style
    IT lVolsBegin,
    //  Risk view
    const RiskView<T>& riskView = RiskView<double>())
{
    //  Estimate ATM so we cut the grid 2 stdevs away to avoid instabilities
    const double atmCall = convert<double>(ivs.call(ivs.spot(), maturity));
    //  Standard deviation, approx. atm call * sqrt(2pi)
    const double std = atmCall * 2.506628274631;

    //  Skip spots below and above 2.5 std
    int il = 0;
    while (spots[il] < ivs.spot() - 2.5 * std) ++il;
    int ih = spots.size() - 1;
    while (spots[ih] > ivs.spot() + 2.5 * std) --ih;

    //  Loop on spots
    for (int i = il; i <= ih; ++i)
    {
        //  Dupire's formula
        lVolsBegin[i] = ivs.localVol(spots[i], maturity, riskView);
    }

    //  Extrapolate flat outside std
    for (int i = 0; i <= il - 1; ++i)
        lVolsBegin[i] = lVolsBegin[il];
    for (int i = ih + 1; i < spots.size(); ++i)
        lVolsBegin[i] = lVolsBegin[ih - 1];
}

//  Returns spots, times and local vols
template<class T = double>
inline tuple<vector<double>, vector<Time>, matrix<T>> 
    dupireCalib(
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
    //  Spots and times
    vector<double> spots = fillData(inclSpots, maxDs);
    vector<Time> times = fillData(inclTimes, maxDt, &vector<Time>(1, maxDt));

    //  Allocate local vols, transposed maturity first
    matrix<T> lVolsT(times.size(), spots.size());

    //  Maturity by maturity
    for (size_t j = 0; j < times.size(); ++j)
    {
        dupireCalibMaturity(
            ivs, 
            times[j], 
            spots, 
            lVolsT[j], 
            riskView);
    }

    //  transpose is defined in matrix.h
    return make_tuple(spots, times, transpose(lVolsT));
}
