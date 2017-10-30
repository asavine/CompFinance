#pragma once

#include "matrix.h"
#include "mcBase.h"
#include "interp.h"

//  Calibration helper
static const struct CalibrationType {} calibrate;

#define EPS 1.0e-08
//  Utility for filling schedules
inline void fillTimeline(
    //  The original timeline
    const vector<Time>&             original,
    //  The maximum spacing allowed
    const Time&                     maxDt,
    //  The filled timeline
    //      Has all original steps 
    //      Plus additional ones so maxDt is not exceeded
    //  Size in unknown to the caller, hence allocation occurs inside
    vector<Time>&                   filled,
    //  True if the point was on the original timeline
    //  We use vector<int> because vector<bool> is broken in C++
    vector<int>&                    common)
{
    //  Position on the start of the product timeline
    auto it = original.begin();

    //  Skip before today
    while (it != original.end() && *it < systemTime) ++it;
    if (it == original.end())
        throw runtime_error
        ("All dates on the timeline are in the past");

    //  Clear
    filled.clear();
    common.clear();

    //  Include today
    filled.push_back(systemTime);

    //  Is today on the product timeline?
    if (*it == systemTime)
    {
        common.push_back(true);
        //  Advance
        ++it;
    }
    else
    {
        common.push_back(false);
    }

    //  Futures dates
    while (it != original.end())
    {
        Time current = filled.back();
        Time next = *it;
        //  Must supplement?
        if (next - current > maxDt)
        {
            //  Number of time steps to add
            int addSteps = int((next - current) / maxDt - EPS) + 1;
            //  Spacing between supplementary points
            Time spacing = (next - current) / addSteps;
            //  Add the steps
            Time t = current + spacing;
            while (t < next)
            {
                filled.push_back(t);
                common.push_back(false);
                t += spacing;
            }
        }
        //  Push the next step on the product timeline and advance
        filled.push_back(*it);
        common.push_back(true);
        ++it;
    }
}

//  Model

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
        fillTimeline(productTimeline, myMaxDt, myTimeline, myCommonSteps);
        
        //  Allocate and compute the local volatilities
        //      pre-interpolated in time
        myInterpVols.resize(myTimeline.size() - 1, mySpots.size());
        for (size_t i = 0; i < myTimeline.size() - 1; ++i)
        {
            const double sqrtdt = sqrt(myTimeline[i+1] - myTimeline[i]);
            for (size_t j = 0; j < mySpots.size(); ++j)
            {
                myInterpVols[i][j] = sqrtdt * linterp(
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
            //  Interpolate volatility
            const T vol = linterp(
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

//  Returns spots, times and local vols
template<class T>
inline tuple<vector<double>, vector<Time>, matrix<T>> 
    dupireCalib(
        //  The IVS we calibrate to
        const IVS<T>& ivs,
        //  The local vol grid
        const double& lowSpot,
        const double& highSpot,
        const size_t numSpots,
        const Time& finalMat,
        const size_t numTimes)
{
    //  Spots and times
    
    const double spot = convert<double>(ivs.spot());

    vector<double> spots(numSpots + 1);
    const double ds = (highSpot - lowSpot) / numSpots;
    spots[0] = lowSpot;
    for (size_t i = 0; i < numSpots; ++i) spots[i + 1] = spots[i] + ds;

    vector<Time> times(numTimes);
    const double dt = finalMat / numTimes;
    times[0] = dt;
    for (size_t i = 0; i < numTimes - 1; ++i) times[i + 1] = times[i] + dt;

    //  Allocate local vols
    matrix<T> lVol(spots.size(), times.size());

    //  Maturity by maturity
    for (j = 0; j < times.size(); ++j)
    {
        //  ATM
        const double atmCall = convert<double>(ivs.call(spot, times[j]);
        //  Standard deviation, approx. atm call * sqrt(2pi)
        const double std = atmCall * 2.506628274631;

        //  Skip spots below and above 3 std
        size_t il = 0;
        while (spots[il] < spot - 3 * std) ++il;
        size_t ih = spots.size() - 1;
        while (spots[ih] > spot + 3 * std) --ih;

        //  Loop on spots
        for (size_t i = il; i <= ih; ++i)
        {
            //  Dupire's formula
            lVol[i][j] = sqrt(2.0 * ivs.cT(strikes[i], mat[j]) 
                / ivs.cKK(strikes[i], mat[j]));
        }

        //  Extrapolate flat outside std
        for (size_t i = 0; i < il - 1; ++i) 
            lVol[i][j] = lVol[il][j];
        for (size_t i = ih + 1; i < spots.size().rows(); ++i)
            lVol[i][j] = lVol[ih - 1][j - j0];

    }

    return make_tuple(spots, times, lVol);
}
