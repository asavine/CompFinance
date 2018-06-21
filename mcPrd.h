
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

#include "mcBase.h"

#define ONE_HOUR 0.000114469

template <class T>
class Europeans : public Product<T>
{
    vector<double>  myStrikes;
    vector<Time>    myMaturities;

public:

    //  Constructor: store data and build timeline
    Europeans(const vector<double>& strikes, const vector<Time>& maturities) : 
        myStrikes(strikes),
        myMaturities(maturities)
    {}

    //  Virtual copy constructor
    unique_ptr<Product<T>> clone() const override
    {
        return unique_ptr<Product<T>>(new Europeans<T>(*this));
    }

    //  Timeline
    const vector<Time>& timeline() const override
    {
        return myMaturities;
    }

    size_t numPayoffs() const override
    {
        return myMaturities.size() * myStrikes.size();
    }

    //  Payoffs, maturity major
    void payoffs(
        //  path, one entry per time step (on the product timeline)
        const vector<scenario<T>>&  path,
        //  pre-allocated space for resulting payoffs
        vector<T>&                  payoffs)
            const override
    {
        const size_t numT = myMaturities.size(), numK = myStrikes.size();
        
        auto payoffIt = payoffs.begin();
        for (size_t i = 0; i < numMat; ++i)
        {
            transform(
                myStrikes.begin(),
                myStrikes.end(),
                payoffIt,
                [spot = path[i].spot](const double& k) {return max<T>(spot - k, convert<T>(0.0)); }
            );

            payoffIt += numK;
        }
    }
};

template <class T>
class UOC : public Product<T>
{
    double          myStrike;
    double          myBarrier;
    Time            myMaturity;
    vector<Time>    myTimeline;

public:

    //  Constructor: store data and build timeline
    //  Timeline = system date to maturity, 
    //  with steps every monitoring frequency
    UOC(const double strike, 
        const double barrier, 
        const Time maturity, 
        const Time monitorFreq)
        : myStrike(strike), 
        myBarrier(barrier), 
        myMaturity(maturity)
    {
        myTimeline.push_back(systemTime);
        Time t = systemTime + monitorFreq;

        while (myMaturity - t > ONE_HOUR)
        {
            myTimeline.push_back(t);
            t += monitorFreq;
        }

        if (myTimeline.back() < myMaturity)
            myTimeline.push_back(myMaturity);
    }

    //  Virtual copy constructor
    unique_ptr<Product<T>> clone() const override
    {
        return unique_ptr<Product<T>>(new UOC<T>(*this));
    }

    //  Timeline
    const vector<Time>& timeline() const override
    {
        return myTimeline;
    }

    //  Payoff
    void payoffs(
        //  path, one entry per time step (on the product timeline)
        const vector<scenario<T>>&  path,
        //  pre-allocated space for resulting payoffs
        vector<T>&                  payoffs)
            const override
    {
        //  We apply the smooth barrier technique to stabilize risks
        //  See Savine's presentation on Fuzzy Logic, Global Derivatives 2016
        //  Or Andreasen and Savine's publication on scripting

        //  We apply a smoothing factor of 1% of the spot both ways, untemplated
        const double smooth = convert<double>(path[0].spot * 0.01);

        //  We start alive
        T alive = convert<T>(1.0);

        //  Go through path, update alive status
        for (const auto& scen: path)
        {
            //  Breached
            if (scen.spot > myBarrier + smooth)
            {
                payoffs[0] = convert<T>(0.0);
                return;
            }

            //  Semi-breached: apply smoothing
            if (scen.spot > myBarrier - smooth)
            {
                alive *= (myBarrier + smooth - scen.spot) / (2 * smooth);
            }
        }

        //  Payoff
        payoffs[0] = alive * max<T>(path.back().spot - myStrike, convert<T>(0.0));
    }
};

template <class T>
class Asians : public Product<T>
{
    vector<double>  myStrikes;
    Time            myMaturity;
    vector<Time>    myTimeline;

public:

    //  Constructor: store data and build timeline
    //  Timeline = system date to maturity, 
    //  with steps every monitoring frequency
    Asians(const vector<double>& strikes,
        const Time maturity,
        const Time monitorFreq)
        : myStrikes(strikes),
        myBarrier(barrier),
        myMaturity(maturity)
    {
        myTimeline.push_back(systemTime);
        Time t = systemTime + monitorFreq;

        while (myMaturity - t > ONE_HOUR)
        {
            myTimeline.push_back(t);
            t += monitorFreq;
        }

        if (myTimeline.back() < myMaturity)
            myTimeline.push_back(myMaturity);
    }

    //  Virtual copy constructor
    unique_ptr<Product<T>> clone() const override
    {
        return unique_ptr<Product<T>>(new Asians<T>(*this));
    }

    //  Timeline
    const vector<Time>& timeline() const override
    {
        return myTimeline;
    }

    size_t numPayoffs() const override
    {
        return myStrikes.size();
    }

    //  Payoff
    void payoffs(
        //  path, one entry per time step (on the product timeline)
        const vector<scenario<T>>&  path,
        //  pre-allocated space for resulting payoffs
        vector<T>&                  payoffs)
        const override
    {
        T ave = 0.0;
        size_t cntr = 0;

        //  Go through path, update average
        for (const auto& scen : path)
        {
            ave += scen.spot;
            ++cntr;
        }

        ave /= cntr;

        transform(
            myStrikes.begin(), 
            myStrikes.end(), 
            payoffs.begin(), 
            [ave](const double& k) {return max<T>(ave - k, convert<T>(0.0)); }
        );
    }
};