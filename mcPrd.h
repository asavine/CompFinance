
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

#include <map>

#include "mcBase.h"

#define ONE_HOUR 0.000114469

template <class T>
class European : public Product<T>
{
    double              myStrike;
    Time                myExerciseDate;
    Time                mySettlementDate;

    vector<Time>        myTimeline;
    vector<simulData>   myDataline;

public:

    //  Constructor: store data and build timeline
    European(const double strike, 
        const Time exerciseDate,
        const Time settlementDate) :
        myStrike(strike),
        myExerciseDate(exerciseDate),
        mySettlementDate(settlementDate)
    {
        myTimeline.push_back(exerciseDate);

        myDataline.resize(1);
        myDataline[0].forwardMats.push_back(settlementDate);
        myDataline[0].discountMats.push_back(settlementDate);
    }

    European(const double strike,
        const Time exerciseDate) : European(strike, exerciseDate, exerciseDate)
    {}

    //  Virtual copy constructor
    unique_ptr<Product<T>> clone() const override
    {
        return unique_ptr<Product<T>>(new European<T>(*this));
    }

    //  Timeline
    const vector<Time>& timeline() const override
    {
        return myTimeline;
    }

    //  Dataline
    const vector<simulData>& dataline() const override
    {
        return myDataline;
    }

    //  Payoffs, maturity major
    void payoffs(
        //  path, one entry per time step (on the product timeline)
        const vector<scenario<T>>&  path,
        //  pre-allocated space for resulting payoffs
        vector<T>&                  payoffs)
            const override
    {
        payoffs[0] = max(path[0].forwards[0] - myStrike, 0.0)
            * path[0].discounts[0]
            / path[0].numeraire; 
    }
};

template <class T>
class UOC : public Product<T>
{
    double              myStrike;
    double              myBarrier;
    Time                myMaturity;
    vector<Time>        myTimeline;
    vector<simulData>   myDataline;

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

        myTimeline.push_back(myMaturity);

        const size_t n = myTimeline.size();
        myDataline.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            myDataline[i].forwardMats.push_back(myTimeline[i]);
        }
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

    //  Dataline
    const vector<simulData>& dataline() const override
    {
        return myDataline;
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
        const double smooth = convert<double>(path[0].forwards[0] * 0.01);

        //  We start alive
        T alive(1.0);

        //  Go through path, update alive status
        for (const auto& scen: path)
        {
            //  Breached
            if (scen.forwards[0] > myBarrier + smooth)
            {
                payoffs[0] = T(0.0);
                return;
            }

            //  Semi-breached: apply smoothing
            if (scen.forwards[0] > myBarrier - smooth)
            {
                alive *= (myBarrier + smooth - scen.forwards[0]) / (2 * smooth);
            }
        }

        //  Payoff
        payoffs[0] = alive * max(path.back().forwards[0] - myStrike, 0.0)
            / path.back().numeraire;
    }
};

template <class T>
class Europeans : public Product<T>
{
    vector<Time>            myMaturities;   //  = timeline
    vector<vector<double>>  myStrikes;      //  a vector of strikes per maturity
    vector<simulData>       myDataline;

public:

    //  Constructor: store data and build timeline
    Europeans(const map<Time, vector<double>>& options) 
    {
        const size_t n = options.size();

        for (const pair<Time, vector<double>>& p : options)
        {
            myMaturities.push_back(p.first);
            myStrikes.push_back(p.second);
        }

        myDataline.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            myDataline[i].forwardMats.push_back(myMaturities[i]);
        }
    }

    //  access to maturities and strikes
    const vector<Time>& maturities() const
    {
        return myMaturities;
    }

    const vector<vector<double>>& strikes() const
    {
        return myStrikes;
    }

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

    //  Dataline
    const vector<simulData>& dataline() const override
    {
        return myDataline;
    }

    size_t numPayoffs() const override
    {
        return accumulate(myStrikes.begin(), myStrikes.end(), 0, 
            [](const int acc, const vector<double>& strikes) 
            { 
                return acc + strikes.size(); 
            }
        );
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
        for (size_t i = 0; i < numT; ++i)
        {
            transform(
                    myStrikes[i].begin(),
                    myStrikes[i].end(),
                    payoffIt,
                    [spot = path[i].forwards[0], num = path[i].numeraire] (const double& k) 
                    {
                        return max(spot - k, 0.0) / num; 
                    }
            );

            payoffIt += myStrikes[i].size();
        }
    }
};
