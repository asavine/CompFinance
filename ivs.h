#pragma once

#include "mcBase.h"
#include "matrix.h"
#include "gaussians.h"

//  Risk view
template <class T>
class RiskView
{
    bool myEmpty;

    vector<double> myStrikes;
    vector<Time> myMats;
    matrix<T> mySpreads;

public:

    //  Default constructor, empty view
    RiskView() : myEmpty(true) {}

    //  Intializes risk view AND put on tape
    RiskView(const vector<double>& strikes, const vector<Time>& mats) :
        myEmpty(false), myStrikes(strikes), myMats(mats), mySpreads(strikes.size(), mats.size())
    {
        for (auto& spr : mySpreads) spr = convert<T>(0.0);
    }

    //  Get spread
    T spread(const double strike, const Time mat) const
    {
        return myEmpty
            ? convert<T>(0.0) 
            : interp2D(myStrikes, myMats, mySpreads, strike, mat, true);
    }

    //  Accessors by const ref
    bool empty() const { return myEmpty; }
    size_t rows() const { return myStrikes.size(); }
    size_t cols() const { return myMats.size(); }
    const vector<double>& strikes() const { return myStrikes; }
    const vector<Time>& mats() const { return myMats; }
    const matrix<T>& risks() const { return mySpreads; }

    //  Iterators
    typedef typename matrix<T>::iterator iterator;
    typedef typename matrix<T>::const_iterator const_iterator;
    iterator begin() { return mySpreads.begin(); }
    iterator end() { return mySpreads.end(); }
    const_iterator begin() const { return mySpreads.begin(); }
    const_iterator end() const { return mySpreads.end(); }
};

class IVS
{
    //  To avoid refer a linear market
    double mySpot;

public:

    IVS(const double spot) : mySpot(spot) {}

    //  Read access to spot
    double spot() const
    {
        return mySpot;
    }

    //  Raw implied vol
    virtual double impliedVol(const double strike, const Time mat) const = 0;

    //  Call price
    template<class T = double>
    T call(
        const double strike, 
        const Time mat, 
        const RiskView<T>& risk = RiskView<double>()) const
    {
        return blackScholes(
            mySpot,
            strike,
            impliedVol(strike, mat) + risk.spread(strike, mat),
            mat);
    }

    //  Local vol, dupire's formula
    template<class T = double>
    T localVol(
        const double strike,
        const double mat,
        const RiskView<T>& risk = RiskView<double>()) const
    {
        //  Derivative to time
        const T c00 = call(strike, mat, risk);
        const T c01 = call(strike, mat - 1.0e-04, risk);
        const T c02 = call(strike, mat + 1.0e-04, risk);
        const T ct = (c02 - c01) * 0.5e04;

        //  Second derivative to strike = density
        const T c10 = call(strike - 1.0e-04, mat, risk);
        const T c20 = call(strike + 1.0e-04, mat, risk);
        const T ckk = (c10 + c20 - 2.0 * c00) * 1.0e08;
        
        //  Dupire's formula
        return sqrt(2.0 * ct / ckk);
    }

    virtual ~IVS() {}
};

//  Concrete IVS just override (raw) call prices

class BachelierIVS : public IVS
{
    double myBachVol;

public:

    BachelierIVS(const double spot, const double logVol)
        : IVS(spot), myBachVol(logVol * spot) {}

    double impliedVol(const double strike, const Time mat) const override
    {
        const double call = bachelier(spot(), strike, myBachVol, mat);
        return blackScholesIvol(spot(), strike, call, mat);
    }
};

class BlackScholesIVS : public IVS
{
    double myVol;

public:

    BlackScholesIVS(const double spot, const double vol)
        : IVS(spot), myVol(vol) {}

    double impliedVol(const double strike, const Time mat) const override
    {
        return myVol;
    }
};

class MertonIVS : public IVS
{
    double myVol;
    double myIntensity, myAverageJmp, myJmpStd;

public:

    MertonIVS(const double spot, const double vol, 
        const double intens, const double aveJmp, const double stdJmp)
        : IVS(spot), 
        myVol(vol),
        myIntensity(intens),
        myAverageJmp(aveJmp),
        myJmpStd(stdJmp)
    {}

    double impliedVol(const double strike, const Time mat) const override
    {
        const double call
            = merton(
                spot(),
                strike,
                myVol,
                mat,
                myIntensity,
                myAverageJmp,
                myJmpStd);

        return blackScholesIvol(spot(), strike, call, mat);
    }
};
