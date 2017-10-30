#pragma once

#include "mcBase.h"
#include "matrix.h"
#include "gaussians.h"


template<class T>
class IVS
{
    //  To avoid refer a linear market
    T mySpot;

    //  Raw call prices
    virtual double callRaw(const double strike, const Time mat) const = 0;

    //  Risk grid
    vector<double> myRiskStrikes;
    vector<Time> myRiskMats;
    matrix<T> myRiskSpreads;
    //  ij = applies to strikes between Ki-1 and Ki
    //  and mats between Tj-1 and Tj

    T riskSpread(const double strike, const Time mat)
    {
        size_t i0;
        auto strIt = lower_bound(myRiskStrikes.begin(), myRiskStrikes.end(), strike);
        if (strIt == myRiskStrikes.end())
        {
            i0 = myRiskStrikes.size();
        }
        else
        {
            i0 = distance(myRiskStrikes.begin(), strIt);
        }

        size_t j0;
        auto matIt = lower_bound(myRiskMats.begin(), myRiskMats.end(), mat);
        if (matIt == myRiskMats.end())
        {
            j0 = myRiskMats.size();
        }
        else
        {
            j0 = distance(myRiskMats.begin(), matIt);
        }

        return myRiskSpreads(i0, j0);
    }

public:

    IVS(const T spot) : mySpot(spot) {}

    T spot() const
    {
        return mySpot;
    }

    const tuple<const vector<double>&, const vector<Time>&, const matrix<double>&>
        riskSpreads() const
    {
        return make_tuple(myRiskStrikes, myRiskMats, myRiskSpreads);
    }

    T call(const double strike, const double mat) const
    {
        return convert<T>(callraw(strike, mat)) + riskSpread(strike, mat);
    }

    //  dC / dT
    T cT(const double strike, const double mat) const
    {
        const double t1 = max(0.0, mat - 0.002739726);   //  1 day
        const double t2 = 0.0, mat + 0.002739726;
        const T c1 = call(strike, t1), c2 = call(strike, t2);
        return (c2 - c1) / (t2 - t1);
    }

    //  d2C / dK2
    T cKK(const double strike, const double mat) const
    {
        const double k0 = strike;
        const double k1 = strike - 1.e-04;
        const double k2 = strike + 1.e-04;

        const T c0 = call(k0 , mat), c1 = call(k1, mat), c2 = call(k2, mat);

        const double dsu = k2 - k0;
        const double dsd = k0 - k1;

        const T dcu = c2 - c0;
        const T dcd = c0 - c1;

        return 2.0 * (dsd * dcu - dsu * dcd) / (dsu * dsd * (dsu + dsd));
    }

    void setRiskGrid(const vector<double>& strikes, const vector<Time>& mats)
    {
        myRiskStrikes.resize(max(1, strikes.size() - 1));
        myRiskMats.resize(max(1, mats.size() - 1));

        if (strikes.size() == 1)
        {
            myRiskStrikes[0] = strikes[0];
        }
        else
        {
            for (size_t i = 0; i < strikes.size() - 1; ++i)
            {
                myRiskStrikes[i] = 0.5 * (strikes[i] + strikes[i + 1]);
            }
        }

        if (mats.size() == 1)
        {
            myRiskMats[0] = mats[0];
        }
        else
        {
            for (size_t i = 0; i < mats.size() - 1; ++i)
            {
                myRiskMats[i] = 0.5 * (mats[i] + mats[i + 1]);
            }
        }

        myRiskSpreads.resize(strikes.size(), mats.size());
        for (auto& spr : myRiskSpreads) spr = convert<T>(0.0);
    }

    virtual ~IVS() {}
};

template <class T>
class BachelierIVS : public IVS<T>
{
    double mySpot, myVol;

    double callRaw(const double strike, const Time mat) const override
    {
        return bachelier(mySpot, strike, myVol, mat);
    }

public:

    BachelierIVS(const T spot, const T vol)
        : IVS(spot), mySpot(convert<double>(spot)), myVol(convert<double>(vol)) {}
};

template <class T>
class BlackScholesIVS : public IVS<T>
{
    double mySpot, myVol;

    double callRaw(const double strike, const Time mat) const override
    {
        return blackScholes(mySpot, strike, myVol, mat);
    }

public:

    BlackScholesIVS(const T spot, const T vol)
        : IVS(spot), mySpot(convert<double>(spot)), myVol(convert<double>(vol)) {}
};

template <class T>
class MertonIVS : public IVS<T>
{
    double mySpot, myVol;
    double myIntensity, myAverageJmp, myJmpStd;

    double callRaw(const double strike, const Time mat) const override
    {
        return merton(mySpot, strike, myVol, mat, myIntensity, myAverageJmp, myJmpStd);
    }

public:

    MertonIVS(const T spot, const T vol, const T intens, const T aveJmp, const T stdJmp)
        : IVS(spot), 
        mySpot(convert<double>(spot)), 
        myVol(convert<double>(vol)),
        myVol(convert<double>(intens)),
        myVol(convert<double>(aveJmp)),
        myVol(convert<double>(stdJmp)),
    {}
};

