#pragma once

#include "mcPrd.h"

/*	A test class, computing statistics on multiple assets on multiple dates
	On each date T, the payoffs are xi and xixj so their expectations give expectations and covariances
	We also export yi and yiyj where y = xi(Tk) - xi(Tk-1) so we get stats on increments too
	Finally, xi(Tk) is the fwd fixed at Tk for the corresponding maturity Tk' >= Tk
*/

template <class T>
class MultiStats : public Product<T>
{
	//	Timeline
	vector<Time>			myFixDates;
	vector<Time>			myFwdDates;

	//	Assets
	size_t					myNumAssets;
	vector<string>			myAssetNames;

	//	Defline and labels
	vector<SampleDef>       myDefline;

    size_t                  myNumPayoffs;
	vector<string>          myLabels;

public:

	//  Constructor: store data and build timeline
	MultiStats(const vector<string>& assets, const vector<Time>& fixDates, const vector<Time>& fwdDates) :
        myNumAssets(assets.size()), myAssetNames(assets), myFixDates(fixDates), myFwdDates(fwdDates)
	{
		//  Defline = num and forward(Tfix, Tfwd) on every fix date, for all assets
		const size_t nTimes = fixDates.size();
		myDefline.resize(nTimes);
		for (size_t i = 0; i < nTimes; ++i)
		{
			myDefline[i].numeraire = false;
            myDefline[i].forwardMats.resize(myNumAssets);
            fill(myDefline[i].forwardMats.begin(), myDefline[i].forwardMats.end(), vector<Time>(1, myFwdDates[i]));
		}

		//  Identify the payoffs
        myNumPayoffs = (nTimes + (nTimes > 1 ? nTimes - 1: 0)) * (myNumAssets + myNumAssets * (myNumAssets + 3) / 2);
        
        //  First, fwd fixings on fix dates
        for (size_t t = 0; t < nTimes; ++t)
        {
            for (size_t a1 = 0; a1 < myNumAssets; ++a1)
            {
			    ostringstream ost;
			    ost.precision(2);
			    ost << fixed;
			    ost << myAssetNames[a1] << " " << myFixDates[t] << " " << myFwdDates[t];
			    myLabels.push_back(ost.str());            
            }

            for (size_t a1 = 0; a1 < myNumAssets; ++a1) for (size_t a2 = 0; a2 <= a1; ++a2)
            {
			    ostringstream ost;
			    ost.precision(2);
			    ost << fixed;
			    ost << myAssetNames[a1] << " " << myAssetNames[a2] << " " << myFixDates[t] << " " << myFwdDates[t];
			    myLabels.push_back(ost.str());                        
            }
        }

        //  Next, differences
        for (size_t t2 = 1; t2 < nTimes; ++t2)
        {
            const size_t t1 = t2 - 1;
            for (size_t a1 = 0; a1 < myNumAssets; ++a1)
            {
			    ostringstream ost;
			    ost.precision(2);
			    ost << fixed;
			    ost << myAssetNames[a1] << " " << myFixDates[t1] << " " << myFwdDates[t1] << " - " << myFixDates[t2] << " " << myFwdDates[t2];
			    myLabels.push_back(ost.str());            
            }

            for (size_t a1 = 0; a1 < myNumAssets; ++a1) for (size_t a2 = 0; a2 <= a1; ++a2)
            {
			    ostringstream ost;
			    ost.precision(2);
			    ost << fixed;
			    ost << myAssetNames[a1] << " " << myAssetNames[a2] << " " << myFixDates[t1] << " " << myFwdDates[t1] << " - " << myFixDates[t2] << " " << myFwdDates[t2];
			    myLabels.push_back(ost.str());                        
            }
        }
	}

	//  Read access to parameters
	const size_t numAssets() const override
	{
		return myNumAssets;
	}

    const vector<string>& assetNames() const override
    {
        return myAssetNames;
    }

    const vector<Time>& fixDates() const 
    {
        return myFixDates;
    }

    const vector<Time>& fwdDates() const 
    {
        return myFwdDates;
    }

	//  Virtual copy constructor
	unique_ptr<Product<T>> clone() const override
	{
		return make_unique<MultiStats<T>>(*this);
	}

	//  Timeline
	const vector<Time>& timeline() const override
	{
		return myFixDates;
	}

	//  Defline
	const vector<SampleDef>& defline() const override
	{
		return myDefline;
	}

	//  Labels
	const vector<string>& payoffLabels() const override
	{
		return myLabels;
	}

	//  Evaluation on scenario
	void payoffs(
		//  path, one entry per time step 
		const Scenario<T>&          path,
		//  pre-allocated space for resulting payoffs
		vector<T>&                  payoffs)
		const override
	{
        const size_t nTimes = myFixDates.size();
        size_t payIdx = 0;

        //  First, fwd fixings on fix dates
        for (size_t t = 0; t < nTimes; ++t)
        {
            for (size_t a1 = 0; a1 < myNumAssets; ++a1)
            {
                payoffs[payIdx++] = path[t].forwards[a1].front();
            }

            for (size_t a1 = 0; a1 < myNumAssets; ++a1) for (size_t a2 = 0; a2 <= a1; ++a2)
            {
                payoffs[payIdx++] = path[t].forwards[a1].front() * path[t].forwards[a2].front();
            }
        }

        //  Next, differences
        for (size_t t2 = 1; t2 < nTimes; ++t2)
        {
            const size_t t1 = t2 - 1;
            for (size_t a1 = 0; a1 < myNumAssets; ++a1)
            {
                payoffs[payIdx++] = path[t2].forwards[a1].front() - path[t1].forwards[a1].front();
            }

            for (size_t a1 = 0; a1 < myNumAssets; ++a1) for (size_t a2 = 0; a2 <= a1; ++a2)
            {
                payoffs[payIdx++] = (path[t2].forwards[a1].front() - path[t1].forwards[a1].front())
                                        * (path[t2].forwards[a2].front() - path[t1].forwards[a2].front());
            }
        }
	}
};

template <class T>
class Baskets : public Product<T>
{
	//	Assets and weights
	size_t					myNumAssets;
	vector<string>			myAssetNames;
	vector<double>			myWeights;

	//	Timeline
	Time					myMaturity;
	//  Vector of strikes 
	vector<double>			myStrikes;
	
	vector<Time>			myTimeline;
	vector<SampleDef>       myDefline;
	vector<string>          myLabels;

public:

	//  Constructor: store data and build timeline
	Baskets(const vector<string>& assets, const vector<double> weights, const Time maturity, const vector<double>& strikes) :
		myNumAssets(assets.size()), myAssetNames(assets), myWeights(weights), myMaturity(maturity), myStrikes(strikes),
		myTimeline(1, maturity), myDefline(1)
	{
		const size_t n = strikes.size();

		//  Defline = num and spot(t) = forward(t,t) on maturity for all assets
		myDefline[0].numeraire = true;
		myDefline[0].forwardMats = vector<vector<Time>>(myNumAssets, { maturity });

		//  Identify the payoffs
		myLabels.reserve(strikes.size());
		for (const double strike : strikes)
		{
			ostringstream ost;
			ost.precision(2);
			ost << fixed;
			ost << "basket strike " << strike;
			myLabels.push_back(ost.str());
		}
	}

    const size_t numAssets() const override
	{
		return myNumAssets;
	}

    const vector<string>& assetNames() const override
    {
        return myAssetNames;
    }

	//  access to weights, maturity and  strikes
	const vector<double>& weights() const
	{
		return myWeights;
	}

	Time maturity() const
	{
		return myMaturity;
	}

	const vector<double>& strikes() const
	{
		return myStrikes;
	}

	//  Virtual copy constructor
	unique_ptr<Product<T>> clone() const override
	{
		return make_unique<Baskets<T>>(*this);
	}

	//  Timeline
	const vector<Time>& timeline() const override
	{
		return myTimeline;
	}

	//  Defline
	const vector<SampleDef>& defline() const override
	{
		return myDefline;
	}

	//  Labels
	const vector<string>& payoffLabels() const override
	{
		return myLabels;
	}

	//  Payoffs, maturity major
	void payoffs(
		//  path, one entry per time step 
		const Scenario<T>&          path,
		//  pre-allocated space for resulting payoffs
		vector<T>&                  payoffs)
		const override
	{
		T basket = inner_product(myWeights.begin(), myWeights.end(), path[0].forwards.begin(), T(0.0),
			plus<T>(), [](const double weight, const vector<T>& fwds) { return weight * fwds[0]; });

		transform(myStrikes.begin(), myStrikes.end(), payoffs.begin(),
			[&basket, num = path[0].numeraire](const double k) {return max(basket - k, 0.0) / num; });
	}
};

template <class T>
class Autocall : public Product<T>
{
	//	Assets and weights
	size_t					myNumAssets;
	vector<string>			myAssetNames;

	//	Timeline
	Time					myMaturity;
    int                     myNumPeriods;
	
	//	References
	vector<double>			myRefs;

    //  Barriers
	double                  myKO;
    double                  myStrike;
    double                  myCpn;
	
    double                  mySmooth;

	vector<Time>			myTimeline;
	vector<SampleDef>       myDefline;
	vector<string>          myLabels;

public:

	//  Constructor: store data and build timeline
	Autocall(const vector<string>& assets, const vector<double> refs, const Time maturity, const int periods, const double ko, const double strike, const double cpn, const double smooth) :
		myNumAssets(assets.size()), myAssetNames(assets), myRefs(refs), myMaturity(maturity), myNumPeriods(periods), myKO(ko), myStrike(strike), myCpn(cpn), mySmooth(max(smooth, EPS)),
        myTimeline(periods), myDefline(periods), myLabels(1)
	{
        //  Timeline and defline
        
        //  Periods
        Time time = systemTime;
        const double dt = maturity / periods;
        for (size_t step=0; step<periods; ++step)
        {
            time += dt;
            myTimeline[step] = time;
            myDefline[step].numeraire = true;
		    myDefline[step].forwardMats = vector<vector<Time>>(myNumAssets, { time });  //  spot(t) = forward(t,t) on maturity for all assets
        }

		//  Identify the payoff
		myLabels[0] = "autocall strike " + to_string(int(100 * myStrike + EPS)) 
            + " KO " +  to_string(int(100 * myKO + EPS)) 
            + " CPN " + to_string(int(100 * myCpn + EPS)) + " "
            + to_string(periods) + " periods of " + to_string(int(12 * maturity / periods + EPS)) + "m";
	}

    const size_t numAssets() const override
	{
		return myNumAssets;
	}

    const vector<string>& assetNames() const override
    {
        return myAssetNames;
    }

	const vector<double>& refs() const
	{
		return myRefs;
	}

	//  access to maturity, strike, KO and cpn

	Time maturity() const { return myMaturity; }
	double strike() const { return myStrike; }
	double ko() const { return myKO; }
	double cpn() const { return myCpn; }

	//  Virtual copy constructor
	unique_ptr<Product<T>> clone() const override
	{
		return make_unique<Autocall<T>>(*this);
	}

	//  Timeline
	const vector<Time>& timeline() const override
	{
		return myTimeline;
	}

	//  Defline
	const vector<SampleDef>& defline() const override
	{
		return myDefline;
	}

	//  Labels
	const vector<string>& payoffLabels() const override
	{
		return myLabels;
	}

	//  Payoff
	void payoffs(
		//  path, one entry per time step 
		const Scenario<T>&          path,
		//  pre-allocated space for resulting payoffs
		vector<T>&                  payoffs)
		const override
	{
        //  Temporaries
        static thread_local vector<T> perfs;
        perfs.resize(myNumAssets);

        //  Periods
        const double dt = myMaturity / myNumPeriods;
        T notionalAlive(1.0);
        payoffs[0] = 0.0;
        for (int step=0; step<myNumPeriods-1; ++step)
        {
            auto& state = path[step];
            transform(state.forwards.begin(), state.forwards.end(), myRefs.begin(), perfs.begin(), [](const vector<T>& fwds, const double ref) { return fwds[0] / ref; });
            T worst = *min_element(perfs.begin(), perfs.end());

            //  receive cpn
            payoffs[0] += notionalAlive * myCpn * dt / state.numeraire;

            //  apply ko smoothly
            T notionalSurviving = notionalAlive * min(1.0, max(0.0, (myKO + mySmooth - worst) / 2 / mySmooth));
            T notionalDead = notionalAlive - notionalSurviving;

            //  receive redemption on dead notional
            payoffs[0] += notionalDead / state.numeraire;

            //  continue with the rest
            notionalAlive = notionalSurviving;
        }

        //  last
        {
            int step = myNumPeriods-1;
            auto& state = path[step];
            transform(state.forwards.begin(), state.forwards.end(), myRefs.begin(), perfs.begin(), [](const vector<T>& fwds, const double ref) { return fwds[0] / ref; });
            T worst = *min_element(perfs.begin(), perfs.end());

            //  receive cpn
            payoffs[0] += notionalAlive * myCpn * dt / state.numeraire;

            //  receive redemption
            payoffs[0] += notionalAlive / state.numeraire;

            // pay put
            payoffs[0] -= notionalAlive * max(myStrike - worst, 0.0) / myStrike / state.numeraire;        
        }
	}
};