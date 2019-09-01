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