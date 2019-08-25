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
	vector<string>          myLabels;

public:

	//  Constructor: store data and build timeline
	TestMulti(const vector<string>& assets, const Time T1, const Time T2, const double K0, const double )
	{
		const size_t n = options.size();

		//  Timeline = each maturity is an event date
		for (const pair<Time, vector<double>>& p : options)
		{
			myMaturities.push_back(p.first);
			myStrikes.push_back(p.second);
		}

		//  Defline = num and spot(t) = forward(t,t) on every step
		myDefline.resize(n);
		for (size_t i = 0; i < n; ++i)
		{
			myDefline[i].numeraire = true;
            myDefline[i].forwardMats.push_back({ myMaturities[i] });
		}

		//  Identify the payoffs
		for (const auto& option : options)
		{
			for (const auto& strike : option.second)
			{
				ostringstream ost;
				ost.precision(2);
				ost << fixed;
				ost << "call " << option.first << " " << strike;
				myLabels.push_back(ost.str());
			}
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
		return make_unique<Europeans<T>>(*this);
	}

	//  Timeline
	const vector<Time>& timeline() const override
	{
		return myMaturities;
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
		const size_t numT = myMaturities.size();

		auto payoffIt = payoffs.begin();
		for (size_t i = 0; i < numT; ++i)
		{
			transform(
				myStrikes[i].begin(),
				myStrikes[i].end(),
				payoffIt,
				[spot = path[i].forwards.front().front(), num = path[i].numeraire]
				(const double& k)
				{
					return max(spot - k, 0.0) / num;
				}
			);

			payoffIt += myStrikes[i].size();
		}
	}
};