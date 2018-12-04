#pragma once

#include "mcBase.h"
#include "interp.h"
#include "gaussians.h"
#include <iostream>

struct Record
{
    int     numArg;       //  ToyNumber of arguments: 0, 1 or 2
    int     idx1;         //  index of first argument on tape
    int     idx2;         //  index of second argument on tape
    double  der1;         //  partial derivative to first argument
    double  der2;         //  partial derivative to second argument
};

//  The tape, declared as a global variable
vector<Record> tape;

//  Custom ToyNumber type
struct ToyNumber
{
    double  value;
    int     idx;

    //  default constructor does nothing
    ToyNumber() {}

    //  constructs with a value and record
    ToyNumber(const double& x) : value(x)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  reference record on tape
        idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 0;
    }

    ToyNumber operator +() const { return *this; }
    ToyNumber operator -() const { return ToyNumber(0.0) - *this; }

    ToyNumber& operator +=(const ToyNumber& rhs) { *this = *this + rhs; return *this; }
    ToyNumber& operator -=(const ToyNumber& rhs) { *this = *this - rhs; return *this; }
    ToyNumber& operator *=(const ToyNumber& rhs) { *this = *this * rhs; return *this; }
    ToyNumber& operator /=(const ToyNumber& rhs) { *this = *this / rhs; return *this; }

    friend ToyNumber operator+(const ToyNumber& lhs, const ToyNumber& rhs)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  compute result
        ToyNumber result;
        result.value = lhs.value + rhs.value;

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 2;
        rec.idx1 = lhs.idx;
        rec.idx2 = rhs.idx;

        //  compute derivatives
        rec.der1 = 1;
        rec.der2 = 1;

        return result;
    }

    friend ToyNumber operator-(const ToyNumber& lhs, const ToyNumber& rhs)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  compute result
        ToyNumber result;
        result.value = lhs.value - rhs.value;

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 2;
        rec.idx1 = lhs.idx;
        rec.idx2 = rhs.idx;

        //  compute derivatives -
        rec.der1 = 1;
        rec.der2 = -1;

        return result;
    }

    friend ToyNumber operator*(const ToyNumber& lhs, const ToyNumber& rhs)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  compute result
        ToyNumber result;
        result.value = lhs.value * rhs.value;

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 2;
        rec.idx1 = lhs.idx;
        rec.idx2 = rhs.idx;

        //  compute derivatives *
        rec.der1 = rhs.value;
        rec.der2 = lhs.value;

        return result;
    }

    friend ToyNumber operator/(const ToyNumber& lhs, const ToyNumber& rhs)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  compute result
        ToyNumber result;
        result.value = lhs.value / rhs.value;

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 2;
        rec.idx1 = lhs.idx;
        rec.idx2 = rhs.idx;

        //  compute derivatives /
        rec.der1 = 1.0 / rhs.value;
        rec.der2 = -lhs.value / (rhs.value * rhs.value);

        return result;
    }

    friend ToyNumber log(const ToyNumber& arg)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  compute result
        ToyNumber result;
        result.value = log(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 1;
        rec.idx1 = arg.idx;

        //  compute derivative
        rec.der1 = 1.0 / arg.value;

        return result;
    }

    friend ToyNumber exp(const ToyNumber& arg)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  compute result
        ToyNumber result;
        result.value = exp(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 1;
        rec.idx1 = arg.idx;

        //  compute derivative
        rec.der1 = result.value;

        return result;
    }

    friend ToyNumber sqrt(const ToyNumber& arg)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  compute result
        ToyNumber result;
        result.value = sqrt(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 1;
        rec.idx1 = arg.idx;

        //  compute derivative
        rec.der1 = 0.5 / result.value;

        return result;
    }

    friend ToyNumber normalDens(const ToyNumber& arg)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  compute result
        ToyNumber result;
        result.value = normalDens(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 1;
        rec.idx1 = arg.idx;

        //  compute derivative
        rec.der1 = -result.value * arg.value;

        return result;
    }

    friend ToyNumber normalCdf(const ToyNumber& arg)
    {
        //  create a new record on tape
        tape.push_back(Record());
        Record& rec = tape.back();

        //  compute result
        ToyNumber result;
        result.value = normalCdf(arg.value);

        //  reference record on tape
        result.idx = tape.size() - 1;

        //  populate record on tape
        rec.numArg = 1;
        rec.idx1 = arg.idx;

        //  compute derivative
        rec.der1 = normalDens(arg.value);

        return result;
    }

    friend bool operator==(const ToyNumber& lhs, const ToyNumber& rhs) { return lhs.value == rhs.value; }
    friend bool operator!=(const ToyNumber& lhs, const ToyNumber& rhs) { return lhs.value != rhs.value; }
    friend bool operator>(const ToyNumber& lhs, const ToyNumber& rhs) { return lhs.value > rhs.value; }
    friend bool operator>=(const ToyNumber& lhs, const ToyNumber& rhs) { return lhs.value >= rhs.value; }
    friend bool operator<(const ToyNumber& lhs, const ToyNumber& rhs) { return lhs.value < rhs.value; }
    friend bool operator<=(const ToyNumber& lhs, const ToyNumber& rhs) { return lhs.value <= rhs.value; }

};

template <class T> inline T blackScholes(
    //  input layer 0
    const T spot, const T rate, const T yield, const T vol, const T strike, const T mat)
{
/* layer 1 */        T df = exp(-rate * mat), fwd = spot * exp((rate - yield) * mat), std = vol * sqrt(mat);
/* layer 2 */        T d = log(fwd / strike) / std;
/* layer 3 */        T d1 = d + 0.5 * std, d2 = d - 0.5 * std;
/* layer 4 */        T p1 = normalCdf(d1), p2 = normalCdf(d2);
/* output layer 5 */ return df * (fwd * p1 - strike * p2);
}

vector<double> calculateAdjoints(ToyNumber& result)
{
    //  initialization
    vector<double> adjoints(tape.size(), 0.0);  //  initialize all to 0
    int N = result.idx;                         //  find N
    adjoints[N] = 1.0;                          //  seed aN = 1
    
    //  backward propagation
    for(int j=N; j>0; --j)  //  iterate backwards over tape
    {
        if (tape[j].numArg > 0)
        {
            adjoints[tape[j].idx1] += adjoints[j] * tape[j].der1;       //  propagate first argument

            if (tape[j].numArg > 1)
            {
                adjoints[tape[j].idx2] += adjoints[j] * tape[j].der2;   //  propagate second argument
            }
        }
    }

    return adjoints;
}

void blackScholesDiff()
{
    ToyNumber spot = 100, rate = 0.02, yield = 0.05, vol = 0.2, strike = 110, mat = 2; // initializes and records inputs
    auto result = blackScholes(spot, rate, yield, vol, strike, mat);                // evaluates and records operations
    cout << "Value = " << result.value << endl;   //  5.03705

    //  propagate adjoints
    vector<double> adjoints = calculateAdjoints(result);

    //  show derivatives
    cout << "Derivative to spot (delta) = " << adjoints[spot.idx] << endl;          //  0.309
    cout << "Derivative to rate (rho) = " << adjoints[rate.idx] << endl;            //  51.772
    cout << "Derivative to dividend yield = " << adjoints[yield.idx] << endl;       //  -61.846
    cout << "Derivative to volatility (vega) = " << adjoints[vol.idx] << endl;      //  46.980
    cout << "Derivative to strike (-digital) = " << adjoints[strike.idx] << endl;   //  -0.235
    cout << "Derivative to maturity (-theta) = " << adjoints[mat.idx] << endl;      //  1.321

    while (true);
}

template <class T>
inline T toyDupireBarrierMc(
    //  Spot
    const T            S0,
    //  Local volatility
    const vector<T>    spots,
    const vector<T>    times,
    const matrix<T>    vols,
    //  Product parameters
    const T            maturity,
    const T            strike,
    const T            barrier,
    //  Number of paths and time steps
    const int          Np,
    const int          Nt,
    //  Smoothing
    const T            epsilon,
    //  Initialized random number generator
    RNG&               random)
{
    //  Initialize
    T result = 0;
    vector<double> gaussianIncrements(Nt); // double because the RNG is not templated (and doesn't need to be, see chapter 12)
    const T dt = maturity / Nt, sdt = sqrt(dt);

    //  Loop over paths
    for (int i = 0; i < Np; ++i)
    {
        //  Generate Nt Gaussian Numbers
        random.nextG(gaussianIncrements);
        //  Step by step
        T spot = S0, time = 0;
        /*  bool alive = true; */ T alive = 1.0; // alive is a real number in (0,1)
        for (size_t j = 0; j < Nt; ++j)
        {
            //  Interpolate volatility
            const T vol = interp2D(spots, times, vols, spot, time);
            time += dt;
            //  Simulate return
            spot *= exp(-0.5 * vol * vol * dt + vol * sdt * gaussianIncrements[j]);
            //  Monitor barrier
            /* if (spot > barrier) { alive = false; break; } */
            if (spot > barrier + epsilon) { alive = 0.0; break; }       //   definitely dead
            else if (spot < barrier - epsilon) { /* do nothing */ }     //   definitely alive
            else /* in between, interpolate */ alive *= 1.0 - (spot - barrier + epsilon) / (2 * epsilon);

        }
        //  Payoff
        /* if (alive && spot > strike) result += spot - strike; */ if (spot > strike) result += alive * (spot - strike); // pay on surviving notional
    }   //  paths

    return result / Np;
}

void dupireRisksMiniBatch(
	const double S0, const vector<double> spots, const vector<double> times, const matrix<double> vols,
	const double maturity, const double strike, const double barrier,
	const int Np, const int Nt, const double epsilon, RNG& random,
	/* results: value and dV/dS, dV/d(local vols) */ double& price, double& delta, matrix<double>& vegas)
{

	//	1. Initialize inputs, record on tape

	ToyNumber nS0(S0), nMaturity(maturity), nStrike(strike), nBarrier(barrier), nEpsilon(epsilon);
	vector<ToyNumber> nSpots(spots.size()), nTimes(times.size());
	matrix<ToyNumber> nVols(vols.rows(), vols.cols());

	for (int i = 0; i < spots.size(); ++i) nSpots[i] = ToyNumber(spots[i]);
	for (int i = 0; i < times.size(); ++i) nTimes[i] = ToyNumber(times[i]);
	for (int i = 0; i < vols.rows(); ++i) for (int j = 0; j < vols.cols(); ++j) nVols[i][j] = ToyNumber(vols[i][j]);

	//	2. Call instrumented evaluation code, which evaluates the barrier option price and records all operations

	ToyNumber nPrice = toyDupireBarrierMc(nS0, nSpots, nTimes, nVols, nMaturity, nStrike, nBarrier, Np, Nt, nEpsilon, random);

	//	3. Adjoint propagation

    //  propagate adjoints
    vector<double> adjoints = calculateAdjoints(nPrice);

	//	4. Pick results

	price = nPrice.value;
	delta = adjoints[nS0.idx];
	for (int i = 0; i < vols.rows(); ++i) for (int j = 0; j < vols.cols(); ++j) vegas[i][j] = adjoints[nVols[i][j].idx];
}

void toyDupireBarrierMcRisks(
	const double S0, const vector<double> spots, const vector<double> times, const matrix<double> vols,
	const double maturity, const double strike, const double barrier,
	const int Np, const int Nt, const double epsilon, RNG& random,
	/* results: value and dV/dS, dV/d(local vols) */ double& price, double& delta, matrix<double>& vegas)
{
	price = 0;
	delta = 0;
	for (int i = 0; i < vegas.rows(); ++i) for (int j = 0; j < vegas.cols(); ++j) vegas[i][j] = 0;

	double batchPrice, batchDelta;
	matrix<double> batchVegas(vegas.rows(), vegas.cols());

	int pathsToGo = Np, pathsPerBatch = 1024;
	while (pathsToGo > 0)
	{
		//	wipe tape
		tape.clear();

		//	do mini batch
		int paths = min(pathsToGo, pathsPerBatch);
		dupireRisksMiniBatch(S0, spots, times, vols, maturity, strike, barrier, paths, Nt, epsilon, random, batchPrice, batchDelta, batchVegas);

		//	update results
		price += batchPrice * paths / Np;
		delta += batchDelta * paths / Np;
		for (int i = 0; i < vegas.rows(); ++i) for (int j = 0; j < vegas.cols(); ++j) 
			vegas[i][j] += batchVegas[i][j] * paths / Np;

		pathsToGo -= paths;
	}
}