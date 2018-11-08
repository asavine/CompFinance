
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

//  Excel export wrappers to functions in main.h

#pragma warning(disable:4996)

#include "ThreadPool.h"
#include "main.h"
#include "toyCode.h"

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

#include "xlcall.h"
#include "xlframework.h"
#include "xlOper.h"

//  Helpers
NumericalParam xl2num(
    const double              useSobol,
    const double              seed1,
    const double              seed2,
    const double              numPath,
    const double              parallel)
{
    NumericalParam num;

    num.numPath = static_cast<int>(numPath + EPS);
    num.parallel = parallel > EPS;
	if (seed1 >= 1)
	{
		num.seed1 = static_cast<int>(seed1 + EPS);
	}
	else
	{
		num.seed1 = 1234;
	}

	if (seed2 >= 1)
	{
		num.seed2 = static_cast<int>(seed2 + EPS);
	}
	else
	{
		num.seed2 = num.seed1 + 1;
	}

	num.useSobol = useSobol > EPS;

    return num;
}

//	Wrappers

//  change number of threads in the pool
extern "C" __declspec(dllexport)
double xRestartThreadPool(
    double              xNthread)
{
    const int numThread = int(xNthread + EPS);
    ThreadPool::getInstance()->stop();
    ThreadPool::getInstance()->start(numThread);

    return numThread;
}

extern "C" __declspec(dllexport)
 LPXLOPER12 xPutDupire(
    //  model parameters
    double              spot,
    FP12*               spots,
    FP12*               times,
    FP12*               vols,
    double              maxDt,
    LPXLOPER12          xid)
{
    FreeAllTempMemory();

    const string id = getString(xid);

    //  Make sure we have an id
    if (maxDt <= 0.0 || id.empty()) return TempErr12(xlerrNA);

    //  Unpack

    if (spots->rows * spots->columns * times->rows * times->columns != vols->rows * vols->columns)
    {
        return TempErr12(xlerrNA);
    }

    vector<double> vspots = to_vector(spots);
    vector<double> vtimes = to_vector(times);
    matrix<double> vvols = to_matrix(vols);

    //  Call and return
    putDupire(spot, vspots, vtimes, vvols, maxDt, id);

    return TempStr12(id);
}

extern "C" __declspec(dllexport)
LPXLOPER12 xPutBlackScholes(
    double              spot,
    double              vol,
    double              qSpot,
    double              rate,
    double              div,
    LPXLOPER12          xid)
{
    FreeAllTempMemory();

    const string id = getString(xid);

    //  Make sure we have an id
    if (id.empty()) return TempErr12(xlerrNA);

    //  Call and return
    putBlackScholes(spot, vol, qSpot > 0, rate, div, id);

    return TempStr12(id);
}

extern "C" __declspec(dllexport)
 LPXLOPER12 xPutEuropean(
    double              strike,
    double              exerciseDate,
    double              settlementDate,
     LPXLOPER12         xid)
{
    FreeAllTempMemory();

    const string id = getString(xid);

    //  Make sure we have an id
    if (id.empty()) return TempErr12(xlerrNA);

    if (settlementDate <= 0) settlementDate = exerciseDate;

    //  Call and return
    putEuropean(strike, exerciseDate, settlementDate, id);

    return TempStr12(id);
}

extern "C" __declspec(dllexport)
LPXLOPER12 xPutBarrier(
    double              strike,
    double              barrier,
    double              maturity,
    double              monitorFreq,
    double              smoothing,
    LPXLOPER12          xid)
{
    FreeAllTempMemory();

    const string id = getString(xid);

    //  Make sure we have an id
    if (id.empty()) return TempErr12(xlerrNA);

    //  Call and return
    putBarrier(strike, barrier, maturity, monitorFreq, smoothing, id);

    return TempStr12(id);
}

extern "C" __declspec(dllexport)
LPXLOPER12 xPutContingent(
    double              coupon,
    double              maturity,
    double              payFreq,
    double              smoothing,
    LPXLOPER12          xid)
{
    FreeAllTempMemory();

    const string id = getString(xid);

    //  Make sure we have an id
    if (id.empty()) return TempErr12(xlerrNA);

    //  Call and return
    putContingent(coupon, maturity, payFreq, smoothing, id);

    return TempStr12(id);
}

extern "C" __declspec(dllexport)
LPXLOPER12 xPutEuropeans(
    FP12*               maturities,
    FP12*               strikes,
    LPXLOPER12          xid)
{
    FreeAllTempMemory();

    const string id = getString(xid);

    //  Make sure we have an id
    if (id.empty()) return TempErr12(xlerrNA);

    //  Make sure we have options
    if (!strikes->rows || !strikes->columns) return TempErr12(xlerrNA);

    //  Maturities and strikes, removing blanks
    vector<double> vmats;
    vector<double> vstrikes;
    {
        size_t rows = maturities->rows;
        size_t cols = maturities->columns;
        
        if (strikes->rows * strikes->columns != rows * cols) return TempErr12(xlerrNA);
        
        double* mats= maturities->array;
        double* strs = strikes->array;

        for (size_t i = 0; i < cols * rows; ++i)
        {
            if (mats[i] > EPS && strs[i] > EPS)
            {
                vmats.push_back(mats[i]);
                vstrikes.push_back(strs[i]);
            }
        }
    }

    //  Make sure we still have options
    if (vmats.empty()) return TempErr12(xlerrNA);

    //  Call and return
    putEuropeans(vmats, vstrikes, id);

    return TempStr12(id);
}

//  Access payoff identifiers and parameters

extern "C" __declspec(dllexport)
LPXLOPER12 xPayoffIds(
    LPXLOPER12          xid)
{
    FreeAllTempMemory();

    const string id = getString(xid);
    //  Make sure we have an id
    if (id.empty()) return TempErr12(xlerrNA);

    const auto* prd = getProduct<double>(id);
    //  Make sure we have a product
    if (!prd) return TempErr12(xlerrNA);

    return from_strVector(prd->payoffLabels());
}

extern "C" __declspec(dllexport)
LPXLOPER12 xParameters(
    LPXLOPER12          xid)
{
    FreeAllTempMemory();

    const string id = getString(xid);
    //  Make sure we have an id
    if (id.empty()) return TempErr12(xlerrNA);

    Model<double>* mdl = const_cast<Model<double>*>(getModel<double>(id));
    //  Make sure we have a model
    if (!mdl) return TempErr12(xlerrNA);

    const auto& paramLabels = mdl->parameterLabels();
    const auto& params = mdl->parameters();
    vector<double> paramsCopy(params.size());
    transform(params.begin(), params.end(), paramsCopy.begin(), [](const double* p) { return *p; });

    return from_labelsAndNumbers(paramLabels, paramsCopy);
}

extern "C" __declspec(dllexport)
LPXLOPER12 xValue(
    LPXLOPER12          modelid,
    LPXLOPER12          productid,
    //  numerical parameters
    double              useSobol,
    double              seed1,
    double              seed2,
    double              numPath,
    double              parallel)
{
    FreeAllTempMemory();

    const string pid = getString(productid);
    //  Make sure we have an id
    if (pid.empty()) return TempErr12(xlerrNA);

    const auto* prd = getProduct<double>(pid);
    //  Make sure we have a product
    if (!prd) return TempErr12(xlerrNA);

    const string mid = getString(modelid);
    //  Make sure we have an id
    if (mid.empty()) return TempErr12(xlerrNA);

    Model<double>* mdl = const_cast<Model<double>*>(getModel<double>(mid));
    //  Make sure we have a model
    if (!mdl) return TempErr12(xlerrNA);

    //  Numerical params
    const auto num = xl2num(useSobol, seed1, seed2, numPath, parallel);
    //  Make sure we have a numPath
    if (!num.numPath) return TempErr12(xlerrNA);

    //  Call and return;
    try 
    {
        auto results = value(mid, pid, num);
        return from_labelsAndNumbers(results.identifiers, results.values);
    }
    catch (const exception&)
    {
        return TempErr12(xlerrNA);
    }
}

extern "C" __declspec(dllexport)
LPXLOPER12 xValueTime(
	LPXLOPER12          modelid,
	LPXLOPER12          productid,
	//  numerical parameters
	double              useSobol,
	double              seed1,
	double              seed2,
	double              numPath,
	double              parallel)
{
	FreeAllTempMemory();

	const string pid = getString(productid);
	//  Make sure we have an id
	if (pid.empty()) return TempErr12(xlerrNA);

	const auto* prd = getProduct<double>(pid);
	//  Make sure we have a product
	if (!prd) return TempErr12(xlerrNA);

	const string mid = getString(modelid);
	//  Make sure we have an id
	if (mid.empty()) return TempErr12(xlerrNA);

	Model<double>* mdl = const_cast<Model<double>*>(getModel<double>(mid));
	//  Make sure we have a model
	if (!mdl) return TempErr12(xlerrNA);

	//  Numerical params
	const auto num = xl2num(useSobol, seed1, seed2, numPath, parallel);
	//  Make sure we have a numPath
	if (!num.numPath) return TempErr12(xlerrNA);

	//  Call and return;
	try
	{
		clock_t t0 = clock();
		auto results = value(mid, pid, num);
		clock_t t1 = clock();
		LPXLOPER12 oper = TempXLOPER12();
		resize(oper, results.values.size() + 1, 1);
		for (int i = 0; i < results.values.size(); ++i) setNum(oper, results.values[i], i, 0);
		setNum(oper, t1 - t0, results.values.size(), 0);

		return oper;
	}
	catch (const exception&)
	{
		return TempErr12(xlerrNA);
	}
}

extern "C" __declspec(dllexport)
LPXLOPER12 xAADrisk(
    LPXLOPER12          modelid,
    LPXLOPER12          productid,
    LPXLOPER12          xRiskPayoff,
    //  numerical parameters
    double              useSobol,
    double              seed1,
    double              seed2,
    double              numPath,
    double              parallel)
{
    FreeAllTempMemory();

    const string pid = getString(productid);
    //  Make sure we have an id
    if (pid.empty()) return TempErr12(xlerrNA);

    const string mid = getString(modelid);
    //  Make sure we have an id
    if (mid.empty()) return TempErr12(xlerrNA);

    //  Numerical params
    const auto num = xl2num(useSobol, seed1, seed2, numPath, parallel);
    //  Make sure we have a numPath
    if (!num.numPath) return TempErr12(xlerrNA);

    //  Risk payoff
    const string riskPayoff = getString(xRiskPayoff);

    try
    {
        auto results = AADriskOne(mid, pid, num, riskPayoff);
		const size_t n = results.risks.size(), N = n + 1;

        LPXLOPER12 oper = TempXLOPER12();
        resize(oper, N, 2);

        setString(oper, "value", 0, 0);
        setNum(oper, results.riskPayoffValue, 0, 1);

        for (size_t i = 0; i < n; ++i)
        {
            setString(oper, results.paramIds[i], i + 1, 0);
            setNum(oper, results.risks[i], i + 1, 1);
        }

        return oper;
    }
    catch (const exception&)
    {
        return TempErr12(xlerrNA);
    }
}

extern "C" __declspec(dllexport)
LPXLOPER12 xAADriskAggregate(
    LPXLOPER12          modelid,
    LPXLOPER12          productid,
    LPXLOPER12          xPayoffs,
    FP12*               xNotionals,
    //  numerical parameters
    double              useSobol,
    double              seed1,
    double              seed2,
    double              numPath,
    double              parallel)
{
    FreeAllTempMemory();

    const string pid = getString(productid);
    //  Make sure we have an id
    if (pid.empty()) return TempErr12(xlerrNA);

    const string mid = getString(modelid);
    //  Make sure we have an id
    if (mid.empty()) return TempErr12(xlerrNA);

    //  Numerical params
    const auto num = xl2num(useSobol, seed1, seed2, numPath, parallel);
    //  Make sure we have a numPath
    if (!num.numPath) return TempErr12(xlerrNA);

    //  Payoffs and notionals, removing blanks
    map<string, double> notionals;
    {
        size_t rows = getRows(xPayoffs);
        size_t cols = getCols(xPayoffs);

        if (rows * cols == 0 ||
            xNotionals->rows * xNotionals->columns != rows * cols) 
                return TempErr12(xlerrNA);

        size_t idx = 0;
        for (size_t i = 0; i < rows; ++i) for (size_t j = 0; j < cols; ++j)
        {
            string payoff = getString(xPayoffs, i, j);
            double notional = xNotionals->array[idx++];
            if (payoff != "" && fabs(notional) > EPS)
            {
                notionals[payoff] = notional;
            }
        }
    }

    try
    {
        auto results = AADriskAggregate(mid, pid, notionals, num);
        const size_t n = results.risks.size(), N = n + 1;

        LPXLOPER12 oper = TempXLOPER12();
        resize(oper, N, 2);

        setString(oper, "value", 0, 0);
        setNum(oper, results.riskPayoffValue, 0, 1);

        for (size_t i = 0; i < n; ++i)
        {
            setString(oper, results.paramIds[i], i + 1, 0);
            setNum(oper, results.risks[i], i + 1, 1);
        }

        return oper;
    }
    catch (const exception&)
    {
        return TempErr12(xlerrNA);
    }
}

unordered_map<string, RiskReports> riskStore;

extern "C" __declspec(dllexport)
LPXLOPER12 xBumprisk(
    LPXLOPER12          modelid,
    LPXLOPER12          productid,
    //  numerical parameters
    double              useSobol,
    double              seed1,
    double              seed2,
    double              numPath,
    double              parallel,
    //  display now or put in memory?
    double              displayNow,
    LPXLOPER12          storeid)
{
    FreeAllTempMemory();

    const string pid = getString(productid);
    //  Make sure we have an id
    if (pid.empty()) return TempErr12(xlerrNA);

    const string mid = getString(modelid);
    //  Make sure we have an id
    if (mid.empty()) return TempErr12(xlerrNA);

    //  Numerical params
    const auto num = xl2num(useSobol, seed1, seed2, numPath, parallel);
    //  Make sure we have a numPath
    if (!num.numPath) return TempErr12(xlerrNA);

    try
    {
        auto results = bumpRisk(mid, pid, num);
        if (displayNow > 0.5)
        {
            return from_labelledMatrix(results.params, results.payoffs, results.risks, "value", results.values);
        }
        else
        {
            const string riskId = getString(storeid);
            if (riskId == "") return TempErr12(xlerrNA);
            riskStore[riskId] = results;
            return storeid;
        }
    }
    catch (const exception&)
    {
        return TempErr12(xlerrNA);
    }
}

extern "C" __declspec(dllexport)
LPXLOPER12 xAADriskMulti(
    LPXLOPER12          modelid,
    LPXLOPER12          productid,
    //  numerical parameters
    double              useSobol,
    double              seed1,
    double              seed2,
    double              numPath,
    double              parallel,
    //  display now or put in memory?
    double              displayNow,
    LPXLOPER12          storeid)
{
    FreeAllTempMemory();

    const string pid = getString(productid);
    //  Make sure we have an id
    if (pid.empty()) return TempErr12(xlerrNA);

    const string mid = getString(modelid);
    //  Make sure we have an id
    if (mid.empty()) return TempErr12(xlerrNA);

    //  Numerical params
    const auto num = xl2num(useSobol, seed1, seed2, numPath, parallel);
    //  Make sure we have a numPath
    if (!num.numPath) return TempErr12(xlerrNA);

    try
    {
        auto results = AADriskMulti(mid, pid, num);
        if (displayNow > 0.5)
        {
            return from_labelledMatrix(results.params, results.payoffs, results.risks, "value", results.values);
        }
        else
        {
            const string riskId = getString(storeid);
            if (riskId == "") return TempErr12(xlerrNA);
            riskStore[riskId] = results;
            return storeid;
        }
    }
    catch (const exception&)
    {
        return TempErr12(xlerrNA);
    }
}

extern "C" __declspec(dllexport)
LPXLOPER12 xDisplayRisk(
    LPXLOPER12          riskid,
    LPXLOPER12          displayid)
{
    FreeAllTempMemory();

    RiskReports* results;
    const string riskId = getString(riskid);
    auto& it = riskStore.find(riskId);
    if (it == riskStore.end()) return TempErr12(xlerrNA);
    else results = &it->second;

    const vector<string> riskIds = to_strVector(displayid);
    
    vector<size_t> riskCols;
    
    for (const auto& id : riskIds)
    {
        auto it2 = find(results->payoffs.begin(), results->payoffs.end(), id);
        if (it2 == results->payoffs.end()) return TempErr12(xlerrNA);
        riskCols.push_back(distance(results->payoffs.begin(), it2));
    }
    if (riskCols.empty()) return TempErr12(xlerrNA);

    const size_t nParam = results->risks.rows();
    const size_t nPayoffs = riskCols.size();

    vector<double> showVals(nPayoffs);
    matrix<double> showRisks(nParam, nPayoffs);
    for (size_t j = 0; j < nPayoffs; ++j)
    {
        for (size_t i = 0; i < nParam; ++i)
        {
            showRisks[i][j] = results->risks[i][riskCols[j]];
        }
        showVals[j] = results->values[riskCols[j]];
    }

    return from_labelledMatrix(results->params, riskIds, showRisks, "value", showVals);
}

extern "C" __declspec(dllexport)
 LPXLOPER12 xDupireCalib(
    //  Merton market parameters
    const double spot,
    const double vol,
    const double jmpIntens,
    const double jmpAverage,
    const double jmpStd,
    //  Discretization
    FP12* spots,
    const double maxDs,
    FP12* times,
    const double maxDt)
{
    FreeAllTempMemory();

    //  Make sure the last input is given
    if (maxDs == 0 || maxDt == 0) return 0;

    //  Unpack
    vector<double> vspots;
    {
        size_t rows = spots->rows;
        size_t cols = spots->columns;
        double* numbers = spots->array;

        vspots.resize(rows*cols);
        copy(numbers, numbers + rows*cols, vspots.begin());
    }

    vector<double> vtimes;
    {
        size_t rows = times->rows;
        size_t cols = times->columns;
        double* numbers = times->array;

        vtimes.resize(rows*cols);
        copy(numbers, numbers + rows*cols, vtimes.begin());
    }

    auto results = dupireCalib(vspots, maxDs, vtimes, maxDt,
        spot, vol, jmpIntens, jmpAverage, jmpStd);

    //  Return
    return from_labelledMatrix(results.spots, results.times, results.lVols);
}

extern "C" __declspec(dllexport)
 LPXLOPER12 xDupireSuperbucket(
    //  Merton market parameters
    const double spot,
    const double vol,
    const double jmpIntens,
    const double jmpAverage,
    const double jmpStd,
    //  risk view
    FP12* strikes,
    FP12* mats,
    //  calibration
    FP12* spots,
    double maxDs,
    FP12* times,
    double maxDtVol,
    //  MC
    double              maxDtSimul,
    //  product 
    LPXLOPER12          productid,
    LPXLOPER12          xPayoffs,
    FP12*               xNotionals,
    //  numerical parameters
    double              useSobol,
    double              seed1,
    double              seed2,
    double              numPath,
    double              parallel,
    //  bump or AAD?
    double              bump)

{
    FreeAllTempMemory();

    //  Make sure the last input is given
    if (numPath <= 0) return TempErr12(xlerrNA);

    //  Make sure we have a product
    const string pid = getString(productid);
    //  Make sure we have an id
    if (pid.empty()) return TempErr12(xlerrNA);

    //  Unpack
    vector<double> vstrikes = to_vector(strikes);
    vector<double> vmats = to_vector(mats);
    vector<double> vspots = to_vector(spots);
    vector<double> vtimes = to_vector(times);

    //  Numerical params
    const auto num = xl2num(useSobol, seed1, seed2, numPath, parallel);
    //  Make sure we have a numPath
    if (!num.numPath) return TempErr12(xlerrNA);

    //  Payoffs and notionals, removing blanks
    map<string, double> notionals;
    {
        size_t rows = getRows(xPayoffs);
        size_t cols = getCols(xPayoffs);;

        if (rows * cols == 0 ||
            xNotionals->rows * xNotionals->columns != rows * cols)
            return TempErr12(xlerrNA);

        size_t idx = 0;
        for (size_t i = 0; i < rows; ++i) for (size_t j = 0; j < cols; ++j)
        {
            string payoff = getString(xPayoffs, i, j);
            double notional = xNotionals->array[idx++];
            if (payoff != "" && fabs(notional) > EPS)
            {
                notionals[payoff] = notional;
            }
        }
    }

    try
    {
        auto results = bump < EPS
            ? dupireSuperbucket(
                spot,
                maxDtSimul,
                pid,
                notionals,
                vspots,
                maxDs,
                vtimes,
                maxDtVol,
                vstrikes,
                vmats,
                vol,
                jmpIntens,
                jmpAverage,
                jmpStd,
                num)
            : dupireSuperbucketBump(
                spot,
                maxDtSimul,
                pid,
                notionals,
                vspots,
                maxDs,
                vtimes,
                maxDtVol,
                vstrikes,
                vmats,
                vol,
                jmpIntens,
                jmpAverage,
                jmpStd,
                num);

        //  Build return
        const size_t n = results.vega.rows(), m = results.vega.cols();
        const size_t N = n + 4, M = m + 2;

        LPXLOPER12 oper = TempXLOPER12();
        resize(oper, N, M);

        setString(oper, "value", 0, 0);
        setNum(oper, results.value, 0, 1);
        setString(oper, "delta", 1, 0);
        setNum(oper, results.delta, 1, 1);
        setString(oper, "vega", 2, 0);
        setString(oper, "mats", 2, 1);
        for (size_t i = 0; i < m; ++i) setNum(oper, results.mats[i], 2, 2 + i);
        setString(oper, "strikes", 3, 0);
        for (size_t i = 0; i < n; ++i)
        {
            setNum(oper, results.strikes[i], 4 + i, 0);
            for (size_t j = 0; j < m; ++j) setNum(oper, results.vega[i][j], 4 + i, 2 + j);
        }

        // Return it
        return oper;
    }
    catch (const exception&)
    {
        return TempErr12(xlerrNA);
    }
}

extern "C" __declspec(dllexport)
double xMerton(double spot, double vol, double mat, double strike, double intens, double meanJmp, double stdJmp)
{
    return merton(spot, strike, vol, mat, intens, meanJmp, stdJmp);
}

extern "C" __declspec(dllexport)
double xBarrierBlackScholes(double spot, double rate, double div, double vol, double mat, double strike, double barrier)
{
    return BlackScholesKO(spot, rate, div, strike, barrier, mat, vol);
}

extern "C" __declspec(dllexport)
double xToyDupireBarrierMc(
    //  model parameters
    double              spot,
    FP12*               spots,
    FP12*               times,
    FP12*               vols,
    double              mat,
    double              strike,
    double              barrier,
    double              paths,
    double              steps,
    double              epsilon,
    double              useSobol,
    double              seed1,
    double              seed2)
{
    FreeAllTempMemory();

    //  Make sure we have paths and steps
    if (paths <= 0.0 || steps <= 0.0) return -1;

    //  Unpack

    if (spots->rows * spots->columns * times->rows * times->columns != vols->rows * vols->columns)
    {
        return -1;
    }

    vector<double> vspots = to_vector(spots);
    vector<double> vtimes = to_vector(times);
    matrix<double> vvols = to_matrix(vols);

    //  Random Number Generator
    unique_ptr<RNG> rng;
    if (useSobol > 0.5) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1 > 0.5 ? int(seed1): 12345, seed2 > 0.5? int(seed2): 123456);
	rng->init(int(steps));

    //  Call and return
    return toyDupireBarrierMc(spot, vspots, vtimes, vvols, mat, strike, barrier, int(paths), int(steps), 100*epsilon, *rng);
}

extern "C" __declspec(dllexport)
LPXLOPER12 xToyDupireBarrierMcRisks(
    //  model parameters
    double              spot,
    FP12*               spots,
    FP12*               times,
    FP12*               vols,
    double              mat,
    double              strike,
    double              barrier,
    double              paths,
    double              steps,
    double              epsilon,
    double              useSobol,
    double              seed1,
    double              seed2)
{
    FreeAllTempMemory();

    //  Make sure we have paths and steps
    if (paths <= 0.0 || steps <= 0.0) TempErr12(xlerrNA);

    //  Unpack

    if (spots->rows * spots->columns * times->rows * times->columns != vols->rows * vols->columns)
    {
        return TempErr12(xlerrNA);
    }

    vector<double> vspots = to_vector(spots);
    vector<double> vtimes = to_vector(times);
    matrix<double> vvols = to_matrix(vols);

    //  Random Number Generator
    unique_ptr<RNG> rng;
    if (useSobol > 0.5) rng = make_unique<Sobol>();
    else rng = make_unique<mrg32k3a>(seed1 > 0.5 ? int(seed1): 12345, seed2 > 0.5? int(seed2): 123456);
	rng->init(int(steps));

    //  Call 
	double price, delta;
	matrix<double> vegas(vvols.rows(), vvols.cols());
    toyDupireBarrierMcRisks(spot, vspots, vtimes, vvols, mat, strike, barrier, int(paths), int(steps), 100*epsilon, *rng,
		price, delta, vegas);

	//	Pack and return
	LPXLOPER12 results = TempXLOPER12();
	resize(results, 2 + vegas.rows()*vegas.cols(), 1);
	setNum(results, price, 0, 0);
	setNum(results, delta, 1, 0);
	size_t r = 2;
	for (int i=0; i<vegas.rows(); ++i) for (int j=0; j<vegas.cols(); ++j)
	{
		setNum(results, vegas[i][j], r, 0);
		++r;
	}

	return results;
}


//	Registers

extern "C" __declspec(dllexport) int xlAutoOpen(void)
{
	XLOPER12 xDLL;

	Excel12f(xlGetName, &xDLL, 0);

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xRestartThreadPool"),
        (LPXLOPER12)TempStr12(L"BB"),
        (LPXLOPER12)TempStr12(L"xRestartThreadPool"),
        (LPXLOPER12)TempStr12(L"numThreads"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Restarts the thread pool with n threads"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xPutBlackScholes"),
        (LPXLOPER12)TempStr12(L"QBBBBBQ"),
        (LPXLOPER12)TempStr12(L"xPutBlackScholes"),
        (LPXLOPER12)TempStr12(L"spot, vol, qSpot, rate, div, id"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Initializes a Black-Scholes in memory"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xPutDupire"),
        (LPXLOPER12)TempStr12(L"QBK%K%K%BQ"),
        (LPXLOPER12)TempStr12(L"xPutDupire"),
        (LPXLOPER12)TempStr12(L"spot, spots, times, vols, maxDt, id"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Initializes a Dupire in memory"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xPutEuropean"),
        (LPXLOPER12)TempStr12(L"QBBBQ"),
        (LPXLOPER12)TempStr12(L"xPutEuropean"),
        (LPXLOPER12)TempStr12(L"strike, exerciseDate, [settlementDate], id"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Initializes a European call in memory"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xPutBarrier"),
        (LPXLOPER12)TempStr12(L"QBBBBBQ"),
        (LPXLOPER12)TempStr12(L"xPutBarrier"),
        (LPXLOPER12)TempStr12(L"strike, barrier, maturity, monitoringFreq, [smoothingFactor], id"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Initializes a European call in memory"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xPutContingent"),
        (LPXLOPER12)TempStr12(L"QBBBBQ"),
        (LPXLOPER12)TempStr12(L"xPutContingent"),
        (LPXLOPER12)TempStr12(L"coupon, maturity, payFreq, [smoothingFactor], id"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Initializes a European call in memory"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xPutEuropeans"),
        (LPXLOPER12)TempStr12(L"QK%K%Q"),
        (LPXLOPER12)TempStr12(L"xPutEuropeans"),
        (LPXLOPER12)TempStr12(L"maturities, strikes, id"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Initializes a collection of European call in memory"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xPayoffIds"),
        (LPXLOPER12)TempStr12(L"QQ"),
        (LPXLOPER12)TempStr12(L"xPayoffIds"),
        (LPXLOPER12)TempStr12(L"id"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"The payoff identifiers in a product"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xParameters"),
        (LPXLOPER12)TempStr12(L"QQ"),
        (LPXLOPER12)TempStr12(L"xParameters"),
        (LPXLOPER12)TempStr12(L"id"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"The parameters of a model"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xValue"),
        (LPXLOPER12)TempStr12(L"QQQBBBBB"),
        (LPXLOPER12)TempStr12(L"xValue"),
        (LPXLOPER12)TempStr12(L"modelId, productId, useSobol, [seed1], [seed2], N, [Parallel]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Monte-Carlo valuation"),
        (LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xValueTime"),
		(LPXLOPER12)TempStr12(L"QQQBBBBB"),
		(LPXLOPER12)TempStr12(L"xValueTime"),
		(LPXLOPER12)TempStr12(L"modelId, productId, useSobol, [seed1], [seed2], N, [Parallel]"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Timed Monte-Carlo valuation"),
		(LPXLOPER12)TempStr12(L""));
	
	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xAADrisk"),
        (LPXLOPER12)TempStr12(L"QQQQBBBBB"),
        (LPXLOPER12)TempStr12(L"xAADrisk"),
        (LPXLOPER12)TempStr12(L"modelId, productId, riskPayoff, useSobol, [seed1], [seed2], N, [Parallel]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"AAD risk report"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xAADriskAggregate"),
        (LPXLOPER12)TempStr12(L"QQQQK%BBBBB"),
        (LPXLOPER12)TempStr12(L"xAADriskAggregate"),
        (LPXLOPER12)TempStr12(L"modelId, productId, payoffs, notionals, useSobol, [seed1], [seed2], N, [Parallel]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"AAD risk report for aggregate book of payoffs"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xBumprisk"),
        (LPXLOPER12)TempStr12(L"QQQBBBBBBQ"),
        (LPXLOPER12)TempStr12(L"xBumprisk"),
        (LPXLOPER12)TempStr12(L"modelId, productId, useSobol, [seed1], [seed2], N, [Parallel], [display?], [storeId]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Bump risk report"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xAADriskMulti"),
        (LPXLOPER12)TempStr12(L"QQQBBBBBBQ"),
        (LPXLOPER12)TempStr12(L"xAADriskMulti"),
        (LPXLOPER12)TempStr12(L"modelId, productId, useSobol, [seed1], [seed2], N, [Parallel], [display?], [storeId]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"AAD risk report for multiple payoffs"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xDisplayRisk"),
        (LPXLOPER12)TempStr12(L"QQQ"),
        (LPXLOPER12)TempStr12(L"xDisplayRisk"),
        (LPXLOPER12)TempStr12(L"reportId, payoffId"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Display risk from stored report"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xDupireSuperbucket"),
        (LPXLOPER12)TempStr12(L"QBBBBBK%K%K%BK%BBQQK%BBBBBB"),
        (LPXLOPER12)TempStr12(L"xDupireSuperbucket"),
        (LPXLOPER12)TempStr12(L"spot, vol, jmpIt, jmpAve, jmpStd, RiskStrikes, riskMats, volSpots, maxDs, volTimes, maxDtVol, maxDtSimul, productId, payoffs, notionals, sobol, s1, s2, numPth, parallel, [bump?]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Computes the risk of a product in Dupire"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xDupireCalib"),
        (LPXLOPER12)TempStr12(L"QBBBBBK%BK%B"),
        (LPXLOPER12)TempStr12(L"xDupireCalib"),
        (LPXLOPER12)TempStr12(L"spot, vol, jumpIntensity, jumpAverage, jumpStd, spots, maxds, times, mxdt"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Calibrates Dupire"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xMerton"),
        (LPXLOPER12)TempStr12(L"BBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xMerton"),
        (LPXLOPER12)TempStr12(L"spot, vol, mat, strike, intens, meanJmp, stdJmp"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Merton"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xBarrierBlackScholes"),
        (LPXLOPER12)TempStr12(L"BBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xBarrierBlackScholes"),
        (LPXLOPER12)TempStr12(L"spot, rate, div, vol, mat, strike, barrier"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Barrier analytic"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xToyDupireBarrierMc"),
        (LPXLOPER12)TempStr12(L"BBK%K%K%BBBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xToyDupireBarrierMc"),
        (LPXLOPER12)TempStr12(L"spot, spots, times, vols, mat, strike, barrier, paths, steps, epsilon, useSobol, [seed1], [seed2]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Toy Dupire Barrier MC"),
        (LPXLOPER12)TempStr12(L""));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xToyDupireBarrierMcRisks"),
        (LPXLOPER12)TempStr12(L"QBK%K%K%BBBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xToyDupireBarrierMcRisks"),
        (LPXLOPER12)TempStr12(L"spot, spots, times, vols, mat, strike, barrier, paths, steps, epsilon, useSobol, [seed1], [seed2]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Toy Dupire Barrier MC AAD risks"),
        (LPXLOPER12)TempStr12(L""));

	/* Free the XLL filename */
	Excel12f(xlFree, 0, 1, (LPXLOPER12)&xDLL);

    /*  Start the thread pool   */
    ThreadPool::getInstance()->start(thread::hardware_concurrency() - 1);

	return 1;
}

extern "C" __declspec(dllexport) int xlAutoClose(void)
{
    /*  Stop the thread pool   */
    ThreadPool::getInstance()->stop();

    return 1;
}