
#include <windows.h>
#include "xlcall.h"
#include "framework.h"

#include "threadPool.h"

#include "mcDriver.h"

//	Wrappers

extern "C" __declspec(dllexport)
inline double xUocDupire(
    //  model parameters
    double              spot,
    FP12*               spots,
    FP12*               times,
    FP12*               vols,   
    double              maxDt,
    //  product parameters
    double              strike,
    double              barrier,
    double              maturity,
    double              monitorFreq,
    //  numerical parameters
    double              useSobol,
    double              useAnti,
    double              seed1,
    double              seed2,
    double              numPath,
    double              parallel)
{
    //  Make sure the last input is given
    if (!numPath) return 0;

    //  MaxDt = 0 crashes
    if (maxDt <= 0) return 0;

    //  Same thing for frequency
    if (monitorFreq <= 0) return 0;
    
    //  Default seeds
    if (seed1 <= 0) seed1 = 12345;
    if (seed2 <= 0) seed2 = 12346;

    //  Unpack

    vector<double> vspots;
    {
        size_t rows = spots->rows;
        size_t cols = spots->columns;
        double* numbers = spots->array;

        vspots.resize(rows);
        copy(numbers, numbers + rows, vspots.begin());
    }

    vector<double> vtimes;
    {
        size_t rows = times->rows;
        size_t cols = times->columns;
        double* numbers = times->array;

        vtimes.resize(cols);
        copy(numbers, numbers + cols, vtimes.begin());
    }

    matrix<double> vvols;
    {
        size_t rows = vols->rows;
        size_t cols = vols->columns;
        double* numbers = vols->array;

        vvols.resize(rows, cols);
        copy(numbers, numbers + rows * cols, vvols.begin());
    }

    //  Call and return

    return uocDupire(spot, vspots, vtimes, vvols, maxDt, strike, barrier, maturity, monitorFreq, parallel>0,
        useSobol > 0, static_cast<int>(numPath), useAnti > 0, static_cast<int>(seed1), static_cast<int>(seed2));

}

extern "C" __declspec(dllexport)
inline FP12* xUocDupireBump(
    //  model parameters
    double              spot,
    FP12*               spots,
    FP12*               times,
    FP12*               vols,
    double              maxDt,
    //  product parameters
    double              strike,
    double              barrier,
    double              maturity,
    double              monitorFreq,
    //  numerical parameters
    double              useSobol,
    double              useAnti,
    double              seed1,
    double              seed2,
    double              numPath,
    double              parallel)
{
    //  Make sure the last input is given
    if (!numPath) return 0;

    //  MaxDt = 0 crashes
    if (maxDt <= 0) return 0;

    //  Same thing for frequency
    if (monitorFreq <= 0) return 0;

    //  Default seeds
    if (seed1 <= 0) seed1 = 12345;
    if (seed2 <= 0) seed2 = 12346;

    //  Unpack

    vector<double> vspots;
    {
        size_t rows = spots->rows;
        size_t cols = spots->columns;
        double* numbers = spots->array;

        vspots.resize(rows);
        copy(numbers, numbers + rows, vspots.begin());
    }

    vector<double> vtimes;
    {
        size_t rows = times->rows;
        size_t cols = times->columns;
        double* numbers = times->array;

        vtimes.resize(cols);
        copy(numbers, numbers + cols, vtimes.begin());
    }

    matrix<double> vvols;
    {
        size_t rows = vols->rows;
        size_t cols = vols->columns;
        double* numbers = vols->array;

        vvols.resize(rows, cols);
        copy(numbers, numbers + rows * cols, vvols.begin());
    }

    //  Call 
    auto res = uocDupireBumpRisk(spot, vspots, vtimes, vvols, maxDt, 
        strike, barrier, maturity, monitorFreq, parallel>0,
        useSobol > 0, static_cast<int>(numPath), useAnti > 0,
        static_cast<int>(seed1), static_cast<int>(seed2));
    //  Build return

    // Allocate result
    // Calculate size
    size_t resultRows = res.vega.rows()+2, resultCols = res.vega.cols(), resultSize = resultRows * resultCols;
    
    // Return an error if size is 0
    if (resultSize <= 0) return nullptr;

    // First, free all memory previously allocated
    // the function is defined in framework.h
    FreeAllTempMemory();
    // Then, allocate the memory for this result
    // We don't need to de-allocate it, that will be done by the next call
    // to this function or another function calling FreeAllTempMemory()
    // Memory size required, details in Dalton's book, section 6.2.2
    size_t memSize = sizeof(FP12) + (resultSize - 1) * sizeof(double);
    // Finally allocate, function definition in framework.h
    FP12* result = (FP12*)GetTempMemory(memSize);
    // Compute result
    result->rows = resultRows;
    result->columns = resultCols;
    for (size_t i = 0; i < resultSize; ++i) result->array[i] = 0.0;
    result->array[0] = res.value;
    result->array[resultCols] = res.delta;
    for (size_t i = 0; i < res.vega.rows(); ++i) 
        for (size_t j = 0; j < res.vega.cols(); ++j)
    {
        result->array[(i+2)*resultCols + j] = res.vega[i][j];
    }

    // Return it
    return result;
}

extern "C" __declspec(dllexport)
    inline FP12* xUocDupireAAD(
        //  model parameters
        double              spot,
        FP12*               spots,
        FP12*               times,
        FP12*               vols,
        double              maxDt,
        //  product parameters
        double              strike,
        double              barrier,
        double              maturity,
        double              monitorFreq,
        //  numerical parameters
        double              useSobol,
        double              useAnti,
        double              seed1,
        double              seed2,
        double              numPath,
        double              parallel)
{
    //  Make sure the last input is given
    if (!numPath) return 0;

    //  MaxDt = 0 crashes
    if (maxDt <= 0) return 0;

    //  Same thing for frequency
    if (monitorFreq <= 0) return 0;

    //  Default seeds
    if (seed1 <= 0) seed1 = 12345;
    if (seed2 <= 0) seed2 = 12346;

    //  Unpack

    vector<double> vspots;
    {
        size_t rows = spots->rows;
        size_t cols = spots->columns;
        double* numbers = spots->array;

        vspots.resize(rows);
        copy(numbers, numbers + rows, vspots.begin());
    }

    vector<double> vtimes;
    {
        size_t rows = times->rows;
        size_t cols = times->columns;
        double* numbers = times->array;

        vtimes.resize(cols);
        copy(numbers, numbers + cols, vtimes.begin());
    }

    matrix<double> vvols;
    {
        size_t rows = vols->rows;
        size_t cols = vols->columns;
        double* numbers = vols->array;

        vvols.resize(rows, cols);
        copy(numbers, numbers + rows * cols, vvols.begin());
    }

    //  Call 
    auto res = uocDupireAADRisk(spot, vspots, vtimes, vvols, maxDt,
        strike, barrier, maturity, monitorFreq, parallel>0,
        useSobol > 0, static_cast<int>(numPath), useAnti > 0,
        static_cast<int>(seed1), static_cast<int>(seed2));

    //  Build return

    // Allocate result
    // Calculate size
    size_t resultRows = res.vega.rows() + 2, resultCols = res.vega.cols(), resultSize = resultRows * resultCols;

    // Return an error if size is 0
    if (resultSize <= 0) return nullptr;

    // First, free all memory previously allocated
    // the function is defined in framework.h
    FreeAllTempMemory();
    // Then, allocate the memory for this result
    // We don't need to de-allocate it, that will be done by the next call
    // to this function or another function calling FreeAllTempMemory()
    // Memory size required, details in Dalton's book, section 6.2.2
    size_t memSize = sizeof(FP12) + (resultSize - 1) * sizeof(double);
    // Finally allocate, function definition in framework.h
    FP12* result = (FP12*)GetTempMemory(memSize);
    // Compute result
    result->rows = resultRows;
    result->columns = resultCols;
    for (size_t i = 0; i < resultSize; ++i) result->array[i] = 0.0;
    result->array[0] = res.value;
    result->array[resultCols] = res.delta;
    for (size_t i = 0; i < res.vega.rows(); ++i) 
        for (size_t j = 0; j < res.vega.cols(); ++j)
    {
        result->array[(i + 2)*resultCols + j] = res.vega[i][j];
    }

    // Return it
    return result;
}

extern "C" __declspec(dllexport)
    inline double xLinterp2D(
        FP12*               xs,
        FP12*               ys,
        FP12*               zs,
        double              x0,
        double              y0)
{
    vector<double> x;
    {
        size_t rows = xs->rows;
        size_t cols = xs->columns;
        double* numbers = xs->array;

        x.resize(rows);
        copy(numbers, numbers + rows, x.begin());
    }

    vector<double> y;
    {
        size_t rows = ys->rows;
        size_t cols = ys->columns;
        double* numbers = ys->array;

        y.resize(cols);
        copy(numbers, numbers + cols, y.begin());
    }

    matrix<double> z;
    {
        size_t rows = zs->rows;
        size_t cols = zs->columns;
        double* numbers = zs->array;

        z.resize(rows, cols);
        copy(numbers, numbers + rows * cols, z.begin());
    }

    //  Call 
    return interp2D(x, y, z, x0, y0);
}

extern "C" __declspec(dllexport)
inline FP12* xDupireCalib(
    //  model parameters
    const double ivsType, //  0: Bach, 1: BS, 2: Merton
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
        ivsType < 0.5 ? 'B' : ivsType < 1.5 ? 'S' : 'M',
        spot, vol, jmpIntens, jmpAverage, jmpStd);

    //  Build return

    // Allocate result
    // Calculate size
    size_t resultRows = results.spots.size() + 1, resultCols = results.times.size() + 1, resultSize = resultRows * resultCols;

    // Return an error if size is 0
    if (resultSize <= 0) return nullptr;

    // First, free all memory previously allocated
    // the function is defined in framework.h
    FreeAllTempMemory();
    // Then, allocate the memory for this result
    // We don't need to de-allocate it, that will be done by the next call
    // to this function or another function calling FreeAllTempMemory()
    // Memory size required, details in Dalton's book, section 6.2.2
    size_t memSize = sizeof(FP12) + (resultSize - 1) * sizeof(double);
    // Finally allocate, function definition in framework.h
    FP12* result = (FP12*)GetTempMemory(memSize);
    // Compute result
    result->rows = resultRows;
    result->columns = resultCols;
    result->array[0] = 0.0;
    for (size_t j = 0; j < results.times.size(); ++j)
    {
        result->array[j + 1] = results.times[j];
    }
    for (size_t i = 0; i < results.spots.size(); ++i)
    {
        result->array[(i + 1)*resultCols] = results.spots[i];
        for (size_t j = 0; j < results.times.size(); ++j)
        {
            result->array[(i + 1)*resultCols + j + 1] = results.lVols[i][j];
        }
    }

    // Return it
    return result;
}

extern "C" __declspec(dllexport)
inline FP12* xDupireSuperbucket(
    //  model parameters
    const double ivsType, //  0: Bach, 1: BS, 2: Merton
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
    const double maxDs,
    FP12* times,
    const double maxDtVol,
    //  MC
    double              maxDtSimul,
    //  product parameters
    double              strike,
    double              barrier,
    double              maturity,
    double              monitorFreq,
    //  numerical parameters
    double              useSobol,
    double              useAnti,
    double              seed1,
    double              seed2,
    double              numPath,
    double              parallel)

{
    //  Make sure the last input is given
    if (numPath <= 0) return 0;

    //  Unpack
    vector<double> vstrikes;
    {
        size_t rows = strikes->rows;
        size_t cols = strikes->columns;
        double* numbers = strikes->array;

        vstrikes.resize(rows*cols);
        copy(numbers, numbers + rows*cols, vstrikes.begin());
    }

    vector<double> vmats;
    {
        size_t rows = mats->rows;
        size_t cols = mats->columns;
        double* numbers = mats->array;

        vmats.resize(rows*cols);
        copy(numbers, numbers + rows*cols, vmats.begin());
    }

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

    auto results = dupireSuperbucket(
        spot,
        strike,
        barrier,
        maturity,
        monitorFreq,
        maxDtSimul,
        parallel>0,
        useSobol > 0,
        static_cast<int>(numPath), 
        useAnti > 0,
        vspots,
        maxDs,
        vtimes,
        maxDtVol,
        vstrikes,
        vmats, 
        ivsType < 0.5 ? 'B' : ivsType < 1.5 ? 'S' : 'M',
        vol,
        jmpIntens,
        jmpAverage,
        jmpStd,
        static_cast<int>(seed1), 
        static_cast<int>(seed2));

    //  Build return

    // Allocate result
    // Calculate size
    size_t resultRows = results.vega.rows() + 3, resultCols = results.vega.cols() + 1, resultSize = resultRows * resultCols;

    // Return an error if size is 0
    if (resultSize <= 0) return nullptr;

    // First, free all memory previously allocated
    // the function is defined in framework.h
    FreeAllTempMemory();
    // Then, allocate the memory for this result
    // We don't need to de-allocate it, that will be done by the next call
    // to this function or another function calling FreeAllTempMemory()
    // Memory size required, details in Dalton's book, section 6.2.2
    size_t memSize = sizeof(FP12) + (resultSize - 1) * sizeof(double);
    // Finally allocate, function definition in framework.h
    FP12* result = (FP12*)GetTempMemory(memSize);
    // Compute result
    result->rows = resultRows;
    result->columns = resultCols;
    for (size_t i = 0; i < resultSize; ++i) result->array[i] = 0.0;
    result->array[0] = results.value;
    result->array[resultCols] = results.delta;
    for (size_t j = 0; j < results.vega.cols(); ++j)
    {
        result->array[2*resultCols + j + 1] = vmats[j];
    }
    for (size_t i = 0; i < results.vega.rows(); ++i)
    {
        result->array[(3 + i) * resultCols] = vstrikes[i];
    }
    for (size_t i = 0; i < results.vega.rows(); ++i) 
        for (size_t j = 0; j < results.vega.cols(); ++j)
    {
        result->array[(3 + i)*resultCols + j + 1] = results.vega[i][j];
    }

    // Return it
    return result;
}

extern "C" __declspec(dllexport)
double xMerton(double spot, double vol, double mat, double strike, double intens, double meanJmp, double stdJmp)
{
    return merton(spot, strike, vol, mat, intens, meanJmp, stdJmp);
}

//	Registers

extern "C" __declspec(dllexport) int xlAutoOpen(void)
{
	XLOPER12 xDLL;

	Excel12f(xlGetName, &xDLL, 0);

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xUocDupire"),
		(LPXLOPER12)TempStr12(L"BBK%K%K%BBBBBBBBBBB"),
		(LPXLOPER12)TempStr12(L"xUocDupire"),
		(LPXLOPER12)TempStr12(L"spot, spots, times, vols, maxDt, K, B, T, freq, useSobol, useAnti, [seed1], [seed2], N, [Parallel]"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Computes the value of an up and out call in Dupire"),
		(LPXLOPER12)TempStr12(L"number 1, number 2"));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xLinterp2D"),
        (LPXLOPER12)TempStr12(L"BK%K%K%BB"),
        (LPXLOPER12)TempStr12(L"xLinterp2D"),
        (LPXLOPER12)TempStr12(L"xs, ys, zs, x0, y0"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Linear interpolation"),
        (LPXLOPER12)TempStr12(L"number 1, number 2"));


    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xUocDupireBump"),
        (LPXLOPER12)TempStr12(L"K%BK%K%K%BBBBBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xUocDupireBump"),
        (LPXLOPER12)TempStr12(L"spot, spots, times, vols, maxDt, K, B, T, freq, useSobol, useAnti, [seed1], [seed2], N, [Parallel]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Computes the risk of an up and out call in Dupire"),
        (LPXLOPER12)TempStr12(L"number 1, number 2"));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xUocDupireAAD"),
        (LPXLOPER12)TempStr12(L"K%BK%K%K%BBBBBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xUocDupireAAD"),
        (LPXLOPER12)TempStr12(L"spot, spots, times, vols, maxDt, K, B, T, freq, useSobol, useAnti, [seed1], [seed2], N, [Parallel]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Computes the risk of an up and out call in Dupire"),
        (LPXLOPER12)TempStr12(L"number 1, number 2"));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xDupireSuperbucket"),
        (LPXLOPER12)TempStr12(L"K%BBBBBBK%K%K%BK%BBBBBBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xDupireSuperbucket"),
        (LPXLOPER12)TempStr12(L"ivs, spot, vol, jmpIt, jmpAve, jmpStd, RiskStrikes, riskMats, volSpots, maxDs, volTimes, maxDtVol, maxDtSimul, K, B, T, Bfreq, sobol, anti, s1, s2, numPth, parallel"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Computes the risk of an up and out call in Dupire"),
        (LPXLOPER12)TempStr12(L"number 1, number 2"));

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xUocDupireSuperbucket"),
        (LPXLOPER12)TempStr12(L"K%BK%K%K%BBBBBBBBBBB"),
        (LPXLOPER12)TempStr12(L"xUocDupireSuperbucket"),
        (LPXLOPER12)TempStr12(L"spot, strikes, mats, calls, maxDt, K, B, T, freq, useSobol, useAnti, [seed1], [seed2], N, [Parallel]"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Computes the risk of an up and out call in Dupire"),
        (LPXLOPER12)TempStr12(L"number 1, number 2"));

/*
    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xDupireCalib"),
        (LPXLOPER12)TempStr12(L"K%BK%K%K%"),
        (LPXLOPER12)TempStr12(L"xDupireCalib"),
        (LPXLOPER12)TempStr12(L"spot, strikes, mats, calls"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Calibrates Dupire"),
        (LPXLOPER12)TempStr12(L"number 1, number 2"));

        */

    Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
        (LPXLOPER12)TempStr12(L"xDupireCalib"),
        (LPXLOPER12)TempStr12(L"K%BBBBBBK%BK%B"),
        (LPXLOPER12)TempStr12(L"xDupireCalib"),
        (LPXLOPER12)TempStr12(L"type, spot, vol, jumpIntensity, jumpAverage, jumpStd, spots, maxds, times, mxdt"),
        (LPXLOPER12)TempStr12(L"1"),
        (LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L""),
        (LPXLOPER12)TempStr12(L"Calibrates Dupire"),
        (LPXLOPER12)TempStr12(L"number 1, number 2"));

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
        (LPXLOPER12)TempStr12(L"number 1, number 2"));

	/* Free the XLL filename */
	Excel12f(xlFree, 0, 1, (LPXLOPER12)&xDLL);

    /*  Start the thread pool   */
    ThreadPool::getInstance()->start();

	return 1;
}

extern "C" __declspec(dllexport) int xlAutoClose(void)
{
    /*  Stop the thread pool   */
    ThreadPool::getInstance()->stop();

    return 1;
}