
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

	/* Free the XLL filename */
	Excel12f(xlFree, 0, 1, (LPXLOPER12)&xDLL);

    /*  Start the thread pool   */
    ThreadPool::getInstance()->start();

	return 1;
}
