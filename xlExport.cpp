
#include <windows.h>
#include "xlcall.h"
#include "framework.h"

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
    double              numPath)
{
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

        vtimes.resize(rows);
        copy(numbers, numbers + rows, vtimes.begin());
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

    return uocDupire(spot, vspots, vtimes, vvols, maxDt, strike, barrier, maturity, monitorFreq, 
        useSobol > 0, static_cast<int>(numPath));

}

//	Registers

extern "C" __declspec(dllexport) int xlAutoOpen(void)
{
	XLOPER12 xDLL;

	Excel12f(xlGetName, &xDLL, 0);

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xUocDupire"),
		(LPXLOPER12)TempStr12(L"BBK%K%K%BBBBBBB"),
		(LPXLOPER12)TempStr12(L"xUocDupire"),
		(LPXLOPER12)TempStr12(L"spot, spots, times, vols, maxDt, K, B, T, freq, useSobol, N"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Computes the value of an up and out call in Dupire"),
		(LPXLOPER12)TempStr12(L"number 1, number 2"));

	/* Free the XLL filename */
	Excel12f(xlFree, 0, 1, (LPXLOPER12)&xDLL);

	return 1;
}
