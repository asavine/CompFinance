#pragma once

#pragma warning(disable:4996)

//  General interface for ranges with strings

#include "xlcall.h"
#include "framework.h"
#include <string>
using namespace std;

//  Additional / alternative functions to the ones in framework.h

LPXLOPER12 TempXLOPER12()
{
    LPXLOPER12 lpx;

    lpx = (LPXLOPER12)GetTempMemory(sizeof(XLOPER12));

    if (!lpx)
    {
        return 0;
    }

    return lpx;
}

LPXLOPER12 TempStr12(const string str)
{
    LPXLOPER12 lpx = TempXLOPER12();
    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeStr;

    lpx->val.str = (wchar_t*)GetTempMemory((str.size() + 1) * sizeof(wchar_t));
    if (!lpx->val.str)
    {
        return 0;
    }

    if (str.size()>0) mbstowcs(lpx->val.str + 1, str.data(), str.size());
    lpx->val.str[0] = static_cast<wchar_t>(str.size());

    return lpx;
}

//  Getters

//  Get argument from Excel
//  Register as type Q
//  Get as LPXLOPER12

//  Number of rows or 0 if invalid
size_t getRows(const LPXLOPER12& oper)
{
    if (!oper) return 0;
    if (oper->xltype != xltypeMulti) return 1;
    return oper->val.array.rows;
}

//  Number of columns or 0 if invalid
size_t getCols(const LPXLOPER12& oper)
{
    if (!oper) return 0;
    if (oper->xltype != xltypeMulti) return 1;
    return oper->val.array.columns;
}

//  String in position (i, j) or "" if invalid
string getString(const LPXLOPER12& oper, const size_t i = 0, const size_t j = 0)
{
    if (!oper) return "";

    LPXLOPER12 strOper;

    if (oper->xltype == xltypeStr && i == 0 && j == 0)
    {
        strOper = oper;
    }
    else if (oper->xltype == xltypeMulti)
    {
        strOper = oper->val.array.lparray + i * getCols(oper) + j;
        if (strOper->xltype != xltypeStr) return "";
    }
    else
    {
        return "";
    }

    string str;
    str.resize(strOper->val.str[0]);
    wcstombs(&str[0], strOper->val.str + 1, strOper->val.str[0]);

    return str;
}

//  Number in position (i, j) or infinity if invalid
double getNum(const LPXLOPER12& oper, const size_t i = 0, const size_t j = 0)
{
    if (!oper) return numeric_limits<double>::infinity();
    if (oper->xltype == xltypeNum && i == 0 && j == 0)
    {
        return oper->val.num;
    }
    else if (oper->xltype == xltypeMulti)
    {
        XLOPER12 nOper = oper->val.array.lparray[i * getCols(oper) + j];
        if (nOper.xltype != xltypeNum) return numeric_limits<double>::infinity();
        return nOper.val.num;

    }
    else
    {
        return numeric_limits<double>::infinity();
    }
}

//  Setter to return result to Excel
//  Register as type Q
//  Return as LPXLOPER12
void resize(LPXLOPER12& oper, const size_t rows, const size_t cols)
{
    oper->xltype = xltypeMulti;
    oper->val.array.rows = rows;
    oper->val.array.columns = cols;
    oper->val.array.lparray = (LPXLOPER12)GetTempMemory(rows*cols*sizeof(XLOPER12));
    LPXLOPER12 nilOper = TempStr12("");
    fill(oper->val.array.lparray, oper->val.array.lparray + rows*cols, *nilOper);
}

//  Set string in position (i, j) 
void setString(LPXLOPER12& oper, const string& str, const size_t i = 0, const size_t j = 0)
{
    XLOPER12& strOper = oper->val.array.lparray[i * getCols(oper) + j];
    strOper = * TempStr12(str.c_str());
}

//  Number in position (i, j) or infinity if invalid
void setNum(LPXLOPER12& oper, const double num, const size_t i = 0, const size_t j = 0)
{
    XLOPER12& nOper = oper->val.array.lparray[i * getCols(oper) + j];
    nOper.xltype = xltypeNum;
    nOper.val.num = num;
}
