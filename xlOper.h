#pragma once

#pragma warning(disable:4996)

//  General interface for ranges with strings

#include "xlcall.h"
#include "xlframework.h"
#include <string>
using namespace std;

#include "matrix.h"

//  Additional / alternative functions to the ones in xlframework.h

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

//  Various utilities

//  FP12* to vector<double>
vector<double> to_vector(const FP12* oper)
{
    size_t rows = oper->rows;
    size_t cols = oper->columns;
    const double* numbers = oper->array;

    vector<double> v(rows * cols);
    copy(numbers, numbers + rows * cols, v.begin());

    return v;
}

//  LPXLOPER to vector<string>
vector<string> to_strVector(const LPXLOPER12& oper)
{
    vector<string> vstr;

    if (oper->xltype == xltypeStr)
    {
        vstr.push_back(getString(oper));
    }
    else if (oper->xltype == xltypeMulti)
    {
        for (size_t i = 0; i < getRows(oper); ++i) for (size_t j = 0; j < getCols(oper); ++j)
        {
            vstr.push_back(getString(oper, i, j));
        }
    }

    return vstr;
}

//  FP12* to matrix<double>
matrix<double> to_matrix(const FP12* oper)
{
    size_t rows = oper->rows;
    size_t cols = oper->columns;
    const double* numbers = oper->array;

    matrix<double> m(rows, cols);
    copy(numbers, numbers + rows * cols, m.begin());

    return m;
}

LPXLOPER12 from_strVector(const vector<string>& strv, const bool rowVector = true)
{
    LPXLOPER12 oper = TempXLOPER12();
    oper->xltype = xltypeMulti;
    oper->val.array.rows = rowVector ? strv.size() : 1;
    oper->val.array.columns = rowVector ? 1 : strv.size();
    oper->val.array.lparray = (LPXLOPER12) GetTempMemory(strv.size() * sizeof(XLOPER12));
    transform(strv.begin(), strv.end(), oper->val.array.lparray,
        [](const string& s)
    {
        return *TempStr12(s);
    });

    return oper;
}

LPXLOPER12 from_labelsAndNumbers(const vector<string>& labels, const vector<double>& numbers)
{
    const size_t n = labels.size();
    if (n == 0 || n != numbers.size()) return TempErr12(xlerrNA);
    LPXLOPER12 oper = TempXLOPER12();
    oper->xltype = xltypeMulti;
    oper->val.array.rows = n;
    oper->val.array.columns = 2;
    oper->val.array.lparray = (LPXLOPER12)GetTempMemory(2 * n * sizeof(XLOPER12));

    for (size_t i = 0; i < n; ++i)
    {
        oper->val.array.lparray[2 * i] = *TempStr12(labels[i]);
        oper->val.array.lparray[2 * i + 1] = *TempNum12(numbers[i]);
    }

    return oper;
}

LPXLOPER12 from_labelledMatrix(const vector<string>& rowLabels, const vector<string>& colLabels, const matrix<double>& mat)
{
    const size_t n = rowLabels.size(), m = colLabels.size();
    if (n == 0 || m == 0 || n != mat.rows() || m != mat.cols()) return TempErr12(xlerrNA);

    LPXLOPER12 oper = TempXLOPER12();
    resize(oper, n + 1, m + 1);
    for (size_t i = 0; i < n; ++i) setString(oper, rowLabels[i], i + 1, 0);
    for (size_t i = 0; i < m; ++i) setString(oper, colLabels[i], 0, i + 1);

    for (size_t i = 0; i < n; ++i)  for (size_t j = 0; j < m; ++j) setNum(oper, mat[i][j], i + 1, j + 1);

    return oper;
}

LPXLOPER12 from_labelledMatrix(const vector<double>& rowLabels, const vector<double>& colLabels, const matrix<double>& mat)
{
    const size_t n = rowLabels.size(), m = colLabels.size();
    if (n == 0 || m == 0 || n != mat.rows() || m != mat.cols()) return TempErr12(xlerrNA);

    LPXLOPER12 oper = TempXLOPER12();
    resize(oper, n + 1, m + 1);
    for (size_t i = 0; i < n; ++i) setNum(oper, rowLabels[i], i + 1, 0);
    for (size_t i = 0; i < m; ++i) setNum(oper, colLabels[i], 0, i + 1);

    for (size_t i = 0; i < n; ++i)  for (size_t j = 0; j < m; ++j) setNum(oper, mat[i][j], i + 1, j + 1);

    return oper;
}

LPXLOPER12 from_labelledMatrix(const vector<string>& rowLabels, const vector<string>& colLabels, const matrix<double>& mat,
	const string& firstLineLabel, const vector<double>& firstLine)
{
	const size_t n = rowLabels.size(), m = colLabels.size();
	if (n == 0 || m == 0 || n != mat.rows() || m != mat.cols() || m != firstLine.size()) return TempErr12(xlerrNA);

	LPXLOPER12 oper = TempXLOPER12();
	resize(oper, n + 2, m + 1);
	setString(oper, firstLineLabel, 1, 0);
	for (size_t i = 0; i < m; ++i) setString(oper, colLabels[i], 0, i + 1);
	for (size_t i = 0; i < m; ++i) setNum(oper, firstLine[i], 1, i + 1);
	for (size_t i = 0; i < n; ++i) setString(oper, rowLabels[i], i + 2, 0);
	for (size_t i = 0; i < n; ++i)  for (size_t j = 0; j < m; ++j) setNum(oper, mat[i][j], i + 2, j + 1);

	return oper;
}