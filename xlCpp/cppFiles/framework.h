///***************************************************************************
// File:        FRAMEWRK.H
//
// Purpose:     Header file for Framework library
//
// Platform:    Microsoft Windows
//
// Comments:
//              Include this file in any source files
//              that use the framework library.
//
// From the Microsoft Excel Developer's Kit, Version 12
// Copyright (c) 1997 - 2007 Microsoft Corporation. All rights reserved.
///***************************************************************************

#pragma once

#include <wchar.h>

#include "memorymanager.h"

extern "C" LPSTR cdecl GetTempMemory(int cBytes)
{
	return MGetTempMemory(cBytes);
}

extern "C" void cdecl FreeAllTempMemory(void)
{
	MFreeAllTempMemory();
}

extern "C" int cdecl Excel12f(int xlfn, LPXLOPER12 pxResult, int count, ...)
{
	int xlret = Excel12v(xlfn, pxResult, count, (LPXLOPER12 *)(&count + 1));
	FreeAllTempMemory();
	return xlret;
}

extern "C" LPXLOPER12 TempStr12(const XCHAR* lpstr)
{
	LPXLOPER12 lpx;
	XCHAR* lps;
	int len;

	len = lstrlenW(lpstr);

	lpx = (LPXLOPER12)GetTempMemory(sizeof(XLOPER12) + (len + 1) * 2);

	if (!lpx)
	{
		return 0;
	}

	lps = (XCHAR*)((CHAR*)lpx + sizeof(XLOPER12));

	lps[0] = (BYTE)len;
	//can't wcscpy_s because of removal of null-termination
	wmemcpy_s(lps + 1, len + 1, lpstr, len);
	lpx->xltype = xltypeStr;
	lpx->val.str = lps;

	return lpx;
}

extern "C" LPXLOPER12 TempErr12(int err)
{
    LPXLOPER12 lpx;

    lpx = (LPXLOPER12)GetTempMemory(sizeof(XLOPER12));

    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeErr;
    lpx->val.err = err;

    return lpx;
}

LPXLOPER12 TempNum12(double d)
{
    LPXLOPER12 lpx;

    lpx = (LPXLOPER12)GetTempMemory(sizeof(XLOPER12));

    if (!lpx)
    {
        return 0;
    }

    lpx->xltype = xltypeNum;
    lpx->val.num = d;

    return lpx;
}