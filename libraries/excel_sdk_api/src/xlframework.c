///***************************************************************************
// File:	FRAMEWRK.C
//
// Purpose:	Framework library for Microsoft Excel.
//
//     This library provides some basic functions
//     that help in writing Excel DLLs. It includes
//     simple functions for managing memory with XLOPER12s,
//     creating temporary XLOPER12s, robustly calling
//     Excel12(), and outputting debugging strings
//     to the debugger for the current application.
//   
//     The main purpose of this library is to help
//     you to write cleaner C code for calling Excel.
//     For example, using the framework library you
//     can write
//   
//         Excel12f(xlcDisplay, 0, 2, TempMissing12(), TempBool12(0));
//   
//     instead of the more verbose
//   
//         XLOPER12 xMissing, xBool;
//         xMissing.xltype = xltypeMissing;
//         xBool.xltype = xltypeBool;
//         xBool.val.xbool = 0;
//         Excel12(xlcDisplay, 0, 2, (LPXLOPER12) &xMissing, (LPXLOPER12) &xBool);
//   
//   
//     The library is non-reentrant.
//   
//     Define _DEBUG to use the debugging functions.
//   
//     Source code is provided so that you may
//     enhance this library or optimize it for your
//     own application.
// 
// Platform:    Microsoft Windows
//
// Functions:		
//              debugPrintf
//              GetTempMemory
//              FreeAllTempMemory
//              Excel
//              Excel12f
//              TempNum
//              TempNum12
//              TempStr
//              TempStrConst
//              TempStr12
//              TempBool
//              TempBool12
//              TempInt
//              TempInt12
//              TempErr
//              TempErr12
//              TempActiveRef
//              TempActiveRef12
//              TempActiveCell
//              TempActiveCell12
//              TempActiveRow
//              TempActiveRow12
//              TempActiveColumn
//              TempActiveColumn12
//              TempMissing
//              TempMissing12
//              InitFramework
//				QuitFramework
//
//***************************************************************************


#include <windows.h>
#include <malloc.h>
#include <wchar.h>
#include "xlcall.h"
#include "xlframework.h"
#include "xlMemoryManager.h"
#include <stdarg.h>


#ifdef _DEBUG

void cdecl debugPrintf(LPSTR lpFormat, ...)
{
	char rgch[256];
	va_list argList;

	va_start(argList, lpFormat);
	wvsprintfA(rgch, lpFormat, argList);
	va_end(argList);
	OutputDebugStringA(rgch);
}

#endif /* DEBUG */

LPSTR GetTempMemory(size_t cBytes)
{
	return MGetTempMemory(cBytes);
}

void FreeAllTempMemory(void)
{
	MFreeAllTempMemory();
}

int cdecl Excel12f(int xlfn, LPXLOPER12 pxResult, int count, ...)
{
	int xlret;
	va_list ppxArgs;

#ifdef _DEBUG
	LPXLOPER12 px;
	int i;

	va_start(ppxArgs, count);

	for (i = 0; i < count; i++)
	{
		px = va_arg(ppxArgs, LPXLOPER12);

		if (px == NULL)
		{
			debugPrintf("Error! A zero parameter is being passed "
				"to Excel12f. Check memory...\r");
			FreeAllTempMemory();
			return xlretFailed;
		}
	}

	va_end(ppxArgs);

#endif /* DEBUG */

	va_start(ppxArgs, count);
	xlret = Excel12v(xlfn, pxResult, count, (LPXLOPER12*)ppxArgs);
	va_end(ppxArgs);

#ifdef _DEBUG

	if (xlret != xlretSuccess)
	{
		debugPrintf("Error! Excel12(");

		if (xlfn & xlCommand)
			debugPrintf("xlCommand | ");
		if (xlfn & xlSpecial)
			debugPrintf("xlSpecial | ");
		if (xlfn & xlIntl)
			debugPrintf("xlIntl | ");
		if (xlfn & xlPrompt)
			debugPrintf("xlPrompt | ");

		debugPrintf("%u) callback failed:", xlfn & 0x0FFF);

		/* More than one error bit may be on */

		if (xlret & xlretAbort)
		{
			debugPrintf(" Macro Halted\r");
		}

		if (xlret & xlretInvXlfn)
		{
			debugPrintf(" Invalid Function Number\r");
		}

		if (xlret & xlretInvCount)
		{
			debugPrintf(" Invalid Number of Arguments\r");
		}

		if (xlret & xlretInvXloper)
		{
			debugPrintf(" Invalid XLOPER12\r");
		}

		if (xlret & xlretStackOvfl)
		{
			debugPrintf(" Stack Overflow\r");
		}

		if (xlret & xlretFailed)
		{
			debugPrintf(" Command failed\r");
		}

		if (xlret & xlretUncalced)
		{
			debugPrintf(" Uncalced cell\r");
		}

	}

#endif /* DEBUG */

	FreeAllTempMemory();

	return xlret;
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

LPXLOPER12 TempStr12(const XCHAR* lpstr)
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

LPXLOPER12 TempBool12(BOOL b)
{
	LPXLOPER12 lpx;

	lpx = (LPXLOPER12)GetTempMemory(sizeof(XLOPER12));

	if (!lpx)
	{
		return 0;
	}

	lpx->xltype = xltypeBool;
	lpx->val.xbool = b ? 1 : 0;

	return lpx;
}

LPXLOPER12 TempInt12(int i)
{
	LPXLOPER12 lpx;

	lpx = (LPXLOPER12)GetTempMemory(sizeof(XLOPER12));

	if (!lpx)
	{
		return 0;
	}

	lpx->xltype = xltypeInt;
	lpx->val.w = i;

	return lpx;
}

LPXLOPER12 TempErr12(int err)
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

LPXLOPER12 TempActiveRef12(RW rwFirst, RW rwLast, COL colFirst, COL colLast)
{
	LPXLOPER12 lpx;
	LPXLMREF12 lpmref;
	int wRet;

	lpx = (LPXLOPER12)GetTempMemory(sizeof(XLOPER12));
	lpmref = (LPXLMREF12)GetTempMemory(sizeof(XLMREF12));

	if (!lpmref)
	{
		return 0;
	}

	/* calling Excel12f() instead of Excel12() would free all temp memory! */
	wRet = Excel12(xlSheetId, lpx, 0);

	if (wRet != xlretSuccess)
	{
		return 0;
	}
	else
	{
		lpx->xltype = xltypeRef;
		lpx->val.mref.lpmref = lpmref;
		lpmref->count = 1;
		lpmref->reftbl[0].rwFirst = rwFirst;
		lpmref->reftbl[0].rwLast = rwLast;
		lpmref->reftbl[0].colFirst = colFirst;
		lpmref->reftbl[0].colLast = colLast;

		return lpx;
	}
}

LPXLOPER12 TempActiveCell12(RW rw, COL col)
{
	return TempActiveRef12(rw, rw, col, col);
}

LPXLOPER12 TempActiveRow12(RW rw)
{
	return TempActiveRef12(rw, rw, 0, 0x00003FFF);
}

LPXLOPER12 TempActiveColumn12(COL col)
{
	return TempActiveRef12(0, 0x000FFFFF, col, col);
}

LPXLOPER12 TempMissing12(void)
{
	LPXLOPER12 lpx;

	lpx = (LPXLOPER12)GetTempMemory(sizeof(XLOPER12));

	if (!lpx)
	{
		return 0;
	}

	lpx->xltype = xltypeMissing;

	return lpx;
}

void FreeXLOperT(LPXLOPER pxloper)
{
	WORD xltype;
	int cxloper;
	LPXLOPER pxloperFree;

	xltype = pxloper->xltype;

	switch (xltype)
	{
	case xltypeStr:
		if (pxloper->val.str != NULL)
			free(pxloper->val.str);
		break;
	case xltypeRef:
		if (pxloper->val.mref.lpmref != NULL)
			free(pxloper->val.mref.lpmref);
		break;
	case xltypeMulti:
		cxloper = pxloper->val.array.rows * pxloper->val.array.columns;
		if (pxloper->val.array.lparray != NULL)
		{
			pxloperFree = pxloper->val.array.lparray;
			while (cxloper > 0)
			{
				FreeXLOperT(pxloperFree);
				pxloperFree++;
				cxloper--;
			}
			free(pxloper->val.array.lparray);
		}
		break;
	case xltypeBigData:
		if (pxloper->val.bigdata.h.lpbData != NULL)
			free(pxloper->val.bigdata.h.lpbData);
		break;
	}
}


void FreeXLOper12T(LPXLOPER12 pxloper12)
{
	DWORD xltype;
	int cxloper12;
	LPXLOPER12 pxloper12Free;

	xltype = pxloper12->xltype;

	switch (xltype)
	{
	case xltypeStr:
		if (pxloper12->val.str != NULL)
			free(pxloper12->val.str);
		break;
	case xltypeRef:
		if (pxloper12->val.mref.lpmref != NULL)
			free(pxloper12->val.mref.lpmref);
		break;
	case xltypeMulti:
		cxloper12 = pxloper12->val.array.rows * pxloper12->val.array.columns;
		if (pxloper12->val.array.lparray)
		{
			pxloper12Free = pxloper12->val.array.lparray;
			while (cxloper12 > 0)
			{
				FreeXLOper12T(pxloper12Free);
				pxloper12Free++;
				cxloper12--;
			}
			free(pxloper12->val.array.lparray);
		}
		break;
	case xltypeBigData:
		if (pxloper12->val.bigdata.h.lpbData != NULL)
			free(pxloper12->val.bigdata.h.lpbData);
		break;
	}
}


BOOL ConvertXLRefToXLRef12(LPXLREF pxref, LPXLREF12 pxref12)
{
	if (pxref->rwLast >= pxref->rwFirst && pxref->colLast >= pxref->colFirst)
	{
		if (pxref->rwFirst >= 0 && pxref->colFirst >= 0)
		{
			pxref12->rwFirst = pxref->rwFirst;
			pxref12->rwLast = pxref->rwLast;
			pxref12->colFirst = pxref->colFirst;
			pxref12->colLast = pxref->colLast;
			return TRUE;
		}
	}
	return FALSE;
}

BOOL ConvertXLRef12ToXLRef(LPXLREF12 pxref12, LPXLREF pxref)
{
	if (pxref12->rwLast >= pxref12->rwFirst && pxref12->colLast >= pxref12->colFirst)
	{
		if (pxref12->rwFirst >= 0 && pxref12->colFirst >= 0)
		{
			if (pxref12->rwLast < rwMaxO8 && pxref12->colLast < colMaxO8)
			{
				pxref->rwFirst = (WORD)pxref12->rwFirst;
				pxref->rwLast = (WORD)pxref12->rwLast;
				pxref->colFirst = (BYTE)pxref12->colFirst;
				pxref->colLast = (BYTE)pxref12->colLast;
				return TRUE;
			}
		}
	}
	return FALSE;
}

BOOL XLOper12ToXLOper(LPXLOPER12 pxloper12, LPXLOPER pxloper)
{
	BOOL fRet;
	BOOL fClean;
	DWORD xltype;
	WORD cref;
	int cxloper12;
	RW crw;
	COL ccol;
	long cbyte;
	XCHAR* st;
	char* ast;
	int cch;
	char cach;
	BYTE* pbyte;
	LPXLMREF pmref;
	LPXLREF12 pref12;
	LPXLREF rgref;
	LPXLREF pref;
	LPXLOPER rgxloperConv;
	LPXLOPER pxloperConv;
	LPXLOPER12 pxloper12Conv;

	fClean = FALSE;
	fRet = TRUE;
	xltype = pxloper12->xltype;

	switch (xltype)
	{
	case xltypeNum:
		pxloper->val.num = pxloper12->val.num;
		break;
	case xltypeBool:
#ifdef __cplusplus
		pxloper->val.xbool = pxloper12->val.xbool;
#else		
		pxloper->val.bool = pxloper12->val.xbool;
#endif
		break;
	case xltypeErr:
		if (pxloper12->val.err > MAXWORD)
		{
			fRet = FALSE;
			// problem... overflow
		}
		else
		{
			pxloper->val.err = (WORD)pxloper12->val.err;
		}
		break;
	case xltypeMissing:
	case xltypeNil:
		break;
	case xltypeInt:
		if ((pxloper12->val.w + MAXSHORTINT + 1) >> 16)
		{
			pxloper->val.num = (float)pxloper12->val.w;
			xltype = xltypeNum;
		}
		else
		{
			pxloper->val.w = (short int)pxloper12->val.w;
		}
		break;
	case xltypeStr:
		st = pxloper12->val.str;
		cch = st[0];
		cach = (BYTE)cch;

		if (cch > cchMaxStz || cch < 0)
		{
			fRet = FALSE;
		}
		else
		{
			ast = malloc((cach + 2) * sizeof(char));
			if (ast == NULL)
			{
				fRet = FALSE;
			}
			else
			{
				WideCharToMultiByte(CP_ACP, 0, st + 1, cch, ast + 1, cach, NULL, NULL);
				ast[0] = cach;
				ast[cach + 1] = '\0';
				pxloper->val.str = ast;
			}
		}
		break;
	case xltypeFlow:
		if (pxloper12->val.flow.rw > rwMaxO8 || pxloper12->val.flow.col > colMaxO8)
		{
			fRet = FALSE;
		}
		else
		{
			pxloper->val.flow.rw = (WORD)pxloper12->val.flow.rw;
			pxloper->val.flow.col = (BYTE)pxloper12->val.flow.col;
			pxloper->val.flow.xlflow = pxloper12->val.flow.xlflow;
			pxloper->val.flow.valflow.idSheet = pxloper12->val.flow.valflow.idSheet;
		}
		break;
	case xltypeRef:
		if (pxloper12->val.mref.lpmref && pxloper12->val.mref.lpmref->count > 0)
		{
			pref12 = pxloper12->val.mref.lpmref->reftbl;
			cref = pxloper12->val.mref.lpmref->count;

			pmref = (LPXLMREF)malloc(sizeof(XLMREF) + sizeof(XLREF) * (cref - 1));
			if (pmref == NULL)
			{
				fRet = FALSE;
			}
			else
			{
				pmref->count = cref;
				rgref = pmref->reftbl;
				pref = rgref;
				while (cref > 0 && !fClean)
				{
					if (!ConvertXLRef12ToXLRef(pref12, pref))
					{
						fClean = TRUE;
						cref = 0;
					}
					else
					{
						pref++;
						pref12++;
						cref--;
					}
				}
				if (fClean)
				{
					free(pmref);
					fRet = FALSE;
				}
				else
				{
					pxloper->val.mref.lpmref = pmref;
					pxloper->val.mref.idSheet = pxloper12->val.mref.idSheet;
				}
			}
		}
		else
		{
			xltype = xltypeMissing;
		}
		break;
	case xltypeSRef:
		if (pxloper12->val.sref.count != 1)
		{
			fRet = FALSE;
		}
		else if (ConvertXLRef12ToXLRef(&pxloper12->val.sref.ref, &pxloper->val.sref.ref))
		{
			pxloper->val.sref.count = 1;
		}
		else
		{
			fRet = FALSE;
		}
		break;
	case xltypeMulti:
		crw = pxloper12->val.array.rows;
		ccol = pxloper12->val.array.columns;
		if (crw > rwMaxO8 || ccol > colMaxO8)
		{
			fRet = FALSE;
		}
		else
		{
			cxloper12 = crw * ccol;
			if (cxloper12 == 0)
			{
				xltype = xltypeMissing;
			}
			else
			{
				rgxloperConv = malloc(cxloper12 * sizeof(XLOPER));
				if (rgxloperConv == NULL)
				{
					fRet = FALSE;
				}
				else
				{
					pxloperConv = rgxloperConv;
					pxloper12Conv = pxloper12->val.array.lparray;
					while (cxloper12 > 0 && !fClean)
					{
						if (!XLOper12ToXLOper(pxloper12Conv, pxloperConv))
						{
							fClean = TRUE;
							cxloper12 = 0;
						}
						else
						{
							pxloperConv++;
							pxloper12Conv++;
							cxloper12--;
						}
					}
					if (fClean)
					{
						fRet = FALSE;
						while (pxloperConv > rgxloperConv)
						{
							FreeXLOperT(pxloperConv);
							pxloperConv--;
						}
						free(rgxloperConv);
					}
					else
					{
						pxloper->val.array.lparray = rgxloperConv;
						pxloper->val.array.rows = crw;
						pxloper->val.array.columns = ccol;
					}
				}
			}
		}
		break;
	case xltypeBigData:
		cbyte = pxloper12->val.bigdata.cbData;
		if (pxloper12->val.bigdata.h.lpbData != NULL && cbyte > 0)
		{
			pbyte = (BYTE*)malloc(cbyte);
			if (pbyte != NULL)
			{
				memcpy_s(pbyte, cbyte, pxloper12->val.bigdata.h.lpbData, cbyte);
				pxloper->val.bigdata.h.lpbData = pbyte;
				pxloper->val.bigdata.cbData = cbyte;
			}
			else
			{
				fRet = FALSE;
			}
		}
		else
		{
			fRet = FALSE;
		}
		break;
	}
	if (fRet)
	{
		pxloper->xltype = (WORD)xltype;
	}
	return fRet;
}

BOOL XLOperToXLOper12(LPXLOPER pxloper, LPXLOPER12 pxloper12)
{
	BOOL fRet;
	BOOL fClean;
	WORD crw;
	WORD ccol;
	WORD cxloper;
	WORD cref12;
	DWORD xltype;
	BYTE* pbyte;
	const char* ast;
	XCHAR* st;
	int cach;
	int cch;
	long cbyte;
	LPXLREF pref;
	LPXLREF12 pref12;
	LPXLREF12 rgref12;
	LPXLMREF12 pmref12;
	LPXLOPER pxloperConv;
	LPXLOPER12 rgxloper12Conv;
	LPXLOPER12 pxloper12Conv;

	fClean = FALSE;
	fRet = TRUE;
	xltype = pxloper->xltype;

	switch (xltype)
	{
	case xltypeNum:
		pxloper12->val.num = pxloper->val.num;
		break;
	case xltypeBool:
#ifdef __cplusplus
		pxloper12->val.xbool = pxloper->val.xbool;
#else		
		pxloper12->val.xbool = pxloper->val.bool;
#endif
		break;
	case xltypeErr:
		pxloper12->val.err = (int)pxloper->val.err;
		break;
	case xltypeMissing:
	case xltypeNil:
		break;
	case xltypeInt:
		pxloper12->val.w = pxloper->val.w;
	case xltypeStr:
		ast = pxloper->val.str;
		if (ast == NULL)
		{
			fRet = FALSE;
		}
		else
		{
			cach = ast[0];
			cch = cach;
			if (cach > cchMaxStz || cach < 0)
			{
				fRet = FALSE;
			}
			else
			{
				st = (XCHAR*)malloc((cch + 2) * sizeof(XCHAR));
				if (st == NULL)
				{
					fRet = FALSE;
				}
				else
				{
					MultiByteToWideChar(CP_ACP, 0, ast + 1, cach, st + 1, cch);
					st[0] = (XCHAR)cch;
					st[cch + 1] = '\0';
					pxloper12->val.str = st;
				}
			}
		}
		break;
	case xltypeFlow:
		pxloper12->val.flow.rw = pxloper->val.flow.rw;
		pxloper12->val.flow.col = pxloper->val.flow.col;
		pxloper12->val.flow.xlflow = pxloper->val.flow.xlflow;
		break;
	case xltypeRef:
		if (pxloper->val.mref.lpmref && pxloper->val.mref.lpmref->count > 0)
		{
			pref = pxloper->val.mref.lpmref->reftbl;
			cref12 = pxloper->val.mref.lpmref->count;

			pmref12 = (LPXLMREF12)malloc(sizeof(XLMREF12) + sizeof(XLREF12) * (cref12 - 1));
			if (pmref12 == NULL)
			{
				fRet = FALSE;
			}
			else
			{
				pmref12->count = cref12;
				rgref12 = pmref12->reftbl;
				pref12 = rgref12;
				while (cref12 > 0 && !fClean)
				{
					if (!ConvertXLRefToXLRef12(pref, pref12))
					{
						fClean = TRUE;
						cref12 = 0;
					}
					else
					{
						pref++;
						pref12++;
						cref12--;
					}
				}
				if (fClean)
				{
					free(pmref12);
					fRet = FALSE;
				}
				else
				{
					pxloper12->val.mref.lpmref = pmref12;
					pxloper12->val.mref.idSheet = pxloper->val.mref.idSheet;
				}
			}
		}
		else
		{
			xltype = xltypeMissing;
		}
		break;
	case xltypeSRef:
		if (pxloper->val.sref.count != 1)
		{
			fRet = FALSE;
		}
		else if (ConvertXLRefToXLRef12(&pxloper->val.sref.ref, &pxloper12->val.sref.ref))
		{
			pxloper12->val.sref.count = 1;
		}
		else
		{
			fRet = FALSE;
		}
		break;
	case xltypeMulti:
		crw = pxloper->val.array.rows;
		ccol = pxloper->val.array.columns;
		if (crw > rwMaxO8 || ccol > colMaxO8)
		{
			fRet = FALSE;
		}
		else
		{
			cxloper = crw * ccol;
			if (cxloper == 0)
			{
				xltype = xltypeMissing;
			}
			else
			{
				rgxloper12Conv = malloc(cxloper * sizeof(XLOPER12));
				if (rgxloper12Conv == NULL)
				{
					fRet = FALSE;
				}
				else
				{
					pxloper12Conv = rgxloper12Conv;
					pxloperConv = pxloper->val.array.lparray;
					while (cxloper > 0 && !fClean)
					{
						if (!XLOperToXLOper12(pxloperConv, pxloper12Conv))
						{
							fClean = TRUE;
							cxloper = 0;
						}
						else
						{
							pxloperConv++;
							pxloper12Conv++;
							cxloper--;
						}
					}
					if (fClean)
					{
						fRet = FALSE;
						while (pxloper12Conv > rgxloper12Conv)
						{
							FreeXLOperT(pxloperConv);
							pxloperConv--;
						}
						free(rgxloper12Conv);
					}
					else
					{
						pxloper12->val.array.lparray = rgxloper12Conv;
						pxloper12->val.array.rows = crw;
						pxloper12->val.array.columns = ccol;
					}
				}
			}
		}
		break;
	case xltypeBigData:
		cbyte = pxloper->val.bigdata.cbData;
		if (pxloper->val.bigdata.h.lpbData != NULL && cbyte > 0)
		{
			pbyte = (BYTE*)malloc(cbyte);
			if (pbyte != NULL)
			{
				memcpy_s(pbyte, cbyte, pxloper->val.bigdata.h.lpbData, cbyte);
				pxloper12->val.bigdata.h.lpbData = pbyte;
				pxloper12->val.bigdata.cbData = cbyte;
			}
			else
			{
				fRet = FALSE;
			}
		}
		else
		{
			fRet = FALSE;
		}
		break;
	}

	if (fRet)
	{
		pxloper12->xltype = xltype;
	}
	return fRet;
}
