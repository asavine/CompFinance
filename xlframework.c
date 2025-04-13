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
//#include "xlcall.cpp"
#include "xlframework.h"
#include "xlMemorymanager.h"
#include <stdarg.h>


#ifdef _DEBUG

///***************************************************************************
// debugPrintf()
//
// Purpose: 
//      sends a string to the debugger for the current application. 
//
// Parameters:
//
//      LPSTR lpFormat  The format definition string
//      ...             The values to print
//
// Returns: 
//
// Comments:
//
///***************************************************************************

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

///***************************************************************************
// GetTempMemory()
//
// Purpose: 
//       Allocates temporary memory. Temporary memory
//       can only be freed in one chunk, by calling
//       FreeAllTempMemory(). This is done by Excel12f().
//
// Parameters:
//
//      size_t cBytes      How many bytes to allocate
//
// Returns: 
//
//      LPSTR           A pointer to the allocated memory,
//                      or 0 if more memory cannot be
//                      allocated. If this fails,
//                      check that MEMORYSIZE is big enough
//                      in MemoryPool.h and that you have
//                      remembered to call FreeAllTempMemory
//
// Comments:
//
//	Algorithm:
//       
//      The memory allocation algorithm is extremely
//      simple: on each call, allocate the next cBytes
//      bytes of a static memory buffer. If the buffer
//      becomes too full, simply fail. To free memory,
//      simply reset the pointer (vOffsetMemBlock)
//      back to zero. This memory scheme is very fast
//      and is optimized for the assumption that the
//      only thing you are using temporary memory
//      for is to hold arguments while you call Excel12f().
//      We rely on the fact that you will free all the
//      temporary memory at the same time. We also
//      assume you will not need more memory than
//      the amount required to hold a few arguments
//      to Excel12f().
// 
//      Note that the memory allocation algorithm
//      now supports multithreaded applications by
//      giving each unique thread its own static
//      block of memory. This is handled in the
//      MemoryManager.cpp file automatically.
//
///***************************************************************************

LPSTR GetTempMemory(size_t cBytes)
{
	return MGetTempMemory(cBytes);
}

///***************************************************************************
// FreeAllTempMemory()
//
// Purpose: 
//
//      Frees all temporary memory that has been allocated
//      for the current thread
//
// Parameters:
//
// Returns: 
//
// Comments:
//
//
///***************************************************************************

void FreeAllTempMemory(void)
{
	MFreeAllTempMemory();
}

///***************************************************************************
// Excel()
//
// Purpose: 
//      A fancy wrapper for the Excel4() function. It also
//      does the following:
//        
//      (1) Checks that none of the LPXLOPER arguments are 0,
//          which would indicate that creating a temporary XLOPER
//          has failed. In this case, it doesn't call Excel
//          but it does print a debug message.
//      (2) If an error occurs while calling Excel,
//          print a useful debug message.
//      (3) When done, free all temporary memory.
//        
//      #1 and #2 require _DEBUG to be defined.
//
// Parameters:
//
//      int xlfn            Function number (xl...) to call
//      LPXLOPER pxResult   Pointer to a place to stuff the result,
//                          or 0 if you don't care about the result.
//      int count           Number of arguments
//      ...                 (all LPXLOPERs) - the arguments.
//
// Returns: 
//
//      A return code (Some of the xlret... values, as defined
//      in XLCALL.H, OR'ed together).
//
// Comments:
//
///***************************************************************************

//int cdecl Excel(int xlfn, LPXLOPER pxResult, int count, ...)
//{
//	int xlret;
//	va_list ppxArgs;
//
//#ifdef _DEBUG
//	LPXLOPER px;
//	int i;
//
//	va_start(ppxArgs, count);
//
//	for (i = 0; i < count; i++)
//	{
//		px = va_arg(ppxArgs, LPXLOPER);
//
//		if (px == NULL)
//		{
//			debugPrintf("Error! A zero parameter is being passed "
//				"to Excel. Check memory...\r");
//			FreeAllTempMemory();
//			return xlretFailed;
//		}
//	}
//
//	va_end(ppxArgs);
//
//#endif /* DEBUG */
//
//	va_start(ppxArgs, count);
//	xlret = Excel4v(xlfn, pxResult, count, (LPXLOPER*)ppxArgs);
//	va_end(ppxArgs);
//
//#ifdef _DEBUG
//
//	if (xlret != xlretSuccess)
//	{
//		debugPrintf("Error! Excel4(");
//
//		if (xlfn & xlCommand)
//			debugPrintf("xlCommand | ");
//		if (xlfn & xlSpecial)
//			debugPrintf("xlSpecial | ");
//		if (xlfn & xlIntl)
//			debugPrintf("xlIntl | ");
//		if (xlfn & xlPrompt)
//			debugPrintf("xlPrompt | ");
//
//		debugPrintf("%u) callback failed:", xlfn & 0x0FFF);
//
//		/* More than one error bit may be on */
//
//		if (xlret & xlretAbort)
//		{
//			debugPrintf(" Macro Halted\r");
//		}
//
//		if (xlret & xlretInvXlfn)
//		{
//			debugPrintf(" Invalid Function Number\r");
//		}
//
//		if (xlret & xlretInvCount)
//		{
//			debugPrintf(" Invalid Number of Arguments\r");
//		}
//
//		if (xlret & xlretInvXloper)
//		{
//			debugPrintf(" Invalid XLOPER\r");
//		}
//
//		if (xlret & xlretStackOvfl)
//		{
//			debugPrintf(" Stack Overflow\r");
//		}
//
//		if (xlret & xlretFailed)
//		{
//			debugPrintf(" Command failed\r");
//		}
//
//		if (xlret & xlretUncalced)
//		{
//			debugPrintf(" Uncalced cell\r");
//		}
//
//	}
//
//#endif /* DEBUG */
//
//	FreeAllTempMemory();
//
//	return xlret;
//}

///***************************************************************************
// Excel12f()
//
// Purpose: 
//      A fancy wrapper for the Excel12() function. It also
//      does the following:
//        
//      (1) Checks that none of the LPXLOPER12 arguments are 0,
//          which would indicate that creating a temporary XLOPER12
//          has failed. In this case, it doesn't call Excel12
//          but it does print a debug message.
//      (2) If an error occurs while calling Excel12,
//          print a useful debug message.
//      (3) When done, free all temporary memory.
//        
//      #1 and #2 require _DEBUG to be defined.
//
// Parameters:
//
//      int xlfn            Function number (xl...) to call
//      LPXLOPER12 pxResult Pointer to a place to stuff the result,
//                          or 0 if you don't care about the result.
//      int count           Number of arguments
//      ...                 (all LPXLOPER12s) - the arguments.
//
// Returns: 
//
//      A return code (Some of the xlret... values, as defined
//      in XLCALL.H, OR'ed together).
//
// Comments:
//
///***************************************************************************

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



///***************************************************************************
// TempNum()
//
// Purpose: 
//      Creates a temporary numeric (IEEE floating point) XLOPER.
//
// Parameters:
//
//      double d        The value
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

//LPXLOPER TempNum(double d)
//{
//	LPXLOPER lpx;
//
//	lpx = (LPXLOPER)GetTempMemory(sizeof(XLOPER));
//
//	if (!lpx)
//	{
//		return 0;
//	}
//
//	lpx->xltype = xltypeNum;
//	lpx->val.num = d;
//
//	return lpx;
//}

///***************************************************************************
// TempNum12()
//
// Purpose: 
//      Creates a temporary numeric (IEEE floating point) XLOPER12.
//
// Parameters:
//
//      double d        The value
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//		
//
///***************************************************************************

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

///***************************************************************************
// TempStr()
//
// Purpose: 
//      Creates a temporary string XLOPER
//
// Parameters:
//
//      LPSTR lpstr     The string, as a null-terminated
//                      C string, with the first byte
//                      undefined. This function will
//                      count the bytes of the string
//                      and insert that count in the
//                      first byte of lpstr. Excel cannot
//                      handle strings longer than 255
//                      characters.
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
//      (1) This function has the side effect of inserting
//          the byte count as the first character of
//          the created string.
//
//      (2) For highest speed, with constant strings,
//          you may want to manually count the length of
//          the string before compiling, and then avoid
//          using this function.
//
//      (3) Behavior is undefined for non-null terminated
//          input or strings longer than 255 characters.
//
// Note: If lpstr passed into TempStr is readonly, TempStr
// will crash your XLL as it will try to modify a read only
// string in place. strings declared on the stack as described below
// are read only by default in VC++
//
// char *str = " I am a string"
//
// Use extreme caution while calling TempStr on such strings. Refer to
// VC++ documentation for complier options to ensure that these strings
// are compiled as read write or use TempStrConst instead.
//
// TempStr is provided mainly for backwards compatability and use of
// TempStrConst is encouraged going forward.
//
///***************************************************************************

//LPXLOPER TempStr(LPSTR lpstr)
//{
//	LPXLOPER lpx;
//
//	lpx = (LPXLOPER)GetTempMemory(sizeof(XLOPER));
//
//	if (!lpx)
//	{
//		return 0;
//	}
//
//	lpstr[0] = (BYTE)strlen(lpstr + 1);
//	lpx->xltype = xltypeStr;
//	lpx->val.str = lpstr;
//
//	return lpx;
//}


///***************************************************************************
// TempStrConst()
//
// Purpose: 
//      Creates a temporary string XLOPER from a
//      const string with a local copy in temp memory
//
// Parameters:
//
//      LPSTR lpstr     The string, as a null-terminated
//                      C string. This function will
//                      count the bytes of the string
//                      and insert that count in the
//                      first byte of the temp string.
//                      Excel cannot handle strings
//                      longer than 255 characters.
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
//      Will take a string of the form "abc\0" and make a
//      temp XLOPER of the form "\003abc"
//
//
///***************************************************************************

LPXLOPER TempStrConst(const LPSTR lpstr)
{
	LPXLOPER lpx;
	LPSTR lps;
	size_t len;

	len = strlen(lpstr);

	lpx = (LPXLOPER)GetTempMemory(sizeof(XLOPER) + len + 1);

	if (!lpx)
	{
		return 0;
	}

	lps = (LPSTR)lpx + sizeof(XLOPER);

	lps[0] = (BYTE)len;
	//can't strcpy_s because of removal of null-termination
	memcpy_s(lps + 1, len + 1, lpstr, len);
	lpx->xltype = xltypeStr;
	lpx->val.str = lps;

	return lpx;
}

///***************************************************************************
// TempStr12()
//
// Purpose: 
//      Creates a temporary string XLOPER12 from a
//      unicode const string with a local copy in
//      temp memory
//
// Parameters:
//
//      XCHAR lpstr     The string, as a null-terminated
//                      unicode string. This function will
//                      count the bytes of the string
//                      and insert that count in the
//                      first byte of the temp string.
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
//      (1) Fix for const string pointers being passed in to TempStr.
//          Note it assumes NO leading space
//
//      (2) Also note that XLOPER12 now uses unicode for the string
//          operators
//
//      (3) Will remove the null-termination on the string
//
//
//
// Note: TempStr12 is different from TempStr and is more like TempStrConst
// in its behavior. We have consciously made this choice and deprecated the
// behavior of TempStr going forward. Refer to the note in comment section
// for TempStr to better understand this design decision.
///***************************************************************************

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

///***************************************************************************
// TempBool()
//
// Purpose: 
//      Creates a temporary logical (true/false) XLOPER.
//
// Parameters:
//
//      int b           0 - for a FALSE XLOPER
//                      Anything else - for a TRUE XLOPER
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

//LPXLOPER TempBool(int b)
//{
//	LPXLOPER lpx;
//
//	lpx = (LPXLOPER)GetTempMemory(sizeof(XLOPER));
//
//	if (!lpx)
//	{
//		return 0;
//	}
//
//	lpx->xltype = xltypeBool;
//#ifdef __cplusplus
//	lpx->val.xbool = b ? 1 : 0;
//#else
//	lpx->val.bool = b ? 1 : 0;
//#endif
//
//	return lpx;
//}

///***************************************************************************
// TempBool12()
//
// Purpose: 
//      Creates a temporary logical (true/false) XLOPER12.
//
// Parameters:
//
//      BOOL b          0 - for a FALSE XLOPER12
//                      Anything else - for a TRUE XLOPER12
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

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

///***************************************************************************
// TempInt()
//
// Purpose: 
//      Creates a temporary integer XLOPER.
//
// Parameters:
//
//      short int i     The integer
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

//LPXLOPER TempInt(short int i)
//{
//	LPXLOPER lpx;
//
//	lpx = (LPXLOPER)GetTempMemory(sizeof(XLOPER));
//
//	if (!lpx)
//	{
//		return 0;
//	}
//
//	lpx->xltype = xltypeInt;
//	lpx->val.w = i;
//
//	return lpx;
//}

///***************************************************************************
// TempInt12()
//
// Purpose: 
//          Creates a temporary integer XLOPER12.
//
// Parameters:
//
//      int i           The integer
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
//      Note that the int oper has increased in size from
//      short int up to int in the 12 opers
//
///***************************************************************************

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

///***************************************************************************
// TempErr()
//
// Purpose: 
//      Creates a temporary error XLOPER.
//
// Parameters:
//
//      WORD err        The error code. One of the xlerr...
//                      constants, as defined in XLCALL.H.
//                      See the Excel user manual for
//                      descriptions about the interpretation
//                      of various error codes.
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

//LPXLOPER TempErr(WORD err)
//{
//	LPXLOPER lpx;
//
//	lpx = (LPXLOPER)GetTempMemory(sizeof(XLOPER));
//
//	if (!lpx)
//	{
//		return 0;
//	}
//
//	lpx->xltype = xltypeErr;
//	lpx->val.err = err;
//
//	return lpx;
//}

///***************************************************************************
// TempErr12()
//
// Purpose: 
//      Creates a temporary error XLOPER12.
//
// Parameters:
//
//      int err         The error code. One of the xlerr...
//                      constants, as defined in XLCALL.H.
//                      See the Excel user manual for
//                      descriptions about the interpretation
//                      of various error codes.
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
//      Note the paramater has changed from a WORD to an int
//      in the new 12 operators
//
///***************************************************************************

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

///***************************************************************************
// TempActiveRef()
//
// Purpose: 
//      Creates a temporary rectangular reference to the active
//      sheet. Remember that the active sheet is the sheet that
//      the user sees in front, not the sheet that is currently
//      being calculated.
//
// Parameters:
//
//      WORD rwFirst    (0 based) The first row in the rectangle.
//      WORD rwLast     (0 based) The last row in the rectangle.
//      BYTE colFirst   (0 based) The first column in the rectangle.
//      BYTE colLast    (0 based) The last column in the rectangle.
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

//LPXLOPER TempActiveRef(WORD rwFirst, WORD rwLast, BYTE colFirst, BYTE colLast)
//{
//	LPXLOPER lpx;
//	LPXLMREF lpmref;
//	int wRet;
//
//	lpx = (LPXLOPER)GetTempMemory(sizeof(XLOPER));
//	lpmref = (LPXLMREF)GetTempMemory(sizeof(XLMREF));
//
//	if (!lpmref)
//	{
//		return 0;
//	}
//
//
//	/* calling Excel() instead of Excel4() would free all temp memory! */
//	wRet = Excel4(xlSheetId, lpx, 0);
//
//	if (wRet != xlretSuccess)
//	{
//		return 0;
//	}
//	else
//	{
//		lpx->xltype = xltypeRef;
//		lpx->val.mref.lpmref = lpmref;
//		lpmref->count = 1;
//		lpmref->reftbl[0].rwFirst = rwFirst;
//		lpmref->reftbl[0].rwLast = rwLast;
//		lpmref->reftbl[0].colFirst = colFirst;
//		lpmref->reftbl[0].colLast = colLast;
//
//		return lpx;
//	}
//}


///***************************************************************************
// TempActiveRef12()
//
// Purpose: 
//      Creates a temporary rectangular reference to the active
//      sheet. Remember that the active sheet is the sheet that
//      the user sees in front, not the sheet that is currently
//      being calculated.
//
// Parameters:
//
//      RW rwFirst      (0 based) The first row in the rectangle.
//      RW rwLast       (0 based) The last row in the rectangle.
//      COL colFirst    (0 based) The first column in the rectangle.
//      COL colLast     (0 based) The last column in the rectangle.
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
//      Note that the formal parameters have changed for Excel 2007
//      The valid size has increased to accomodate the increase 
//      in Excel 2007 workbook sizes
//
///***************************************************************************

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

///***************************************************************************
// TempActiveCell()
//
// Purpose: 
//      Creates a temporary reference to a single cell on the active
//      sheet. Remember that the active sheet is the sheet that
//      the user sees in front, not the sheet that is currently
//      being calculated.
//
// Parameters:
//
//      WORD rw         (0 based) The row of the cell.
//      BYTE col        (0 based) The column of the cell.
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

//LPXLOPER TempActiveCell(WORD rw, BYTE col)
//{
//	return TempActiveRef(rw, rw, col, col);
//}

///***************************************************************************
// TempActiveCell12()
//
// Purpose: 
//      Creates a temporary reference to a single cell on the active
//      sheet. Remember that the active sheet is the sheet that
//      the user sees in front, not the sheet that is currently
//      being calculated.
//
// Parameters:
//
//      RW rw           (0 based) The row of the cell.
//      COL col         (0 based) The column of the cell.
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
//      Paramter types changed to RW and COL to accomodate the increase
//      in sheet sizes introduced in Excel 2007
//
///***************************************************************************

LPXLOPER12 TempActiveCell12(RW rw, COL col)
{
	return TempActiveRef12(rw, rw, col, col);
}

///***************************************************************************
// TempActiveRow()
//
// Purpose: 
//      Creates a temporary reference to an entire row on the active
//      sheet. Remember that the active sheet is the sheet that
//      the user sees in front, not the sheet that is currently
//      being calculated.
//
// Parameters:
//
//      RW rw           (0 based) The row.
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

//LPXLOPER TempActiveRow(WORD rw)
//{
//	return TempActiveRef(rw, rw, 0, 0xFF);
//}

///***************************************************************************
// TempActiveRow12()
//
// Purpose: 
//      Creates a temporary reference to an entire row on the active
//      sheet. Remember that the active sheet is the sheet that
//      the user sees in front, not the sheet that is currently
//      being calculated.
//
// Parameters:
//
//      RW rw           (0 based) The row.
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
//      Paramter type change to RW to accomodate the increase in sheet
//      sizes introduced in Excel 2007
//
///***************************************************************************

LPXLOPER12 TempActiveRow12(RW rw)
{
	return TempActiveRef12(rw, rw, 0, 0x00003FFF);
}

///***************************************************************************
// TempActiveColumn()
//
// Purpose: 
//      Creates a temporary reference to an entire column on the active
//      sheet. Remember that the active sheet is the sheet that
//      the user sees in front, not the sheet that is currently
//      being calculated.
//
// Parameters:
//
//      LPSTR s         First string
//      LPSTR t         Second string
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

//LPXLOPER TempActiveColumn(BYTE col)
//{
//	return TempActiveRef(0, 0xFFFF, col, col);
//}

///***************************************************************************
// TempActiveColumn12()
//
// Purpose: 
//      Creates a temporary reference to an entire column on the active
//      sheet. Remember that the active sheet is the sheet that
//      the user sees in front, not the sheet that is currently
//      being calculated.
//
// Parameters:
//
//      COL col	        (0 based) The column.
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
//      Paramter type change to COL to accomodate the increase in sheet
//      sizes introduced in Excel 2007
//
///***************************************************************************

LPXLOPER12 TempActiveColumn12(COL col)
{
	return TempActiveRef12(0, 0x000FFFFF, col, col);
}


///***************************************************************************
// TempMissing()
//
// Purpose: 
//      This is used to simulate a missing argument when
//      calling Excel(). It creates a temporary
//      "missing" XLOPER.
//
// Parameters:
//
// Returns: 
//
//      LPXLOPER        The temporary XLOPER, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

//LPXLOPER TempMissing(void)
//{
//	LPXLOPER lpx;
//
//	lpx = (LPXLOPER)GetTempMemory(sizeof(XLOPER));
//
//	if (!lpx)
//	{
//		return 0;
//	}
//
//	lpx->xltype = xltypeMissing;
//
//	return lpx;
//}

///***************************************************************************
// TempMissing12()
//
// Purpose: 
//      This is used to simulate a missing argument when
//      calling Excel12f(). It creates a temporary
//      "missing" XLOPER12.
//
// Parameters:
//
// Returns: 
//
//      LPXLOPER12      The temporary XLOPER12, or 0
//                      if GetTempMemory() failed.
//
// Comments:
//
///***************************************************************************

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

///***************************************************************************
// InitFramework()
//
// Purpose: 
//      Initializes all the framework functions.
//
// Parameters:
//
// Returns: 
//
// Comments:
//
///***************************************************************************

void InitFramework(void)
{
	FreeAllTempMemory();
}

///***************************************************************************
// QuitFramework()
//
// Purpose: 
//      Cleans up for all framework functions.
//
// Parameters:
//
// Returns: 
//
// Comments:
//
///***************************************************************************
void QuitFramework(void)
{
	FreeAllTempMemory();
}

///***************************************************************************
// FreeXLOperT()
//
// Purpose: 
//      Will free any malloc'd memory associated with the given
//      LPXLOPER, assuming it has any memory associated with it
//
// Parameters:
//
//      LPXLOPER pxloper    Pointer to the XLOPER whose associated
//                          memory we want to free
//
// Returns: 
//
// Comments:
//
///***************************************************************************

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


///***************************************************************************
// FreeXLOper12T()
//
// Purpose: 
//      Will free any malloc'd memory associated with the given
//      LPXLOPER12, assuming it has any memory associated with it
//
// Parameters:
//
//      LPXLOPER12 pxloper12    Pointer to the XLOPER12 whose
//         	                    associated memory we want to free
//
// Returns: 
//
// Comments:
//
///***************************************************************************


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


///***************************************************************************
// ConvertXLRefToXLRef12()
//
// Purpose: 
//      Will attempt to convert an XLREF into the given XREF12
//
// Parameters:
//
//      LPXLREF pxref       Pointer to the XLREF to copy
//      LPXLREF12 pxref12   Pointer to the XLREF12 to copy into
//
// Returns: 
//
//      BOOL                TRUE if the conversion succeeded, FALSE otherwise
//
// Comments:
//
///***************************************************************************


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

///***************************************************************************
// ConvertXLRef12ToXLRef()
//
// Purpose: 
//      Will attempt to convert an XLREF12 into the given XLREF
//
// Parameters:
//
//      LPXLREF12 pxref12   Pointer to the XLREF12 to copy
//      LPXLREF pxref       Pointer to the XLREF to copy into
//
// Returns: 
//
//      BOOL                TRUE if the conversion succeeded, FALSE otherwise
//
// Comments:
//
///***************************************************************************

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

///***************************************************************************
// XLOper12ToXLOper()
//
// Purpose: 
//      Conversion routine used to convert from the new XLOPER12
//      to the old XLOPER.
//
// Parameters:
//
//      LPXLOPER12 pxloper12    Pointer to the XLOPER12 to copy
//      LPXLOPER pxloper        Pointer to the XLOPER to copy into
//
// Returns: 
//
//      BOOL                    TRUE if the conversion succeeded, FALSE otherwise
//
// Comments:
//
//      - The caller is responsible for freeing any memory associated with
//        the copy if the conversion is a success; FreeXOperT can be
//        used, or it may be done by hand.
//
//      - If the conversion fails, any memory the method needed to malloc
//        up until the point of failure in the conversion will be freed
//        by the method itself during cleanup.
//
///***************************************************************************

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


///***************************************************************************
// XLOperToXLOper12()
//
// Purpose: 
//      Conversion routine used to convert from the old XLOPER
//      to the new XLOPER12.
//
// Parameters:
//
//      LPXLOPER pxloper        Pointer to the XLOPER to copy
//      LPXLOPER12 pxloper12    Pointer to the XLOPER12 to copy into
//
// Returns: 
//
//      BOOL                    TRUE if the conversion succeeded, FALSE otherwise
//
// Comments:
//
//      - The caller is responsible for freeing any memory associated with
//        the copy if the conversion is a success; FreeXLOper12T can be
//        used, or it may be done by hand.
//
//      - If the conversion fails, any memory the method needed to malloc
//        up until the point of failure in the conversion will be freed
//        by the method itself during cleanup.
//
///***************************************************************************


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
