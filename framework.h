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
