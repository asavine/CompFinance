///***************************************************************************
// File:		MemoryPool.h
//
// Purpose:		Class definition for the memory pool class used by the
//				memory manager.  Each pool is a block of memory set
//				aside for a specific thread for use in creating temporary
//				XLOPER/XLOPER12's in the framework
// 
// Platform:    Microsoft Windows
//
///***************************************************************************

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

//
// Total amount of memory to allocate for all temporary XLOPERs
//

#define MEMORYSIZE  4194304

class MemoryPool
{
public:
	MemoryPool(void);
	~MemoryPool(void);
	void ClearPool(void);
	LPSTR GetTempMemory(int cBytes);
	void FreeAllTempMemory();

	DWORD m_dwOwner;			// ID of ownning thread
	char* m_rgchMemBlock;		// Memory for temporary XLOPERs
	int m_ichOffsetMemBlock;	// Offset of next memory block to allocate
};
