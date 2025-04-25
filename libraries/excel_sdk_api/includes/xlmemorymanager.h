///***************************************************************************
// File:        MemoryManager.h
//
// Purpose:     Class definition for the memory manager used in the framework
//              library.  Includes exported methods for accessing the class
//              in C.
// 
// Platform:    Microsoft Windows
//
///***************************************************************************

#pragma once

#ifdef __cplusplus

#include "xlMemoryPool.h"

//
// Total number of memory allocation pools to manage
//

#define MEMORYPOOLS 4

class MemoryManager
{
public:
	MemoryManager(void);
	~MemoryManager(void);

	static MemoryManager* GetManager();

	LPSTR CPP_GetTempMemory(size_t cByte);
	void CPP_FreeAllTempMemory();

private:
	MemoryPool* CreateNewPool(DWORD dwThreadID);
	MemoryPool* GetMemoryPool(DWORD dwThreadID);
	void GrowPools();

	size_t m_impCur;		// Current number of pools
	size_t m_impMax;		// Max number of mem pools
	MemoryPool* m_rgmp;	// Storage for the memory pools
};

#endif //__cplusplus

//
// Defines functions for accessing class from a C projects
//

#ifdef __cplusplus
extern "C"
{
#endif
	LPSTR MGetTempMemory(size_t cByte);
	void MFreeAllTempMemory();
#ifdef __cplusplus
}
#endif 
