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

#ifdef __cplusplus

#include "MemoryPool.h"

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

	LPSTR CPP_GetTempMemory(int cByte);
	void CPP_FreeAllTempMemory();

private:
	MemoryPool* CreateNewPool(DWORD dwThreadID);
	MemoryPool* GetMemoryPool(DWORD dwThreadID);
	void GrowPools();

	int m_impCur;		// Current number of pools
	int m_impMax;		// Max number of mem pools
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
	LPSTR MGetTempMemory(int cByte);
	void MFreeAllTempMemory();
#ifdef __cplusplus
}
#endif 
