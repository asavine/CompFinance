///***************************************************************************
// File:        MemoryManager.cpp
//
// Purpose:     The memory manager class is an update to the memory manager
//              in the previous release of the framework.  This class provides
//              each thread with an array of bytes to use as temporary memory.
//              The size of the array, and the methods for dealing with the
//              memory explicitly, is in the class MemoryPool.  
//
//              MemoryManager handles assigning of threads to pools, and the
//              creation of new pools when a thread asks for memory the first
//              time.  Using a singleton class, the manager provides an interface
//              to C code into the manager.  The number of unique pools starts
//              as MEMORYPOOLS, defined in MemoryManager.h.  When a new thread
//              needs a pool, and the current set of pools are all assigned,
//              the number of pools increases by a factor of two.
// 
// Platform:    Microsoft Windows
//
///***************************************************************************

#include "MemoryManager.h"

//
// Singleton instance of the class
//
MemoryManager* vpmm;

//
// Interface for C callers to ask for memory
//
// See MemoryPool.h for more details
//
LPSTR MGetTempMemory(int cByte)
{
	return MemoryManager::GetManager()->CPP_GetTempMemory(cByte);
}

//
// Interface for C callers to allow their memory to be reused
//
// See MemoryPool.h for more details
//
void MFreeAllTempMemory()
{
	MemoryManager::GetManager()->CPP_FreeAllTempMemory();
}

//
// Returns the singleton class, or creates one if it doesn't exit
//
MemoryManager* MemoryManager::GetManager()
{
	if (!vpmm)
	{
		vpmm = new MemoryManager();
	}
	return vpmm;
}

//
// Default constructor
//
MemoryManager::MemoryManager(void)
{
	m_impCur = 0;
	m_impMax = MEMORYPOOLS;
	m_rgmp = new MemoryPool[MEMORYPOOLS];
}

//
// Destructor.  Because of the way memory pools get copied,
// this function needs to call an additional function to clear
// up the MemoryPool memory - the deconstructor on MemoryPool
// does not actually delete its memory
//
MemoryManager::~MemoryManager(void)
{
	MemoryPool* pmp = m_rgmp;
	int i;

	for (i = 0; i < m_impMax; i++)
	{
		if (pmp->m_rgchMemBlock)
		{
			pmp->ClearPool();
		}
		pmp++;
	}
	delete [] m_rgmp;
}

//
// Method that will query the correct memory pool of the calling
// thread for a set number of bytes.  Returns 0 if there was a
// failure in getting the memory.
//
LPSTR MemoryManager::CPP_GetTempMemory(int cByte)
{
	DWORD dwThreadID;
	MemoryPool* pmp;

	dwThreadID = GetCurrentThreadId(); //the id of the calling thread
	pmp = GetMemoryPool(dwThreadID);

	if (!pmp) //no more room for pools
	{
		return 0;
	}

	return pmp->GetTempMemory(cByte);
}

//
// Method that tells the pool owned by the calling thread that
// it is free to reuse all of its memory
//
void MemoryManager::CPP_FreeAllTempMemory()
{
	DWORD dwThreadID;
	MemoryPool* pmp;

	dwThreadID = GetCurrentThreadId(); //the id of the calling thread
	pmp = GetMemoryPool(dwThreadID);

	if (!pmp) //no more room for pools
	{
		return;
	}

	pmp->FreeAllTempMemory();
}

//
// Method iterates through the memory pools in an attempt to find
// the pool that matches the given thread ID. If a pool is not found,
// it creates a new one
//
MemoryPool* MemoryManager::GetMemoryPool(DWORD dwThreadID)
{
	int imp; //loop var
	MemoryPool* pmp; //current pool

	pmp = m_rgmp;

	for (imp = 0; imp < m_impCur; imp++)
	{
		if (pmp->m_dwOwner == dwThreadID)
		{
			return pmp;
		}

		pmp++;
	}

	return CreateNewPool(dwThreadID); //didn't find the owner, make a new one
}

//
// Will assign an unused pool to a thread; should all pools be assigned,
// it will grow the number of pools available.
//
MemoryPool* MemoryManager::CreateNewPool(DWORD dwThreadID)
{
	if (m_impCur >= m_impMax)
	{
		GrowPools();
	}
	m_rgmp[m_impCur++].m_dwOwner = dwThreadID;

	return m_rgmp+m_impCur-1;
}

//
// Increases the number of available pools by a factor of two. All of
// the old pools have their memory pointed to by the new pools. The
// memory for the new pools that get replaced is first freed. The reason
// ~MemoryPool() can't free its array is in this method - they would be
// deleted when the old array of pools is freed at the end of the method,
// despite the fact they are now being pointed to by the new pools.
//
void MemoryManager::GrowPools()
{
	MemoryPool* rgmpTemp;
	MemoryPool* pmpDst;
	MemoryPool* pmpSrc;
	int i;

	pmpDst = rgmpTemp = new MemoryPool[2*m_impMax];
	pmpSrc = m_rgmp;

	for (i = 0; i < m_impCur; i++)
	{
		delete [] pmpDst->m_rgchMemBlock;
		pmpDst->m_rgchMemBlock = pmpSrc->m_rgchMemBlock;
		pmpDst->m_dwOwner = pmpSrc->m_dwOwner;

		pmpDst++;
		pmpSrc++;
	}
	delete [] m_rgmp;
	m_rgmp = rgmpTemp;
}
