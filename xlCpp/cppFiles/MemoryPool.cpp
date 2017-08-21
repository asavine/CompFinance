///***************************************************************************
// File:        MemoryPool.cpp
//
// Purpose:     A memory pool is an array of characters that is pre-allocated,
//              and used as temporary memory by the caller. The allocation
//              algorithm is very simple. When a thread asks for some memory,
//              the index into the array moves forward by that many bytes, and
//              a pointer is returned to the previous index before the pointer
//              was advanced. When a call comes to free all of the memory, the
//              pointer is set back to the beginning of the array.
//
//              Each pool has MEMORYSIZE bytes of storage space available to it
// 
// Platform:    Microsoft Windows
//
///***************************************************************************

#include "MemoryPool.h"

//
// Constructor creates the memory to be used by the pool
// and starts the index at the beginning.
//
MemoryPool::MemoryPool(void)
{
	m_rgchMemBlock = new char[MEMORYSIZE];
	m_ichOffsetMemBlock = 0;
	m_dwOwner = (DWORD)-1;
}

//
// An empty destructor - see reasoning below
//
MemoryPool::~MemoryPool(void)
{
}

//
// Unable to delete the memory block when we delete the pool,
// as it may be still be in use due to a GrowPools() call; this
// method will actually delete the pool's memory
//

void MemoryPool::ClearPool(void)
{
	delete [] m_rgchMemBlock;
}

//
// Advances the index forward by the given number of bytes.
// Should there not be enough memory, or the number of bytes
// is not allowed, this method will return 0. Can be called
// and used exactly as malloc().
//
LPSTR MemoryPool::GetTempMemory(int cBytes)
{
	LPSTR lpMemory;

	if (m_ichOffsetMemBlock + cBytes > MEMORYSIZE || cBytes <= 0)
	{
		return 0;
	}
	else
	{
		lpMemory = (LPSTR) m_rgchMemBlock + m_ichOffsetMemBlock;
		m_ichOffsetMemBlock += cBytes;

		return lpMemory;
	}
}

//
// Frees all the temporary memory by setting the index for
// available memory back to the beginning
//
void MemoryPool::FreeAllTempMemory()
{
	m_ichOffsetMemBlock = 0;
}
