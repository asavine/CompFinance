#pragma once

#include <cstdlib>
#include <exception>

#include <vector>
#include <list>

#include "AADNode.h"

using namespace std;

#define DEFAULT_BLOCK_SIZE 524288
#define DEFAULT_INDEX_SIZE 32768

using Ptr = char*;

class MemoryBlock 
{
    //	Pointers to the first, 
    //      next available and last memory address
    Ptr myBegin;
    Ptr myNext;
    Ptr	myEnd;

public:

    //	Constructor allocates memory, throws if unsuccessful
    MemoryBlock(const size_t size = DEFAULT_BLOCK_SIZE)
    {
        myBegin = (Ptr)malloc(size);
        if (!myBegin) throw bad_alloc();

        myNext = myBegin;
        myEnd = myBegin + size;
    }

    //	Provides the requested amount of memory if available or nullptr
    Ptr requestMemory(const size_t size)
    {
        if (myNext + size > myEnd) return nullptr;

        Ptr prevNext = myNext;
        myNext += size;
        return prevNext;
    }

    //	Disable copy construction and assignment
    MemoryBlock(const MemoryBlock& rhs) = delete;
    MemoryBlock& operator=(const MemoryBlock& rhs) = delete;

    //	But enable move semantics
    MemoryBlock(MemoryBlock&& rhs)
    {
        myBegin = rhs.myBegin;
        myNext = rhs.myNext;
        myEnd = rhs.myEnd;
        rhs.myBegin = rhs.myNext = rhs.myEnd = nullptr;
    }
    MemoryBlock& operator=(MemoryBlock&& rhs)
    {
        if (this == &rhs) return *this;
        myBegin = rhs.myBegin;
        myNext = rhs.myNext;
        myEnd = rhs.myEnd;
        rhs.myBegin = rhs.myNext = rhs.myEnd = nullptr;
        return *this;
    }

    //	Release the block
    void free()
    {
        if (myBegin) ::free(myBegin);
        myBegin = nullptr;
    }

    //	Destructor frees the allocated memory
    ~MemoryBlock() { free(); }

    //	Rewind the block without freeing the memory
    void rewind()
    {
        myNext = myBegin;
    }

    //  Put a mark on the current position

private:
    
    Ptr myMark;

public:
    
    void mark()
    {
        myMark = myNext;
    }

    //  Rewind to mark
    void rewindToMark()
    {
        myNext = myMark;
    }
};

class Tape
{
    //	Memory blocks for the storage of tape entries
    list<MemoryBlock> myBlocks;
    //  Size of each block
    size_t myBlockSize;		

    //	In order to iterate through the tape in forward or reverse order, 
    //      we need to keep a list of entries as they are constructed
    //	We keep pointers to entries in a list of vectors (for efficiency) 
    //      and provide iterators

    //	Memory blocks for indices
    //	List of vectors of pointers to allocated bits
    list<vector<Node*>> myPointers;
    //	Size of each vector 
    size_t myVectorSize;

    //  Current state = position 
    struct State
    {
        //  block
        list<MemoryBlock>::iterator block;
        //  vector
        list<vector<Node*>>::iterator vec;
        //  position in vector
        //  index of the next free slot in the current vector
        size_t idx;
    };
    State myState;

    //	Private methods

    //	Create a new block
    void newBlock()
    {
        myBlocks.push_back(MemoryBlock(myBlockSize));
        
        //  current = last
        myState.block = myBlocks.end();
        --myState.block;
    }

    //	Move to next block, if none create one
    void nextBlock()
    {
        ++myState.block;
        if (myState.block == myBlocks.end()) newBlock();
    }

    //	Create a new vector of pointers
    void newVector()
    {
        myPointers.push_back(vector<Node*>(myVectorSize));
        
        //  current = last
        myState.vec = myPointers.end();
        --myState.vec;

        //  set index (in vector) to 0
        myState.idx = 0;
    }

    //	Move to next vector, if none create one
    void nextVector()
    {
        ++myState.vec;
        myState.idx = 0;
        if (myState.vec == myPointers.end()) newVector();
    }

public:

    Tape(const size_t blockSize = DEFAULT_BLOCK_SIZE,
        const size_t vectorSize = DEFAULT_INDEX_SIZE)
        :
        myBlockSize(blockSize), myVectorSize(vectorSize)
    {
        //  Create one block on construction
        myBlocks.push_back(MemoryBlock(blockSize));
        myPointers.push_back(vector<Node*>(vectorSize));

        //  Set state
        myState.block = myBlocks.begin();
        myState.vec = myPointers.begin();
        myState.idx = 0;
    }

    //  Get memory, this is what we call in place of new/malloc/calloc
    //  Note the is no free - all tape is released at once
    //	We assume that the requested size is always smaller than the block size, otherwise will fail
    Ptr allocate(const size_t size)
    {
        //  Get memory from current block
        Ptr mem = myState.block->requestMemory(size);

        if (!mem) 
        {
            //	Failed : try next block, create one if necessary, then retry
            nextBlock();
            mem = myState.block->requestMemory(size);
            //	If allocation failed, will have thrown by now, 
            //      the only way it could have failed is if size > myBlockSize, 
            //      which we assume never happens
        }

        //	Register pointer to the node so we can iterate
        (*myState.vec)[myState.idx++] = reinterpret_cast<Node*>(mem);
        if (myState.idx == myVectorSize) nextVector();

        return mem;
    }

    //	Rewind the tape without freeing the memory
    void rewind()
    {
        //	Rewind all blocks
        for (auto& block : myBlocks) block.rewind();

        //	Reset state
        myState.block = myBlocks.begin();
        myState.vec = myPointers.begin();
        myState.idx = 0;
    }

    //  Mark the current position
private:
    State myMark;
public:
    void mark()
    {
        myMark = myState;
        //  Put a mark on the current block
        myState.block->mark();
    }

    //	Rewind to mark
    void rewindToMark()
    {
        //	Rewind all blocks after mark
        
        //  position on the block next to mark
        auto it = myMark.block;
        ++it;
        //  rewind
        while (it != myBlocks.end())
        {
            it->rewind();
            ++it;
        }

        //  Rewind the marked block
        myMark.block->rewindToMark();

        //	Reset state
        myState = myMark;
    }

    //	Free all the allocated memory
    void clear() 
    {

        myBlocks.clear();
        myPointers.clear();

        //  Create one block so the tape is usable
        myBlocks.push_back(MemoryBlock(myBlockSize));
        myPointers.push_back(vector<Node*>(myVectorSize));

        myState.block = myBlocks.begin();
        myState.vec = myPointers.begin();
        myState.idx = 0;
    }

    //  Tape iterators

    class iterator
    {
        //  Vector and index
        list<vector<Node*>>::iterator myVector;
        size_t myIdx;

    public:

        //	Default constructor 
        iterator() {}

        //	Constructor 
        iterator(const list<vector<Node*>>::iterator vec, 
            const size_t idx)
            : myVector(vec), myIdx(idx) {}

        //	Pre-increment (we do not provide post)
        iterator& operator++()
        {
            ++myIdx;
            if (myIdx >= myVector->size()) 
            {
                ++myVector;
                myIdx = 0;
            }

            return *this;
        }

        //	Pre-decrement 
        iterator& operator--()
        {
            if (myIdx > 0) 
            {
                --myIdx;
            }
            else 
            {
                --myVector;
                myIdx = myVector->size() - 1;
            }

            return *this;
        }

        //	Access to the node by dereferencing
        Node& operator*()
        {
            return *(*myVector)[myIdx];
        }
        const Node& operator*() const
        {
            return *(*myVector)[myIdx];
        }
        Node* operator->()
        {
            return (*myVector)[myIdx];
        }
        const Node* operator->() const
        {
            return (*myVector)[myIdx];
        }

        //	Check equality
        bool operator ==(const iterator& rhs) const
        {
            return (myIdx == rhs.myIdx && myVector == rhs.myVector);
        }
        bool operator !=(const iterator& rhs) const
        {
            return (myIdx != rhs.myIdx || myVector != rhs.myVector);
        }
    };

    //  Basic iterator facilities
    iterator begin()
    {
        return iterator(myPointers.begin(), 0);
    }

    iterator end()
    {
        return iterator(myState.vec, myState.idx);
    }

    //  Iterator on last
    iterator back()
    {
        iterator it = end();
        --it;
        return it;
    }

    //  Iterator on mark
    iterator markIt()
    {
        return iterator(myMark.vec, myMark.idx);
    }

    //  Find specific node, searching from the end
    iterator find(Node* node)
    {
        //	Search from the end
        iterator it = end();
        iterator b = begin();

        while (it != b) 
        {
            --it;
            if (&*it == node) return it;
        }

        if (&*it == node) return it;

        return end();
    }
};