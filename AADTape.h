#pragma once

#include "blocklist.h"
#include "AADNode.h"

constexpr size_t BLOCKSIZE  = 16384;		//	Number of nodes
constexpr size_t ADJSIZE    = 32768;		//	Number of adjoints
constexpr size_t DATASIZE   = 65536;		//	Data in bytes

class Tape
{
	//	Working with multiple results / adjoints?
	static bool							multi;

	//  Storage for adjoints
    blocklist<double, ADJSIZE>			myAdjointsMulti;
    
	//  Storage for derivatives and child adjoint pointers
	blocklist<double, DATASIZE>			myDers;
	blocklist<double*, DATASIZE>		myArgPtrs;

    //  Storage for the nodes
	blocklist<Node, BLOCKSIZE>		    myNodes;

	//	Padding so tapes in a vector don't interfere
    char                                myPad[64];

    friend auto setNumResultsForAAD(const bool, const size_t);
    friend struct numResultsResetterForAAD;
	friend class Number;

public:

    //  Build note in place and return a pointer
	//	N : number of childs (arguments)
    template <size_t N>
    Node* recordNode()
    {
        //  Construct the node in place on tape
        Node* node = myNodes.emplace_back(N);
        
        //  Store and zero the adjoint(s)
        if (multi)
        {
            node->pAdjoints = myAdjointsMulti.emplace_back_multi(Node::numAdj);
            fill(node->pAdjoints, node->pAdjoints + Node::numAdj, 0.0);
        }

		//	Store the derivatives and child adjoint pointers unless leaf
		if constexpr(N)
		{
			node->pDerivatives = myDers.emplace_back_multi<N>();
			node->pAdjPtrs = myArgPtrs.emplace_back_multi<N>();

		}

        return node;
    }

	void resetAdjoints()
	{
		if (multi)
		{
			myAdjointsMulti.memset(0);
		}
		else
		{
			for (Node& node : myNodes)
			{
				node.mAdjoint = 0;
			}
		}
	}

    //  Clear
    void clear()
    {
        myAdjointsMulti.clear();
		myDers.clear();
		myArgPtrs.clear();
        myNodes.clear();
    }

    //  Rewind
    void rewind()
    {

#ifdef _DEBUG

		clear();

#else
		if (multi)
		{
			myAdjointsMulti.rewind();
		}
		myDers.rewind();
		myArgPtrs.rewind();
		myNodes.rewind();

#endif

    }

    //  Set mark
    void mark()
    {
        if (multi)
        {
            myAdjointsMulti.setmark();
        }
		myDers.setmark();
		myArgPtrs.setmark();
		myNodes.setmark();
    }

    //  Rewind to mark
    void rewindToMark()
    {
        if (multi)
        {
            myAdjointsMulti.rewind_to_mark();
        }
		myDers.rewind_to_mark();
		myArgPtrs.rewind_to_mark();
		myNodes.rewind_to_mark();
    }

    //  Iterators
    
    using iterator = blocklist<Node, BLOCKSIZE>::iterator;

    auto begin()
    {
        return myNodes.begin();
    }

    auto end()
    {
        return myNodes.end();
    }

    auto markIt()
    {
        return myNodes.mark();
    }

    auto find(Node* node)
    {
        return myNodes.find(node);
    }
};