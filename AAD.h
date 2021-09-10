
/*
Written by Antoine Savine in 2018

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications 
is freely granted to any person or company who purchased a copy of the book

Modern Computational Finance: AAD and Parallel Simulations
Antoine Savine
Wiley, 2018

As long as this comment is preserved at the top of the file
*/

#pragma once

#include <memory>

//  So we can instrument Gaussians like standard math functions
#include "gaussians.h"

//  Use traditional AAD of chapter 10 (false)
//      or expression templated (AADET) of chapter 15 (true)
#define AADET   true

#if AADET

#include "AADExpr.h"

#else

#include "AADNumber.h"

#endif

//  Routines for multi-dimensional AAD (chapter 14)
//  Set static context for multi-dimensional AAD

//	RAII: reset dimension 1 on destruction
struct numResultsResetterForAAD
{
	~numResultsResetterForAAD()
	{
		Tape::multi = false;
		Node::numAdj = 1;
	}
};

//  Routine: set dimension and get RAII resetter
inline auto setNumResultsForAAD(const bool multi = false, const size_t numResults = 1)
{
	Tape::multi = multi;
	Node::numAdj = numResults;
	return make_unique<numResultsResetterForAAD>();
}

//  Other utilities

//	Put collection on tape
template <class IT>
inline void putOnTape(IT begin, IT end)
{
    for_each(begin, end, [](Number& n) {n.putOnTape(); });
}

//	Convert collection between double and Number or reverse
template<class It1, class It2>
inline void convertCollection(It1 srcBegin, It1 srcEnd, It2 destBegin)
{
    using destType = remove_reference_t<decltype(*destBegin)>;
    transform(srcBegin, srcEnd, destBegin, 
        [](const auto& source) { return destType(source); });
}
