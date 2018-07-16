
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

//  AAD with expression templates
#define AADET   true        

//  So we can instrument Gaussians like standard math functions
#include "gaussians.h"

#if AADET

#include "AADExpr.h"

#else

#include "AADNumber.h"

#endif

//	Globally set number of results (adjoints) on node and tape
//	Get an object that resets to 1 on destruction

struct numResultsResetterForAAD
{
	~numResultsResetterForAAD()
	{
		Tape::multi = false;
		Node::numAdj = 1;
	}
};

inline auto setNumResultsForAAD(const bool multi = false, const size_t numResults = 1)
{
	Tape::multi = multi;
	Node::numAdj = numResults;
	return make_unique<numResultsResetterForAAD>();
}

//	Put collection on tape

template <class IT>
inline void putOnTape(IT begin, IT end)
{
    for_each(begin, end, [](Number& n) {n.putOnTape(); });
}

//	Converters

template<class To, class From>
struct Convert;

template<>
struct Convert<double, double>
{
    static inline double convert(const double from)
    {
        return from;
    }
};

template<class T>
struct Convert<T, T>
{
    static inline T convert(const T from)
    {
        return from;
    }
};

template<>
struct Convert<Number, double>
{
    static inline Number convert(const double from)
    {
        return Number(from);
    }
};

template<class T>
struct Convert<double, T>
{
    static inline double convert(const T from)
    {
        return from.value();
    }
};

template<class To, class From>
inline To convert(const From from)
{
    return Convert<To, From>::convert(from);
}

template<class It1, class It2>
inline void convertCollection(It1 srcBegin, It1 srcEnd, It2 destBegin)
{
    using destType = remove_reference_t<decltype(*destBegin)>;

    while (srcBegin != srcEnd)
    {
        *destBegin++ = convert<destType>(*srcBegin++);
    }
}
