#pragma once

#define AADET true

#if AADET

#include "AADExpr.h"

#else

#include "AADNumber.h"

#endif

template <class IT>
inline void putOnTape(IT begin, IT end)
{
    for (auto it = begin; it != end; ++it)
        it->putOnTape();
}

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
