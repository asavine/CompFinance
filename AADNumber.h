
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

//  Traditional AAD implementation of chapter 10
//  (With multi-dimensional additions of chapter 14)

//  The custom number type

#include <algorithm>
#include "AADTape.h"

class Number
{
    //  Value and node are the only data members
    double myValue;
    Node* myNode;

    //  Create node on tape
	template <size_t N>
	void createNode()
    {
		myNode = tape->recordNode<N>();
    }

	//	Access node (friends only)
    //  Note const incorectness
	Node& node() const 
    {

#ifdef _DEBUG

        //  Help identify errors when arguments are not on tape

		//	Find node on tape
		auto it = tape->find(myNode);

		//	Not found
		if (it == tape->end())
		{
			throw runtime_error("Put a breakpoint here");
		}

#endif
        //  Const incorrectness
		return const_cast<Node&>(*myNode);
    }

    //	Convenient access to node data for friends

    double& derivative() { return myNode->pDerivatives[0]; }
    double& lDer() { return myNode->pDerivatives[0]; }
    double& rDer() { return myNode->pDerivatives[1]; }

    double*& adjPtr() { return myNode->pAdjPtrs[0]; }
    double*& leftAdj() { return myNode->pAdjPtrs[0]; }
    double*& rightAdj() { return myNode->pAdjPtrs[1]; }

	//	Private constructors for operator overloading
	
	//	Unary
	Number(Node& arg, const double val) :
		myValue(val)
    {
        createNode<1>();

		myNode->pAdjPtrs[0] = Tape::multi
			? arg.pAdjoints 
			: &arg.mAdjoint;
    }

	//	Binary
	Number(Node& lhs, Node& rhs, const double val) :
		myValue(val)
	{
        createNode<2>();
        
        if (Tape::multi)
		{
			myNode->pAdjPtrs[0] = lhs.pAdjoints;
			myNode->pAdjPtrs[1] = rhs.pAdjoints;
		}
		else
    {
			myNode->pAdjPtrs[0] = &lhs.mAdjoint;
			myNode->pAdjPtrs[1] = &rhs.mAdjoint;
		}
    }

public:

    //  Static access to tape
    static thread_local Tape* tape;

    //  Public constructors for leaves

    Number() {}

    //  Put on tape on construction
    explicit Number(const double val) :
		myValue(val)
    {
        createNode<0>();
    }

    //  Put on tape on assignment
    Number& operator=(const double val)
    {
        myValue = val;
		createNode<0>();

        return *this;
    }

    //  Explicitly put existing Number on tape
    void putOnTape()
    {
		createNode<0>();
    }

    //  Explicit coversion to double
    explicit operator double& () { return myValue; }
    explicit operator double() const { return myValue; }

    //  Accessors: value and adjoint

    double& value()
    {
        return myValue;
    }
    double value() const
    {
        return myValue;
    }
    //  Single dimensional
    double& adjoint()
    {
        return myNode->adjoint();
    }
    double adjoint() const
    {
        return myNode->adjoint();
    }
    //  Multi dimensional
	double& adjoint(const size_t n)
	{
		return myNode->adjoint(n);
	}
	double adjoint(const size_t n) const
	{
		return myNode->adjoint(n);
    }

    //  Reset all adjoints on the tape
	//		note we don't use this method
    void resetAdjoints()
    {
		tape->resetAdjoints();
    }

	//  Propagation

    //  Propagate adjoints
    //      from and to both INCLUSIVE
    static void propagateAdjoints(
        Tape::iterator propagateFrom,
        Tape::iterator propagateTo)
    {
        auto it = propagateFrom;
        while (it != propagateTo)
        {
            it->propagateOne();
            --it;
        }
        it->propagateOne();
    }

    //  Convenient overloads

    //  Set the adjoint on this node to 1,
    //  Then propagate from the node
    void propagateAdjoints(
        //  We start on this number's node
        Tape::iterator propagateTo)
    {
        //  Set this adjoint to 1
        adjoint() = 1.0;
        //  Find node on tape
        auto propagateFrom = tape->find(myNode);
        propagateAdjoints(propagateFrom, propagateTo);
    }

    //  These 2 set the adjoint to 1 on this node
    void propagateToStart()
    {
        propagateAdjoints(tape->begin());
    }
    void propagateToMark()
    {
        propagateAdjoints(tape->markIt());
    }

    //  This one only propagates
    //  Note: propagation starts at mark - 1
    static void propagateMarkToStart()
    {
        propagateAdjoints(prev(tape->markIt()), tape->begin());
    }

    //  Multi dimensional case:
	//  Propagate adjoints from and to both INCLUSIVE
	static void propagateAdjointsMulti(
		Tape::iterator propagateFrom,
		Tape::iterator propagateTo)
	{
		auto it = propagateFrom;
		while (it != propagateTo)
		{
			it->propagateAll();
			--it;
		}
		it->propagateAll();
    }

    //  Operator overloading

    inline friend Number operator+(const Number& lhs, const Number& rhs)
    {
		const double e = lhs.value() + rhs.value();
        //  Eagerly evaluate and put on tape
        Number result(lhs.node(), rhs.node(), e);
        //  Eagerly compute derivatives
        result.lDer() = 1.0;
		result.rDer() = 1.0;

        return result;
    }
    inline friend Number operator+(const Number& lhs, const double& rhs)
    {
		const double e = lhs.value() + rhs;
        //  Eagerly evaluate and put on tape
		Number result(lhs.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = 1.0;

        return result;

    }
    inline friend Number operator+(const double& lhs, const Number& rhs)
    {
        return rhs + lhs;
    }

    inline friend Number operator-(const Number& lhs, const Number& rhs)
    {
		const double e = lhs.value() - rhs.value();
        //  Eagerly evaluate and put on tape
		Number result(lhs.node(), rhs.node(), e);
        //  Eagerly compute derivatives
		result.lDer() = 1.0;
		result.rDer() = -1.0;

        return result;
    }
    inline friend Number operator-(const Number& lhs, const double& rhs)
    {
		const double e = lhs.value() - rhs;
        //  Eagerly evaluate and put on tape
		Number result(lhs.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = 1.0;

        return result;

    }
    inline friend Number operator-(const double& lhs, const Number& rhs)
    {
		const double e = lhs - rhs.value();
        //  Eagerly evaluate and put on tape
		Number result(rhs.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = -1.0;

        return result;
    }

    inline friend Number operator*(const Number& lhs, const Number& rhs)
    {
		const double e = lhs.value() * rhs.value();
        //  Eagerly evaluate and put on tape
		Number result(lhs.node(), rhs.node(), e);
        //  Eagerly compute derivatives
		result.lDer() = rhs.value();
		result.rDer() = lhs.value();

        return result;
    }
    inline friend Number operator*(const Number& lhs, const double& rhs)
    {
		const double e = lhs.value() * rhs;
        //  Eagerly evaluate and put on tape
		Number result(lhs.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = rhs;

        return result;

    }
    inline friend Number operator*(const double& lhs, const Number& rhs)
    {
        return rhs * lhs;
    }

    inline friend Number operator/(const Number& lhs, const Number& rhs)
    {
		const double e = lhs.value() / rhs.value();
        //  Eagerly evaluate and put on tape
		Number result(lhs.node(), rhs.node(), e);
        //  Eagerly compute derivatives
        const double invRhs = 1.0 / rhs.value();
		result.lDer() = invRhs;
		result.rDer() = -lhs.value() * invRhs * invRhs;

        return result;
    }
    inline friend Number operator/(const Number& lhs, const double& rhs)
    {
		const double e = lhs.value() / rhs;
        //  Eagerly evaluate and put on tape
		Number result(lhs.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = 1.0 / rhs;

        return result;

    }
    inline friend Number operator/(const double& lhs, const Number& rhs)
    {
		const double e = lhs / rhs.value();
		//  Eagerly evaluate and put on tape
		Number result(rhs.node(), e);
		//  Eagerly compute derivatives
		result.derivative() = -lhs / rhs.value() / rhs.value();

        return result;
    }

	inline friend Number pow(const Number& lhs, const Number& rhs)
	{
		const double e = pow(lhs.value(), rhs.value());
		//  Eagerly evaluate and put on tape
		Number result(lhs.node(), rhs.node(), e);
		//  Eagerly compute derivatives
		result.lDer() = rhs.value() * e / lhs.value();
		result.rDer() = log(lhs.value()) * e;

		return result;
	}
	inline friend Number pow(const Number& lhs, const double& rhs)
	{
		const double e = pow(lhs.value(), rhs);
		//  Eagerly evaluate and put on tape
		Number result(lhs.node(), e);
		//  Eagerly compute derivatives
		result.derivative() = rhs * e / lhs.value();

		return result;
	}
	inline friend Number pow(const double& lhs, const Number& rhs)
	{
		const double e = pow(lhs, rhs.value());
		//  Eagerly evaluate and put on tape
		Number result(rhs.node(), e);
		//  Eagerly compute derivatives
		result.derivative() = log(lhs) * e;

		return result;
	}

	inline friend Number max(const Number& lhs, const Number& rhs)
	{
		const bool lmax = lhs.value() > rhs.value();
		//  Eagerly evaluate and put on tape
		Number result(lhs.node(), rhs.node(), lmax? lhs.value() : rhs.value());
		//  Eagerly compute derivatives
		if (lmax)
		{
			result.lDer() = 1.0;
			result.rDer() = 0.0;
		}
		else
		{
			result.lDer() = 0.0;
			result.rDer() = 1.0;
		}

		return result;
	}
	inline friend Number max(const Number& lhs, const double& rhs)
	{
		const bool lmax = lhs.value() > rhs;
		//  Eagerly evaluate and put on tape
		Number result(lhs.node(), lmax ? lhs.value() : rhs);
		//  Eagerly compute derivatives
		result.derivative() = lmax ? 1.0 : 0.0;

		return result;
	}
	inline friend Number max(const double& lhs, const Number& rhs)
	{
		const bool rmax = rhs.value() > lhs;
		//  Eagerly evaluate and put on tape
		Number result(rhs.node(), rmax ? rhs.value() : lhs);
		//  Eagerly compute derivatives
		result.derivative() = rmax ? 1.0 : 0.0;

		return result;
	}

	inline friend Number min(const Number& lhs, const Number& rhs)
	{
		const bool lmin = lhs.value() < rhs.value();
		//  Eagerly evaluate and put on tape
		Number result(lhs.node(), rhs.node(), lmin ? lhs.value() : rhs.value());
		//  Eagerly compute derivatives
		if (lmin)
		{
			result.lDer() = 1.0;
			result.rDer() = 0.0;
		}
		else
		{
			result.lDer() = 0.0;
			result.rDer() = 1.0;
		}

		return result;
	}
	inline friend Number min(const Number& lhs, const double& rhs)
	{
		const bool lmin = lhs.value() < rhs;
		//  Eagerly evaluate and put on tape
		Number result(lhs.node(), lmin ? lhs.value() : rhs);
		//  Eagerly compute derivatives
		result.derivative() = lmin ? 1.0 : 0.0;

		return result;
	}
	inline friend Number min(const double& lhs, const Number& rhs)
	{
		const bool rmin = rhs.value() < lhs;
        //  Eagerly evaluate and put on tape
		Number result(rhs.node(), rmin ? rhs.value() : lhs);
        //  Eagerly compute derivatives
		result.derivative() = rmin ? 1.0 : 0.0;

        return result;
    }

    Number& operator+=(const Number& arg)
    {
        *this = *this + arg;
        return *this;
    }
    Number& operator+=(const double& arg)
    {
        *this = *this + arg;
        return *this;
    }

    Number& operator-=(const Number& arg)
    {
        *this = *this - arg;
        return *this;
    }
    Number& operator-=(const double& arg)
    {
        *this = *this - arg;
        return *this;
    }

    Number& operator*=(const Number& arg)
    {
        *this = *this * arg;
        return *this;
    }
    Number& operator*=(const double& arg)
    {
        *this = *this * arg;
        return *this;
    }

    Number& operator/=(const Number& arg)
    {
        *this = *this / arg;
        return *this;
    }
    Number& operator/=(const double& arg)
    {
        *this = *this / arg;
        return *this;
    }

    //  Unary +/-
    Number operator-() const
    {
        return 0.0 - *this;
    }
    Number operator+() const 
    {
        return *this;
    }

    //  Overloading continued, unary functions

    inline friend Number exp(const Number& arg)
    {
        const double e = exp(arg.value());
        //  Eagerly evaluate and put on tape
        Number result(arg.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = e;

        return result;
    }

    inline friend Number log(const Number& arg)
    {
		const double e = log(arg.value());
        //  Eagerly evaluate and put on tape
		Number result(arg.node(), e);
        //  Eagerly compute derivatives
        result.derivative() = 1.0 / arg.value();

        return result;
    }

    inline friend Number sqrt(const Number& arg)
    {
        const double e = sqrt(arg.value());
        //  Eagerly evaluate and put on tape
		Number result(arg.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = 0.5 / e;

        return result;
    }

    inline friend Number fabs(const Number& arg)
    {
        const double e = fabs(arg.value());
        //  Eagerly evaluate and put on tape
		Number result(arg.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = arg.value() > 0.0 ? 1.0 : -1.0;

        return result;
    }

    inline friend Number normalDens(const Number& arg)
    {
        const double e = normalDens(arg.value());
        //  Eagerly evaluate and put on tape
		Number result(arg.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = - arg.value() * e;

        return result;
    }

    inline friend Number normalCdf(const Number& arg)
    {
        const double e = normalCdf(arg.value());
        //  Eagerly evaluate and put on tape
		Number result(arg.node(), e);
        //  Eagerly compute derivatives
		result.derivative() = normalDens(arg.value());

        return result;
    }

    //  Finally, comparison

    inline friend bool operator==(const Number& lhs, const Number& rhs)
    {
        return lhs.value() == rhs.value();
    }
    inline friend bool operator==(const Number& lhs, const double& rhs)
    {
        return lhs.value() == rhs;
    }
    inline friend bool operator==(const double& lhs, const Number& rhs)
    {
        return lhs == rhs.value();
    }

    inline friend bool operator!=(const Number& lhs, const Number& rhs)
    {
        return lhs.value() != rhs.value();
    }
    inline friend bool operator!=(const Number& lhs, const double& rhs)
    {
        return lhs.value() != rhs;
    }
    inline friend bool operator!=(const double& lhs, const Number& rhs)
    {
        return lhs != rhs.value();
    }

    inline friend bool operator<(const Number& lhs, const Number& rhs)
    {
        return lhs.value() < rhs.value();
    }
    inline friend bool operator<(const Number& lhs, const double& rhs)
    {
        return lhs.value() < rhs;
    }
    inline friend bool operator<(const double& lhs, const Number& rhs)
    {
        return lhs < rhs.value();
    }

    inline friend bool operator>(const Number& lhs, const Number& rhs)
    {
        return lhs.value() > rhs.value();
    }
    inline friend bool operator>(const Number& lhs, const double& rhs)
    {
        return lhs.value() > rhs;
    }
    inline friend bool operator>(const double& lhs, const Number& rhs)
    {
        return lhs > rhs.value();
    }

    inline friend bool operator<=(const Number& lhs, const Number& rhs)
    {
        return lhs.value() <= rhs.value();
    }
    inline friend bool operator<=(const Number& lhs, const double& rhs)
    {
        return lhs.value() <= rhs;
    }
    inline friend bool operator<=(const double& lhs, const Number& rhs)
    {
        return lhs <= rhs.value();
    }

    inline friend bool operator>=(const Number& lhs, const Number& rhs)
    {
        return lhs.value() >= rhs.value();
    }
    inline friend bool operator>=(const Number& lhs, const double& rhs)
    {
        return lhs.value() >= rhs;
    }
    inline friend bool operator>=(const double& lhs, const Number& rhs)
    {
        return lhs >= rhs.value();
    }
};


