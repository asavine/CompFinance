
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

//  Implementation of AAD with expression templates 
//  (AADET, chapter 15)

//  Defines expressions and the Number type

#include <algorithm>
#include "AADTape.h"

//  Base CRTP expression class 
//      Note: overloaded operators catch all expressions and nothing else
template <class E>
struct Expression
{
    //  CRTP "virtualization"
    double value() const 
    { 
        return static_cast<const E*>(this)->value(); 
    }
	
    //  Just another interface
    explicit operator double() const 
    { 
        return value(); 
    }
};

//  Note that Number is a leaf expression
//  Defined in the bottom of the file

//  Binary expression
//  LHS : the expression on the left
//  RHS : the expression on the right
//  OP : the binary operator
template <class LHS, class RHS, class OP>
class BinaryExpression 
    //  CRTP
    : public Expression<BinaryExpression<LHS, RHS, OP>>
{
    const double myValue;

    const LHS lhs;
    const RHS rhs;

public:

    //  Constructor out of 2 expressions
    //  Note: eager evaluation on construction
    explicit BinaryExpression(
        const Expression<LHS>& l,
        const Expression<RHS>& r)
        : myValue(OP::eval(l.value(), r.value())), 
        lhs(static_cast<const LHS&>(l)), 
        rhs(static_cast<const RHS&>(r))
    {}

    //  Value accessors
    double value() const { return myValue; }

    //	Expression template magic
    //  Expressions know
    //  AT COMPILE TIME
    //  the number of active inputs in their sub-expressions
    enum { numNumbers = LHS::numNumbers + RHS::numNumbers };

    //  Push adjoint down the expression
    //  N : total number of active inputs in the expression
    //  n : number of active inputs already processed
    template <size_t N, size_t n>
    void pushAdjoint(
        //  Node for the complete expression being processed
        Node&		exprNode,
        //  Adjoint cumulated for this binary node, or 1 if top
        const double	adjoint)       
        const
    {
        //  Push on the left
        if (LHS::numNumbers > 0)
        {
            lhs.pushAdjoint<N, n>(
                exprNode, 
                adjoint * OP::leftDerivative(lhs.value(), rhs.value(), value()));
        }

        //  Push on the right
        if (RHS::numNumbers > 0)
        {
            //  Note left push processed LHS::numNumbers numbers
            //  So the next number to be processed is n + LHS::numNumbers
            rhs.pushAdjoint<N, n + LHS::numNumbers>(
                exprNode, 
                adjoint * OP::rightDerivative(lhs.value(), rhs.value(), value()));
        }
    }
};

//  "Concrete" binaries, we only need to define operations and derivatives
struct OPMult
{
    static const double eval(const double l, const double r) 
    { 
        return l * r; 
    }
    
    static const double leftDerivative
        (const double l, const double r, const double v) 
    { 
        return r; 
    }
    
    static const double rightDerivative
        (const double l, const double r, const double v) 
    { 
        return l; 
    }
};

struct OPAdd
{
    static const double eval(const double l, const double r)
    { 
        return l + r; 
    }
    
    static const double leftDerivative
        (const double l, const double r, const double v)
    { 
        return 1.0; 
    }
    
    static const double rightDerivative
        (const double l, const double r, const double v)
    { 
        return 1.0; 
    }
};

struct OPSub
{
    static const double eval(const double l, const double r)
    {
        return l - r;
    }

    static const double leftDerivative
    (const double l, const double r, const double v)
    {
        return 1.0;
    }

    static const double rightDerivative
    (const double l, const double r, const double v)
    {
        return -1.0;
    }
};

struct OPDiv
{
    static const double eval(const double l, const double r)
    {
        return l / r;
    }

    static const double leftDerivative
    (const double l, const double r, const double v)
    {
        return 1.0 / r;
    }

    static const double rightDerivative
    (const double l, const double r, const double v)
    {
        return -l / r / r;
    }
};

struct OPPow
{
    static const double eval(const double l, const double r)
    {
        return pow(l, r);
    }

    static const double leftDerivative
    (const double l, const double r, const double v)
    {
        return r*v / l;
    }

    static const double rightDerivative
    (const double l, const double r, const double v)
    {
        return log(l)*v;
    }
};

struct OPMax
{
    static const double eval(const double l, const double r)
    {
        return max(l, r);
    }

    static const double leftDerivative
    (const double l, const double r, const double v)
    {
        return l > r ? 1.0 : 0.0;
    }

    static const double rightDerivative
    (const double l, const double r, const double v)
    {
        return r > l? 1.0 : 0.0;
    }
};

struct OPMin
{
    static const double eval(const double l, const double r)
    {
        return min(l, r);
    }

    static const double leftDerivative
    (const double l, const double r, const double v)
    {
        return l < r ? 1.0 : 0.0;
    }

    static const double rightDerivative
    (const double l, const double r, const double v)
    {
        return r < l ? 1.0 : 0.0;
    }
};

//  Operator overloading for binary expressions
//      build the corresponding expressions

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPMult> operator*(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs) 
{ 
    return BinaryExpression<LHS, RHS, OPMult>(lhs, rhs); 
}

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPAdd> operator+(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPAdd>(lhs, rhs);
}

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPSub> operator-(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPSub>(lhs, rhs);
}

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPDiv> operator/(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPDiv>(lhs, rhs);
}

template <class LHS, class RHS>
 BinaryExpression<LHS, RHS, OPPow> pow(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPPow>(lhs, rhs);
}

template <class LHS, class RHS>
BinaryExpression<LHS, RHS, OPMax> max(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPMax>(lhs, rhs);
}

template <class LHS, class RHS>
BinaryExpression<LHS, RHS, OPMin> min(
    const Expression<LHS>& lhs, const Expression<RHS>& rhs)
{
    return BinaryExpression<LHS, RHS, OPMin>(lhs, rhs);
}

//  Unary expressions : Same logic with one argument

//  The CRTP class
template <class ARG, class OP>
class UnaryExpression
    //  CRTP
    : public Expression<UnaryExpression<ARG, OP>>
{
    const double myValue;

    const ARG arg;

    //  For binary operators with a double on one side
    //      we store the double
    const double dArg = 0.0;   

public:

    //  Constructor 
    //  Note : eager evaluation on construction
    explicit UnaryExpression(
        const Expression<ARG>& a)
        : myValue(OP::eval(a.value(), 0.0)), arg(static_cast<const ARG&>(a)) {}

    //  Special constructor for binary expressions with a double on one side
    explicit UnaryExpression(
        const Expression<ARG>& a,
        const double b)
        : myValue(OP::eval(a.value(), b)), arg(static_cast<const ARG&>(a)), dArg(b) {}

    //  Value accessors
    double value() const { return myValue; }

    //	Expression template magic
    enum { numNumbers = ARG::numNumbers };

    //  Push adjoint down the expression 
    template <size_t N, size_t n>
    void pushAdjoint(
        //  Node for the complete expression being processed
        Node&		exprNode,
        //  Adjoint cumulated on the node, 1 if top
        const double	adjoint)       
        const
    {
        //  Push into argument
        if (ARG::numNumbers > 0)
        {
            arg.pushAdjoint<N, n>(
                exprNode, 
                adjoint * OP::derivative(arg.value(), value(), dArg));
        }
    }
};

//  The unary operators

struct OPExp
{
    static const double eval(const double r, const double d) 
    { 
        return exp(r); 
    }
    
    static const double derivative
        (const double r, const double v, const double d)
    { 
        return v; 
    }
};

struct OPLog
{
    static const double eval(const double r, const double d)
    { 
        return log(r); 
    }
    
    static const double derivative
        (const double r, const double v, const double d)
    { 
        return 1.0 / r; 
    }
};

struct OPSqrt
{
    static const double eval(const double r, const double d)
    { 
        return sqrt(r); 
    }

    static const double derivative
        (const double r, const double v, const double d)
    { 
        return 0.5 / v; 
    }
};

struct OPFabs
{
    static const double eval(const double r, const double d)
    {
        return fabs(r);
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return r > 0.0 ? 1.0 : -1.0;
    }
};

struct OPNormalDens
{
    static const double eval(const double r, const double d)
    {
        return normalDens(r);
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return - r * v;
    }
};

struct OPNormalCdf
{
    static const double eval(const double r, const double d)
    {
        return normalCdf(r);
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return normalDens(r);
    }
};

//  Binary operators with a double on one side 

//  * double or double *
struct OPMultD
{
    static const double eval(const double r, const double d)
    {
        return r * d;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return d;
    }
};

//  + double or double +
struct OPAddD
{
    static const double eval(const double r, const double d)
    {
        return r + d;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return 1.0;
    }
};

//  double -
struct OPSubDL
{
    static const double eval(const double r, const double d)
    {
        return d - r;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return -1.0;
    }
};

//  - double
struct OPSubDR
{
    static const double eval(const double r, const double d)
    {
        return r - d;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return 1.0;
    }
};

//  double /
struct OPDivDL
{
    static const double eval(const double r, const double d)
    {
        return d / r;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return -d / r / r;
    }
};

//  / double
struct OPDivDR
{
    static const double eval(const double r, const double d)
    {
        return r / d;
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return 1.0 / d;
    }
};

//  pow (d,)
struct OPPowDL
{
    static const double eval(const double r, const double d)
    {
        return pow(d, r);
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return log(d) * v;
    }
};

//  pow (,d)
struct OPPowDR
{
    static const double eval(const double r, const double d)
    {
        return pow(r, d);
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return d * v / r;
    }
};

//  max (d,)
struct OPMaxD
{
    static const double eval(const double r, const double d)
    {
        return max(r, d);
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return r > d ? 1.0 : 0.0;
    }
};

//  min (d,)
struct OPMinD
{
    static const double eval(const double r, const double d)
    {
        return min(r, d);
    }

    static const double derivative
    (const double r, const double v, const double d)
    {
        return r < d ? 1.0 : 0.0;
    }
};

//  And overloading

template <class ARG>
 UnaryExpression<ARG, OPExp> exp(const Expression<ARG>& arg)
{ 
    return UnaryExpression<ARG, OPExp>(arg);
}

template <class ARG>
 UnaryExpression<ARG, OPLog> log(const Expression<ARG>& arg)
{
    return UnaryExpression<ARG, OPLog>(arg);
}

template <class ARG>
 UnaryExpression<ARG, OPSqrt> sqrt(const Expression<ARG>& arg)
{
    return UnaryExpression<ARG, OPSqrt>(arg);
}

template <class ARG>
UnaryExpression<ARG, OPFabs> fabs(const Expression<ARG>& arg)
{
    return UnaryExpression<ARG, OPFabs>(arg);
}

template <class ARG>
UnaryExpression<ARG, OPNormalDens> normalDens(const Expression<ARG>& arg)
{
    return UnaryExpression<ARG, OPNormalDens>(arg);
}

template <class ARG>
UnaryExpression<ARG, OPNormalCdf> normalCdf(const Expression<ARG>& arg)
{
    return UnaryExpression<ARG, OPNormalCdf>(arg);
}

//  Overloading continued,
//      binary operators with a double on one side 

template <class ARG>
 UnaryExpression<ARG, OPMultD> operator*(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPMultD>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPMultD> operator*(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPMultD>(lhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPAddD> operator+(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPAddD>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPAddD> operator+(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPAddD>(lhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPSubDL> operator-(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPSubDL>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPSubDR> operator-(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPSubDR>(lhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPDivDL> operator/(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPDivDL>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPDivDR> operator/(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPDivDR>(lhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPPowDL> pow(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPPowDL>(rhs, d);
}

template <class ARG>
 UnaryExpression<ARG, OPPowDR> pow(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPPowDR>(lhs, d);
}


template <class ARG>
UnaryExpression<ARG, OPMaxD> max(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPMaxD>(rhs, d);
}

template <class ARG>
UnaryExpression<ARG, OPMaxD> max(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPMaxD>(lhs, d);
}

template <class ARG>
UnaryExpression<ARG, OPMinD> min(
    const double d, const Expression<ARG>& rhs)
{
    return UnaryExpression<ARG, OPMinD>(rhs, d);
}

template <class ARG>
UnaryExpression<ARG, OPMinD> min(
    const Expression<ARG>& lhs, const double d)
{
    return UnaryExpression<ARG, OPMinD>(lhs, d);
}

//  Comparison, same as traditional

template<class E, class F>
 bool operator==(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() == rhs.value();
}
template<class E>
 bool operator==(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() == rhs;
}
template<class E>
 bool operator==(const double& lhs, const Expression<E>& rhs)
{
    return lhs == rhs.value();
}

template<class E, class F>
 bool operator!=(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() != rhs.value();
}
template<class E>
 bool operator!=(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() != rhs;
}
template<class E>
 bool operator!=(const double& lhs, const Expression<E>& rhs)
{
    return lhs != rhs.value();
}

template<class E, class F>
 bool operator<(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() < rhs.value();
}
template<class E>
 bool operator<(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() < rhs;
}
template<class E>
 bool operator<(const double& lhs, const Expression<E>& rhs)
{
    return lhs < rhs.value();
}

template<class E, class F>
 bool operator>(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() > rhs.value();
}
template<class E>
 bool operator>(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() > rhs;
}
template<class E>
 bool operator>(const double& lhs, const Expression<E>& rhs)
{
    return lhs > rhs.value();
}

template<class E, class F>
 bool operator<=(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() <= rhs.value();
}
template<class E>
 bool operator<=(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() <= rhs;
}
template<class E>
 bool operator<=(const double& lhs, const Expression<E>& rhs)
{
    return lhs <= rhs.value();
}

template<class E, class F>
 bool operator>=(const Expression<E>& lhs, const Expression<F>& rhs)
{
    return lhs.value() >= rhs.value();
}
template<class E>
 bool operator>=(const Expression<E>& lhs, const double& rhs)
{
    return lhs.value() >= rhs;
}
template<class E>
 bool operator>=(const double& lhs, const Expression<E>& rhs)
{
    return lhs >= rhs.value();
}

//  Finally, unary +/- operators

template <class RHS>
UnaryExpression<RHS, OPSubDL> operator-
(const Expression<RHS>& rhs)
{
    return 0.0 - rhs;
}

template <class RHS>
Expression<RHS> operator+
(const Expression<RHS>& rhs)
{
    return rhs;
}

//  The Number type, also an expression

class Number : public Expression<Number>
{
    //  The value and node for this number, same as traditional
    double		myValue;
    Node*	    myNode;

    //  Node creation on tape

    template <size_t N>
    Node* createMultiNode()
    {
        return tape->recordNode<N>();
    }

    //  Flattening:
    //      This is where, on assignment or construction from an expression,
    //      that derivatives are pushed through the expression's DAG 
    template<class E>
    void fromExpr(
        //  RHS expression, will be flattened into this Number
        const Expression<E>& e)
    {
        //  Build expression node on tape
        auto* node = createMultiNode<E::numNumbers>();
        
        //  Push adjoints through expression with adjoint = 1 on top
        static_cast<const E&>(e).pushAdjoint<E::numNumbers, 0>(*node, 1.0);

        //  Set my node
        myNode = node;
    }

public:

    //  Expression template magic
    enum { numNumbers = 1 };

    //  Push adjoint
    //  Numbers are expression leaves, 
    //      pushAdjoint() receives their adjoint in the expression
    //  Numbers don't "push" anything, they register their derivatives on tape
    template <size_t N, size_t n>
    void pushAdjoint(
        //  Node for the complete expression
        Node&		    exprNode,	
        //  Adjoint accumulated for this number, in the expression
        const double	adjoint)	
        const
    {
        //  adjoint = d (expression) / d (thisNumber) : register on tape
        //  note n: index of this number on the node on tape

        //  Register adjoint
        exprNode.pAdjPtrs[n] = Tape::multi? myNode->pAdjoints : &myNode->mAdjoint;
		
        //  Register derivative
        exprNode.pDerivatives[n] = adjoint;
    }

    //  Static access to tape, same as traditional
    static thread_local Tape* tape;

    //  Constructors

    Number() {}

    explicit Number(const double val) : myValue(val)
    {
        //  Create leaf
        myNode = createMultiNode<0>();
    }

    Number& operator=(const double val)
    {
        myValue = val;
        //  Create leaf
        myNode = createMultiNode<0>();
        return *this;
    }

    //  No need for copy and assignment
    //  Default ones do the right thing:
    //      copy value and pointer to node on tape

    //  Construct or assign from expression
    
    template <class E>
    Number(const Expression<E>& e) : myValue(e.value())
    {
        //  Flatten RHS expression
        fromExpr<E>(static_cast<const E&>(e));
    }

    template <class E>
    Number& operator=
        (const Expression<E>& e) 
    {
        myValue = e.value();
        //  Flatten RHS expression
        fromExpr<E>(static_cast<const E&>(e));
        return *this;
    }

    //  Explicit coversion to double
    explicit operator double& () { return myValue; }
    explicit operator double () const { return myValue; }

    //  All the normal accessors and propagators, same as traditional 
    
    //  Put on tape
    void putOnTape()
    {
        myNode = createMultiNode<0>();
    }

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
        auto it = tape->find(myNode);
        //  Reverse and propagate until we hit the stop
        while (it != propagateTo)
        {
            it->propagateOne();
            --it;
        }
        it->propagateOne();
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

    //  Multi-adjoint propagation

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

    //  Unary operators

    template <class E>
    Number& operator+=(const Expression<E>& e)
    {
        *this = *this + e;
        return *this;
    }

    template <class E>
    Number& operator*=(const Expression<E>& e)
    {
        *this = *this * e;
        return *this;
    }

    template <class E>
    Number& operator-=(const Expression<E>& e)
    {
        *this = *this - e;
        return *this;
    }

    template <class E>
    Number& operator/=(const Expression<E>& e)
    {
        *this = *this / e;
        return *this;
    }
};

