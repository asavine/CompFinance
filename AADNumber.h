#pragma once

#include "AADTape.h"

class Number
{
    double myValue;
    Node* myNode;

    //  Node creation on tape

    template <class NodeType>
    inline void createNode()
    {
        //  Placement syntax to allocate in place on tape
        myNode = new (tape->allocate(sizeof(NodeType))) NodeType;
    }

    //  Access to node for friends
    template <class NodeType = Node>
    inline NodeType* node() const
    {
        return reinterpret_cast<NodeType*>(myNode);
    }

    //  To help overloading
    static const struct LeafType {} leaf;
    static const struct UnaryType {} unary;
    static const struct BinaryType {} binary;

    Number(const double val, UnaryType) :
        myValue(val)
    {
        createNode<UnaryNode>();
    }

    Number(const double val, BinaryType) :
        myValue(val)
    {
        createNode<BinaryNode>();
    }

public:

    //  Static access to tape
    static thread_local Tape* tape;

    //  Constructors

    Number() {}

    explicit Number(const double val) :
        myValue(val)
    {
        createNode<Leaf>();
    }

    //  Assignments

    Number& operator=(const double val)
    {
        myValue = val;
        createNode<Leaf>();

        return *this;
    }

    //  Put on tape
    void putOnTape()
    {
        createNode<Leaf>();
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
    double& adjoint()
    {
        return myNode->adjoint;
    }
    double adjoint() const
    {
        return myNode->adjoint;
    }

/*
    //  This may be dangerous, 
    //  but necessary for templated code to compile
    operator double() const
    {
        return myValue;
    }
*/

    //  Propagation

    //  Reset all adjoints on the tape
    static void resetAdjoints()
    {
        for (Node& node : *tape) node.adjoint = 0.0;
    }
    //  Propagate adjoints
    //      from and to both INCLUSIVE
    static void propagateAdjoints(
        Tape::iterator propagateFrom,
        Tape::iterator propagateTo)
    {
        auto it = propagateFrom;
        while (it != propagateTo)
        {
            it->propagate();
            --it;
        }
        it->propagate();
    }

    //  Convenient overloads

    //  Set the adjoint on this node to 1,
    //  Then propagate from the node
    void propagateAdjoints(
        //  We start on this number's node
        Tape::iterator propagateTo,
        //  reset adjoints first?
        const bool reset = false)
    {
        //  Reset
        if (reset) resetAdjoints();
        //  Set this adjoint to 1
        adjoint() = 1.0;
        //  Find node on tape
        auto it = tape->find(myNode);
        //  Reverse and propagate until we hit the stop
        while (it != propagateTo)
        {
            it->propagate();
            --it;
        }
        it->propagate();
    }

    //  These 2 set the adjoint to 1 on this node
    void propagateToStart(
        const bool reset = false)
    {
        propagateAdjoints(tape->begin(), reset);
    }
    void propagateToMark(
        const bool reset = false)
    {
        propagateAdjoints(tape->markIt(), reset);
    }

    //  This one leaves the adjoints untouched
    static void propagateMarkToStart()
    {
        propagateAdjoints(tape->markIt(), tape->begin());
    }

    //  Operator overloading

    inline friend Number operator+(const Number& lhs, const Number& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs.value() + rhs.value(), binary);
        //  Set arguments
        result.node<BinaryNode>()->arguments[0] = lhs.node();
        result.node<BinaryNode>()->arguments[1] = rhs.node();
        //  Eagerly compute derivatives
        result.node<BinaryNode>()->derivatives[0] = 1.0;
        result.node<BinaryNode>()->derivatives[1] = 1.0;

        return result;
    }
    inline friend Number operator+(const Number& lhs, const double& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs.value() + rhs, unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = lhs.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = 1.0;

        return result;

    }
    inline friend Number operator+(const double& lhs, const Number& rhs)
    {
        return rhs + lhs;
    }

    inline friend Number operator-(const Number& lhs, const Number& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs.value() - rhs.value(), binary);
        //  Set arguments
        result.node<BinaryNode>()->arguments[0] = lhs.node();
        result.node<BinaryNode>()->arguments[1] = rhs.node();
        //  Eagerly compute derivatives
        result.node<BinaryNode>()->derivatives[0] = 1.0;
        result.node<BinaryNode>()->derivatives[1] = -1.0;

        return result;
    }
    inline friend Number operator-(const Number& lhs, const double& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs.value() - rhs, unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = lhs.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = 1.0;

        return result;

    }
    inline friend Number operator-(const double& lhs, const Number& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs - rhs.value(), unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = rhs.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = -1.0;

        return result;
    }

    inline friend Number operator*(const Number& lhs, const Number& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs.value() * rhs.value(), binary);
        //  Set arguments
        result.node<BinaryNode>()->arguments[0] = lhs.node();
        result.node<BinaryNode>()->arguments[1] = rhs.node();
        //  Eagerly compute derivatives
        result.node<BinaryNode>()->derivatives[0] = rhs.value();
        result.node<BinaryNode>()->derivatives[1] = lhs.value();

        return result;
    }
    inline friend Number operator*(const Number& lhs, const double& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs.value() * rhs, unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = lhs.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = rhs;

        return result;

    }
    inline friend Number operator*(const double& lhs, const Number& rhs)
    {
        return rhs * lhs;
    }

    inline friend Number operator/(const Number& lhs, const Number& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs.value() / rhs.value(), binary);
        //  Set arguments
        result.node<BinaryNode>()->arguments[0] = lhs.node();
        result.node<BinaryNode>()->arguments[1] = rhs.node();
        //  Eagerly compute derivatives
        const double invRhs = 1.0 / rhs.value();
        result.node<BinaryNode>()->derivatives[0] = 1.0 * invRhs;
        result.node<BinaryNode>()->derivatives[1] = -lhs.value() * invRhs * invRhs;

        return result;
    }
    inline friend Number operator/(const Number& lhs, const double& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs.value() / rhs, unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = lhs.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = 1.0 / rhs;

        return result;

    }
    inline friend Number operator/(const double& lhs, const Number& rhs)
    {
        //  Eagerly evaluate and put on tape
        Number result(lhs / rhs.value(), unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = rhs.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = -lhs / rhs.value() / rhs.value();
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
    Number operator-()
    {
        return 0.0 - *this;
    }
    Number operator+()
    {
        return *this;
    }

    //  Unary functions
    inline friend Number exp(const Number& arg)
    {
        const double e = exp(arg.value());
        //  Eagerly evaluate and put on tape
        Number result(e, unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = arg.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = e;

        return result;
    }

    inline friend Number log(const Number& arg)
    {
        //  Eagerly evaluate and put on tape
        Number result(log(arg.value()), unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = arg.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = 1.0 / arg.value();

        return result;
    }

    inline friend Number sqrt(const Number& arg)
    {
        const double e = sqrt(arg.value());
        //  Eagerly evaluate and put on tape
        Number result(e, unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = arg.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = 0.5 / e;

        return result;
    }

    //  Binary functions
    inline friend Number pow(const Number& lhs, const Number& rhs)
    {
        const double e = pow(lhs.value(), rhs.value());
        //  Eagerly evaluate and put on tape
        Number result(e, unary);
        //  Set arguments
        result.node<BinaryNode>()->arguments[0] = lhs.node();
        result.node<BinaryNode>()->arguments[1] = rhs.node();
        //  Eagerly compute derivatives
        result.node<BinaryNode>()->derivatives[0] = rhs.value() * e / lhs.value();
        result.node<BinaryNode>()->derivatives[1] = log(lhs.value()) * e;

        return result;
    }
    inline friend Number pow(const Number& lhs, const double& rhs)
    {
        const double e = pow(lhs.value(), rhs);
        //  Eagerly evaluate and put on tape
        Number result(e, unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = lhs.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = rhs * e / lhs.value();

        return result;
    }
    inline friend Number pow(const double& lhs, const Number& rhs)
    {
        const double e = pow(lhs, rhs.value());
        //  Eagerly evaluate and put on tape
        Number result(e, unary);
        //  Set arguments
        result.node<UnaryNode>()->argument = rhs.node();
        //  Eagerly compute derivatives
        result.node<UnaryNode>()->derivative = log(lhs) * e;

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

template <class IT>
inline void putOnTape(IT begin, IT end)
{
    for (auto it = begin; it != end; ++it)
        it->putOnTape();
}

template<class To, class From>
struct Convert;

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

template<>
struct Convert<double, Number>
{
    static inline double convert(const Number from)
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
    using destType = remove_reference<decltype(*destBegin)>::type;

    while (srcBegin != srcEnd)
    {
        *destBegin++ = convert<destType>(*srcBegin++);
    } 
}

