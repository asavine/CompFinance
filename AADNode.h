
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

#include <exception>
using namespace std;

struct Node
{
    double adjoint = 0.0;

    virtual void propagate() {}
    virtual ~Node() {}
};

struct Leaf : public Node
{
    
};

struct UnaryNode : public Node
{
    double derivative;
    Node* argument;

    void propagate() override
    {
        if (adjoint == 0.0) return;
        argument->adjoint += derivative * adjoint;
    }
};

struct BinaryNode : public Node
{
    double derivatives[2];
    Node* arguments[2];

    void propagate() override
    {
        if (adjoint == 0.0) return;
        arguments[0]->adjoint += derivatives[0] * adjoint;
        arguments[1]->adjoint += derivatives[1] * adjoint;
    }
};

template <size_t N>
struct MultiNode : public Node
{
    double derivatives[N];
    Node* arguments[N];

    void propagate() override
    {
        if (adjoint == 0.0) return;

        //  The compiler should unroll the loop
        for (int i = 0; i < N; ++i)
        {
            arguments[i]->adjoint += derivatives[i] * adjoint;
        }
    }
};

//  Specialization for leaf
template<>
struct MultiNode<0> : public Node
{

};

