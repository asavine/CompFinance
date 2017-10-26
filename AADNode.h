#pragma once

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



