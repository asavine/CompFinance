#pragma once

#include <vector>
#include <memory>

using namespace std;

using Time = double;

extern Time systemTime;

//  The scenario is all the model data needs on a given event date
//  In this simple implementation, a scenario is just a spot price
template <class T>
struct scenario
{
    T   spot;
};

template <class T>
class Product
{
public:

    virtual const vector<Time>& timeline() const = 0;

    virtual T payoff(const vector<scenario<T>>& path) const = 0;

    virtual unique_ptr<Product<T>> clone() const = 0;

    virtual ~Product() {}
};

template <class T>
class Model
{
public:

    virtual void init(const vector<Time>& productTimeline) = 0;

    virtual size_t simDim() const = 0;

    virtual void generatePath(const vector<double>& gaussVec, vector<scenario<T>>& path) const = 0;

    virtual unique_ptr<Model<T>> clone() const = 0;

    virtual ~Model() {}
};

class RNG
{
public:
    
    virtual void init(const size_t simDim) = 0;

    virtual const vector<double>& nextG() = 0;

    virtual unique_ptr<RNG> clone() const = 0;

    virtual ~RNG() {}
};

//	MC simulator: free function that conducts simulations and returns a vector of nPath payoffs
template <class T>
inline vector<T> mcSimul(
    const Product<T>& prd,
    //  Note model and rng are not const because they set themselves up for the simulation
    //  This is where people would use the mutable keyword, something we learned never to indulge in
    Model<T>& mdl,          
    RNG& rng,			            
    const size_t nPath)                      
{
    vector<T> res(nPath);	                            //	Allocate results
    mdl.init(prd.timeline());                           //  Init the simulation timeline
    rng.init(mdl.simDim());                             //  Init the RNG
    vector<scenario<T>> path(prd.timeline().size());    //  Allocate path

    //	Iterate through paths	
    for (size_t i = 0; i<nPath; i++)
    {
        const auto& gaussRanVec = rng.nextG();  //  Next Gaussian vector, dimension D
        mdl.generatePath(gaussRanVec, path);    //  Generate path, consume Gaussian vector
        res[i] = prd.payoff(path);              //	Compute result
    }

    return res;	//	C++11: move
}