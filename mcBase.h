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

    //  Access to the product timeline
    virtual const vector<Time>& timeline() const = 0;

    //  Compute payoff given a path (on the product timeline)
    virtual T payoff(const vector<scenario<T>>& path) const = 0;

    virtual unique_ptr<Product<T>> clone() const = 0;

    virtual ~Product() {}
};

template <class T>
class Model
{
public:

    //  Initialize with product timeline
    virtual void init(const vector<Time>& productTimeline) = 0;

    //  Access to the MC dimension
    virtual size_t simDim() const = 0;

    //  Generate a path consuming a vector[simDim()] of independent Gaussians
    //  The path vector filled by the function must be pre-allocated
    virtual void generatePath(const vector<double>& gaussVec, vector<scenario<T>>& path) const = 0;

    virtual unique_ptr<Model<T>> clone() const = 0;

    virtual ~Model() {}
};

class RNG
{
public:
    
    //  Initialise with dimension simDim
    virtual void init(const size_t simDim) = 0;

    //  Compute the next vector[simDim] of independent Gaussians
    //  The vector is filled by the function and must be pre-allocated
    virtual void nextG(vector<double>& gaussVec) = 0;

    virtual unique_ptr<RNG> clone() const = 0;

    virtual ~RNG() {}

    //  Access dimension
    virtual size_t simDim() const = 0;

    //  Skip ahead
    virtual void skipTo(const long b)
    {
        vector<double> dummy(simDim());
        for (int i = 0; i < b; ++i) nextG(dummy);
    }
};

//	MC simulator: free function that conducts simulations 
//      and returns a vector of nPath payoffs
template <class T>
inline vector<T> mcSimul(
    const Product<T>& prd,
    const Model<T>& mdl,          
    const RNG& rng,			            
    const size_t nPath,
    const bool antithetic)                      
{
    //  Work with copies of the model and RNG
    //      which are modified when we set up the simulation
    //  Copies are OK at high level
    auto cMdl = mdl.clone();
    auto cRng = rng.clone();

    vector<T> res(nPath);	                           //	Allocate results
    cMdl->init(prd.timeline());                        //  Init the simulation timeline
    cRng->init(cMdl->simDim());                        //  Init the RNG
    vector<double> gaussVec(cMdl->simDim());           //  Allocate Gaussian vector
    vector<scenario<T>> path(prd.timeline().size());   //  Allocate path

    //	Iterate through paths	
    bool antiPath = false;
    for (size_t i = 0; i<nPath; i++)
    {
        if (!antithetic)
        {
            //  Next Gaussian vector, dimension D
            cRng->nextG(gaussVec);                        
        }
        else
        {
            //  Antithetic logic
            if (!antiPath)
            {
                cRng->nextG(gaussVec);
                antiPath = true;
            }
            else
            {
                for (auto& gauss : gaussVec) gauss = -gauss;
                antiPath = false;
            }

        }
        cMdl->generatePath(gaussVec, path);       //  Generate path, consume Gaussian vector
        res[i] = prd.payoff(path);              //	Compute result
    }

    return res;	//	C++11: move
}

#include "threadPool.h"

#define BATCHSIZE 64
//	MC simulator: free function that conducts simulations 
//      and returns a vector of nPath payoffs
template <class T>
inline vector<T> mcParallelSimul(
    const Product<T>& prd,
    const Model<T>& mdl,
    const RNG& rng,
    const size_t nPath,
    const bool antithetic)
{
    auto cMdl = mdl.clone();
    auto cRng = rng.clone();

    vector<T> res(nPath);	                           
    cMdl->init(prd.timeline());                        
    cRng->init(cMdl->simDim());                        

    //  Allocate space for Gaussian vectors and paths, 
    //      one for each thread
    ThreadPool *pool = ThreadPool::getInstance();
    const size_t nThread = pool->numThreads();
    vector<vector<double>> gaussVecs(nThread+1);    //  +1 for main
    vector<vector<scenario<T>>> paths(nThread+1);
    for (auto& vec : gaussVecs) vec.resize(cMdl->simDim());
    for (auto& vec : paths) vec.resize(prd.timeline().size());
    
    //  Reserve memory for futures
    vector<TaskHandle> futures;
    futures.reserve(nPath / BATCHSIZE + 1); 

    //  Start
    //  Same as mcSimul() except we send tasks to the pool 
    //  instead of executing them

    size_t firstPath = 0;
    size_t pathsLeft = nPath;
    while (pathsLeft > 0)
    {
        size_t pathsInTask = min(pathsLeft, BATCHSIZE);

        futures.push_back( pool->spawnTask ( [&, firstPath, pathsInTask]()
        {
            //  Inside the parallel task, 
            //      pick the right pre-allocated vectors
            const size_t threadNum = pool->threadNum();
            vector<double>& gaussVec = gaussVecs[threadNum];
            vector<scenario<T>>& path = paths[threadNum];

            //  Get a RNG and position it correctly
            auto taskRng = cRng->clone();
            taskRng->skipTo(antithetic? firstPath / 2: firstPath);

            //  And conduct the simulations, exactly same as sequential
            bool antiPath = false;
            for (size_t i = 0; i < pathsInTask; i++)
            {
                if (!antithetic)
                {
                    //  Next Gaussian vector, dimension D
                    taskRng->nextG(gaussVec);
                }
                else
                {
                    //  Antithetic logic
                    if (!antiPath)
                    {
                        taskRng->nextG(gaussVec);
                        antiPath = true;
                    }
                    else
                    {
                        for (auto& gauss : gaussVec) gauss = -gauss;
                        antiPath = false;
                    }

                }
                cMdl->generatePath(gaussVec, path);       
                res[firstPath + i] = prd.payoff(path);
            }

            //  Remember tasks must return bool
            return true;
        }));

        pathsLeft -= pathsInTask;
        firstPath += pathsInTask;
    }

    //  Wait and help
    for (auto& future : futures) pool->activeWait(future);

    return res;	//	C++11: move
}