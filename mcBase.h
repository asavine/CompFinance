#pragma once

#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>

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

    //  Access to all parameters by copy
    virtual vector<T> parameters() const = 0;

    //  For AAD enabled models
    //  Put parameters on tape 
    virtual void putOnTape() {}
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
inline vector<double> mcSimul(
    const Product<double>& prd,
    const Model<double>& mdl,
    const RNG& rng,			            
    const size_t nPath,
    const bool antithetic)                      
{
    //  Work with copies of the model and RNG
    //      which are modified when we set up the simulation
    //  Copies are OK at high level
    auto cMdl = mdl.clone();
    auto cRng = rng.clone();

    //	Allocate results
    vector<double> res(nPath);
    //  Init the simulation timeline
    cMdl->init(prd.timeline());              
    //  Init the RNG
    cRng->init(cMdl->simDim());                        
    //  Allocate Gaussian vector
    vector<double> gaussVec(cMdl->simDim());           
    //  Allocate path
    vector<scenario<double>> path(prd.timeline().size());   

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
        //  Generate path, consume Gaussian vector
        cMdl->generatePath(gaussVec, path);     
        //	Compute result
        res[i] = prd.payoff(path);              
    }

    return res;	//	C++11: move
}

#include "threadPool.h"

#define BATCHSIZE 64
//	MC simulator: free function that conducts simulations 
//      and returns a vector of nPath payoffs
inline vector<double> mcParallelSimul(
    const Product<double>& prd,
    const Model<double>& mdl,
    const RNG& rng,
    const size_t nPath,
    const bool antithetic)
{
    auto cMdl = mdl.clone();
    auto cRng = rng.clone();

    vector<double> res(nPath);	                           
    cMdl->init(prd.timeline());                        
    cRng->init(cMdl->simDim());                        

    //  Allocate space for Gaussian vectors and paths, 
    //      one for each thread
    ThreadPool *pool = ThreadPool::getInstance();
    const size_t nThread = pool->numThreads();
    vector<vector<double>> gaussVecs(nThread+1);    //  +1 for main
    vector<vector<scenario<double>>> paths(nThread+1);
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
        size_t pathsInTask = min<size_t>(pathsLeft, BATCHSIZE);

        futures.push_back( pool->spawnTask ( [&, firstPath, pathsInTask]()
        {
            //  Inside the parallel task, 
            //      pick the right pre-allocated vectors
            const size_t threadNum = pool->threadNum();
            vector<double>& gaussVec = gaussVecs[threadNum];
            vector<scenario<double>>& path = paths[threadNum];

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

#include "AAD.h"

//	MC simulator: free function that conducts simulations 
//  Note we return:
//  - a vector of pathwise payoffs as usual
//  - a clone of the original model, 
//      with derivatives accumulated inside the adjoints of the model parameters 
struct AADSimulResults
{
    vector<double> payoffs;
    unique_ptr<Model<Number>> model;
};
//  Also note: tape must be wiped afterwards
inline AADSimulResults
mcSimulAAD(
    const Product<Number>& prd,
    const Model<Number>& mdl,
    const RNG& rng,
    const size_t nPath,
    const bool antithetic)
{
    //  Work with copies of the model and RNG
    //      which are modified when we set up the simulation
    //  Copies are OK at high level
    auto cMdl = mdl.clone();
    auto cRng = rng.clone();

    //  AAD - 1
    //  Access to tape
    Tape& tape = *Number::tape;
    //  Rewind tape
    tape.rewind();
    //  Put parameters on tape
    //  note that also initializes all adjoints
    cMdl->putOnTape();
    //  Init the simulation timeline
    //  CAREFUL: simulation timeline must be on tape
    //  Hence moved here
    cMdl->init(prd.timeline());                        
    //  Mark the tape straight after initialization
    tape.mark();
    //

    //  Other init / allocate
    vector<Number> nPayoffs(nPath);	                    //  Allocate results
    cRng->init(cMdl->simDim());                         //  Init the RNG
    vector<double> gaussVec(cMdl->simDim());            //  Allocate Gaussian vector
    vector<scenario<Number>> path(prd.timeline().size());   //  Allocate path

    //	Iterate through paths	
    bool antiPath = false;
    for (size_t i = 0; i<nPath; i++)
    {
        //  AAD - 2
        //  Rewind tape to mark
        //  parameters stay on tape but the rest is wiped
        tape.rewindToMark();
        //

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
        //  Generate path, consume Gaussian vector
        cMdl->generatePath(gaussVec, path);     
        //	Compute result
        nPayoffs[i] = prd.payoff(path);              

        //  AAD - 3
        //  Propagate adjoints
        nPayoffs[i].propagateToMark();
        //
    }

    //  AAD - 4
    //  Mark = limit between pre-calculations and path-wise operations
    //  Operations above mark have been propagated and accumulated
    //  We conduct one propagation mark to start
    Number::propagateMarkToStart();

    //  Results
    AADSimulResults results;
    
    //  Pathwise payoffs
    results.payoffs.resize(nPath);
    //  Copy and convert
    convertCollection(nPayoffs.begin(), nPayoffs.end(), results.payoffs.begin());
    
    //  Model with accumulated adjoints
    results.model = move(cMdl);

    return results;
}

//  Extract initialization, except RNG
template<class T>
void initSimul(
    //  Inputs
    const Product<T>& prd,
    const Model<T>& mdl,
    //  Stuff to init
    unique_ptr<Model<T>>& mdlClone,
    //  Stuff to allocate
    vector<double>& gaussVec,
    vector<scenario<T>>& path)
{
    //  Work with copies of the model and RNG
    //      which are modified when we set up the simulation
    //  Copies are OK at high level
    mdlClone = mdl.clone();

    //  Access to tape
    Tape& tape = *Number::tape;
    //  Rewind tape
    tape.rewind();
    //  Put parameters on tape
    //  note that also initializes all adjoints
    mdlClone->putOnTape();
    //  Init the simulation timeline
    //  CAREFUL: simulation timeline must be on tape
    //  Hence moved here
    mdlClone->init(prd.timeline());
    //  Mark the tape straight after parameters
    tape.mark();
    //

    //  Other init / allocate
    //  Allocate Gaussian vector
    gaussVec.resize(mdlClone->simDim());     
    //  Allocate path
    path.resize(prd.timeline().size());   

}

//	MC simulator: free function that conducts simulations 
//  Note we return:
//  - a vector of pathwise payoffs as usual
//  - a clone of the original model, 
//      with derivatives accumulated inside the adjoints of the model parameters 
//  Also note: tape must be wiped afterwards
inline AADSimulResults
mcParallelSimulAAD(
    const Product<Number>& prd,
    const Model<Number>& mdl,
    const RNG& rng,
    const size_t nPath,
    const bool antithetic)
{
    //  Allocate results
    vector<Number> nPayoffs(nPath);	                       

    //  We need one of all these for each thread
    //  0: main thread
    //  1 to n : worker threads

    ThreadPool *pool = ThreadPool::getInstance();
    const size_t nThread = pool->numThreads();

    //  Clones
    vector<unique_ptr<Model<Number>>> cMdl(nThread + 1);

    //  Space for Gaussian vectors and paths 
    vector<vector<double>> gaussVecs(nThread + 1);   
    vector<vector<scenario<Number>>> paths(nThread + 1);

    //  Tapes for the worker threads
    //  The main thread has one of its own
    vector<Tape> tapes(nThread);

    //  Initialized indicators
    vector<int> init(nThread + 1, false);

    //  Initialize main thread
    initSimul(prd, mdl, cMdl[0], gaussVecs[0], paths[0]);

    //  Init the RNG
    auto cRng = rng.clone();
    cRng->init(cMdl[0]->simDim());

    //  Mark as initialized
    init[0] = true;

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
        size_t pathsInTask = min<size_t>(pathsLeft, BATCHSIZE);

        futures.push_back(pool->spawnTask([&, firstPath, pathsInTask]()
        {
            const size_t threadNum = pool->threadNum();

            //  Use this thread's tape
            //  Thread local magic: each thread its own pointer
            //  Note main thread = 0 is not reset
            if (threadNum > 0) Number::tape = &tapes[threadNum - 1];

            //  Initialize once on each thread
            if (!init[threadNum])
            {
                //  Initialize
                initSimul(prd, mdl, 
                    cMdl[threadNum], 
                    gaussVecs[threadNum], 
                    paths[threadNum]);

                //  Mark as initialized
                init[threadNum] = true;
            }

            //  Get a RNG and position it correctly
            auto taskRng = cRng->clone();
            taskRng->skipTo(antithetic ? firstPath / 2 : firstPath);

            //  And conduct the simulations, exactly same as sequential
            bool antiPath = false;
            for (size_t i = 0; i < pathsInTask; i++)
            {
                //  Rewind tape to mark
                Number::tape->rewindToMark();

                if (!antithetic)
                {
                    //  Next Gaussian vector, dimension D
                    taskRng->nextG(gaussVecs[threadNum]);
                }
                else
                {
                    //  Antithetic logic
                    if (!antiPath)
                    {
                        taskRng->nextG(gaussVecs[threadNum]);
                        antiPath = true;
                    }
                    else
                    {
                        for (auto& gauss : gaussVecs[threadNum]) 
                            gauss = -gauss;
                        antiPath = false;
                    }

                }
                cMdl[threadNum]->generatePath(
                    gaussVecs[threadNum], 
                    paths[threadNum]);
                nPayoffs[firstPath + i] = prd.payoff(paths[threadNum]);

                //  Propagate adjoints
                nPayoffs[firstPath + i].propagateToMark();
            }

            //  Remember tasks must return bool
            return true;
        }));

        pathsLeft -= pathsInTask;
        firstPath += pathsInTask;
    }

    //  Wait and help
    for (auto& future : futures) pool->activeWait(future);
    
    //  Mark = limit between pre-calculations and path-wise operations
    //  Operations above mark have been propagated and accumulated
    //  We conduct one propagation mark to start
    //  On the main thread's tape
    Number::propagateMarkToStart();
    //  And on the worker thread's tapes
    Tape* mainThreadPtr = Number::tape;
    for (size_t i = 0; i < nThread; ++i)
    {
        if (init[i+1])
        {
            //  Set tape pointer
            Number::tape = &tapes[i];
            //  On that tape, propagate
            Number::propagateMarkToStart();
        }
    }
    //  Reset tape to main thread's
    Number::tape = mainThreadPtr;

    //  At this point, we have nThread + 1 models
    //  Each model i has parameters that accumlated adjoints
    //      on its own tape for all paths on thread i
    //  We accumulate them all on the main thread's model
    vector<Number> params0 = cMdl[0]->parameters();
    for (size_t i = 0; i < nThread; ++i)
    {
        if (init[i+1])
        {
            vector<Number> params = cMdl[i+1]->parameters();
            for (size_t j = 0; j < params0.size(); ++j)
            {
                params0[j].adjoint() += params[j].adjoint();
            }
        }
    }

    //  Results
    AADSimulResults results;

    //  Pathwise payoffs
    results.payoffs.resize(nPath);
    //  Copy and convert
    convertCollection(nPayoffs.begin(), nPayoffs.end(), results.payoffs.begin());

    //  Model with accumulated adjoints
    results.model = move(cMdl[0]);

    return results;
}



