#pragma once

# pragma once

#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

#include <algorithm>
#include <list>
#include <queue>
#include <functional>
#include <memory>
#include <future>
#include <string>
#include <sstream>

using namespace std;

//	Semaphore

template <class LOCKABLE>
class Semaphore 
{

	LOCKABLE myMutex;
	condition_variable myCV;

	int myN;

private:

	Semaphore(const Semaphore& rhs);
	Semaphore& operator=(const Semaphore& rhs);

public:

	Semaphore(int n) : myN(n) {}

	void post(int p = 1) 
	{ 
		lock_guard<LOCKABLE> l(myMutex);
		myN += p; 
		for (int i = 0; i<p; ++i) myCV.notify_one(); 
	}

	void wait() 
	{ 
		unique_lock<LOCKABLE> l(myMutex);
		
		while (!myN) 
		{ 
			myCV.wait(l); 
		} 
		
		myN--; 
	}

	//	Lockable interface
	void lock() { wait(); }
	void unlock() { post(); }
};

//	Tests

class semaphore_tester 
{

protected:

	mutex coutm;
	Semaphore<mutex> sem;

	void thread_func(int sleep_milli);

public:

	semaphore_tester(int n) : sem(n) {}
};

//	Initialize to locked then unlock
class semaphore_tester1 : public semaphore_tester 
{

public:

	semaphore_tester1() : semaphore_tester(0) {}			//	Initialized to locked
	void go1(int nthreads);									//	Do some work, then unlock for 1 thread
	void go2(int nthreads, int m);							//	Do some wotk, then unlock for m threads
};

//	Let n threads through
class semaphore_tester2 : public semaphore_tester 
{

public:

	semaphore_tester2(int n) : semaphore_tester(n) {}		//	Initialized to open for n threads
	void go(int nthreads);									//	Go with nthreads
};

//	Barrier

template <class LOCKABLE>
class Barrier 
{
	int myN0;													//	Initial capacity
	int myN;													//	Current capacity = n0 - threads_waiting

	int myStep;													//	or generation, starts at 0, incremented when all threads arrive

	Barrier(const Barrier& rhs);
	Barrier& operator=(const Barrier& rhs);

	LOCKABLE myMutex;
	condition_variable myCV;

public:

	Barrier(const int n) : myN0(n), myN(n), myStep(0) {}

	void wait() 
	{
		unique_lock<LOCKABLE> l(myMutex);

		int old_step = myStep;
		--myN;

		if (!myN) 
		{
			myN = myN0;
			++myStep;
			myCV.notify_all();
		}
		else 
		{
			while (myStep == old_step) myCV.wait(l);
		}
	}

	//	Reduce the barrier by 1
	void opt_out() 
	{
		lock_guard<LOCKABLE> l(myMutex);

		--myN0;
		--myN;

		if (!myN) 
		{

			myN = myN0;
			++myStep;
			myCV.notify_all();
		}
	}
};

class barrier_tester 
{

	int myNthread;
	Barrier<mutex> myBarrier;
	void thread_func(const int waitms);

public:
	barrier_tester(const int nthread) : myNthread(nthread), myBarrier(nthread) {}
	void go();
};

//	Shared mutex and lock

class shared_mutex 
{

	mutex myMutex;													//	Internal mutex: only for protection of the shared mutex
																	//	Actual locks work with the condition variables

	condition_variable myCVshared, myCVexcl;						//	Condition variables for the exclusive and shared locks

	int myActShared, myActExcl;										//	Number of exclusive and shared locks currently locked
	int myWaitShared, myWaitExcl;									//	Number of exclusive and shared locks waiting

																	//	Number of shared locks allowed through since at least one exclusive lock is waiting
	int myNsharedAllowed;											//		Default: hardware concurrency
																	//		Reset to 0 when an exclusive lock is granted
																	//	Number of esclusive locks allowed through since at least one shared lock is waiting
	int myNexclAllowed;												//		Default: 1
																	//		Reset to 0 when a shared lock is granted

																	//	Maximum shared/exclusive locks to be allowed while others are waiting
																	//		Once maximum is reached, the other type of lock takes priority and nshared_allowed/nexcl_allowed is reset to 0 accordingly
	int myMaxNsharedAllowed, myMaxNexclAllowed;

public:

	shared_mutex(int max_nshared_allowed = thread::hardware_concurrency(), int max_nexcl_allowed = 1)
		: myMaxNsharedAllowed(max_nshared_allowed > 0 ? max_nshared_allowed : 2), myMaxNexclAllowed(max_nexcl_allowed),
		myNsharedAllowed(0), myNexclAllowed(0),
		myActShared(0), myActExcl(0), myWaitShared(0), myWaitExcl(0) {}

	void lock();
	void lock_shared();

	void unlock();
	void unlock_shared();
};

//	Shared lock template, to be used like unique_lock but with shared mutexes in order to obtain a shared lock

template <class M>
class shared_lock 
{

	M* mySM;

public:

	shared_lock() { m = nullptr };
	shared_lock(M& m) { mySM = &m; mySM->lock_shared(); }
	shared_lock(M& m, defer_lock_t) : { mySM = &m; }

	void lock() { if (mySM) mySM->lock_shared(); }
	void unlock() { if (mySM) mySM->unlock_shared(); }

	~shared_lock() { unlock(); }
};

//	Testing

class shared_lock_driver 
{

	//	The shared mutex
	shared_mutex sm;
	//	This is for synchronization of cout
	mutex cout_m;

	list<thread> threads;

	void launch_writer_thread(int len);
	void launch_reader_thread(int len);

public:

	void go();

	//	Join all on destruction
	~shared_lock_driver() { for_each(threads.begin(), threads.end(), mem_fn(&thread::join)); cout << " All threads clear " << endl; }
};






