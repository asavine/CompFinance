#pragma once

//	ThreadPool.h

#include <future>
#include <thread>
#include "ConcurrentQueue.h"

using namespace std;

typedef packaged_task<bool(void)> Task;
typedef future<bool> TaskHandle;

class ThreadPool 
{
	//	The one and only instance
	static ThreadPool myInstance;

	//	The task queue
	ConcurrentQueue<Task> myQueue;

	//	The threads
	vector<thread> myThreads;

	//	Interruption indicator
	bool myInterrupt;

	//	Thread number
	static thread_local size_t myTLSNum;

	//	The function that is executed on every thread
	void threadFunc(const size_t num)
	{
		myTLSNum = num;

		Task t;

		//	"Infinite" loop, only broken on destruction
		while (!myInterrupt) 
		{
			//	Pop and executes tasks
			myQueue.pop(t);
			t();			
		}
	}

public:

	//	Access the instance
	static ThreadPool* getInstance() { return &myInstance; }

	//	Number of threads
	size_t numThreads() const { return myThreads.size(); }

	//	The number of the caller thread
	static size_t threadNum() { return myTLSNum; }

	//	Constructor
	ThreadPool(const size_t nThread = thread::hardware_concurrency() - 1)
		: myInterrupt(false)
	{
		myThreads.reserve(nThread);
		
		//	Launch threads on threadFunc and keep handles in a vector
		for (size_t i = 0; i<nThread; i++)
			myThreads.push_back(thread(&ThreadPool::threadFunc, this, i+1));
	}

	//	Destructor
	~ThreadPool()
	{
		//	Interrupt mode
		myInterrupt = true;

		//	Interrupt all waiting threads
		myQueue.interrupt();

		//	Wait for them all to join
		for_each(myThreads.begin(), myThreads.end(), mem_fn(&thread::join));
	}

	//	Forbid copies etc
	ThreadPool(const ThreadPool& rhs) = delete;
	ThreadPool& operator=(const ThreadPool& rhs) = delete;
	ThreadPool(ThreadPool&& rhs) = delete;
	ThreadPool& operator=(ThreadPool&& rhs) = delete;

	//	Spawn task
	template<typename Callable>
	TaskHandle spawnTask(Callable c)
	{
		Task t(move(c));
		TaskHandle f = t.get_future();
		myQueue.push(move(t));
		return f;
	}

	//	Run queued tasks synchronously 
	//	while waiting on a future, 
	//	return true if at least one task was run
	bool activeWait(const TaskHandle& f)
	{
		Task t;
		bool b = false;

		//	Check if the future is ready without blocking
		//	The only syntax C++11 provides for that is
		//	wait 0 seconds and return status
		while (f.wait_for(0s) != future_status::ready)
		{
			//	Non blocking
			if (myQueue.tryPop(t)) 
			{
				t();
				b = true;
			}
			else //	Nothing in the queue: go to sleep
			{
				f.wait();
			}
		}

		return b;
	}
};