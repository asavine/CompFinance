
/*
Written by Antoine Savine in 2018

This code is the strict IP of Antoine Savine

License to use and alter this code for personal and commercial applications
is freely granted to any person of company who purchased a copy of the book

Modern Computational Finance: AAD and Parallel Simulations
Antoine Savine
Wiley, 2018

As long as this comment is preserved at the top of the file
*/

#pragma once

#include <queue>
#include <mutex>

using namespace std;

template <class T>
class ConcurrentQueue
{

	queue<T> myQueue;
	mutable mutex myMutex;
	condition_variable myCV;
	bool myInterrupt;

public:

	ConcurrentQueue() : myInterrupt(false) {}
	~ConcurrentQueue() { interrupt(); }

	bool empty() const
	{
		//	Lock
		lock_guard<mutex> lk(myMutex);
		//	Access underlying queue
		return myQueue.empty();
	}	//	Unlock

		//	Pop into argument
	bool tryPop(T& t)
	{
		//	Lock
		lock_guard<mutex> lk(myMutex);
		if (myQueue.empty()) return false;
		//	Move from queue
		t = move(myQueue.front());
		//	Combine front/pop
		myQueue.pop();

		return true;
	}	//	Unlock

		//	Pass t byVal or move with push( move( t))
	void push(T t)
	{
		{
			//	Lock
			lock_guard<mutex> lk(myMutex);
			//	Move into queue
			myQueue.push(move(t));
		}	//	Unlock before notification

			//	Unlock before notification 
		myCV.notify_one();
	}

	//	Wait if empty
	bool pop(T& t)
	{
		//	(Unique) lock
		unique_lock<mutex> lk(myMutex);

		//	Wait if empty, release lock until notified 
		while (!myInterrupt && myQueue.empty()) myCV.wait(lk);

		//	Re-acquire lock, resume 

		//  Check for interruption
		if (myInterrupt) return false;

		//	Combine front/pop 
		t = move(myQueue.front());
		myQueue.pop();

		return true;

	}	//	Unlock

	void interrupt()
	{
        {
            lock_guard<mutex> lk(myMutex);
            myInterrupt = true;
        }
		myCV.notify_all();
	}
};