
#include "SyncConstructs.h"

#include <cstdlib>

#include <iostream>
#include <chrono>

using namespace std;

//	Semaphore

void semaphore_tester::thread_func(int sleep_milli) 
{
	lock_guard<Semaphore<mutex>> l(sem);
	coutm.lock(); cout << ".Thread " << this_thread::get_id() << " entering" << endl; coutm.unlock();
	this_thread::sleep_for(chrono::milliseconds(sleep_milli));
	coutm.lock(); cout << "...Thread " << this_thread::get_id() << " exiting" << endl; coutm.unlock();
}

void semaphore_tester1::go1(int nthreads) 
{
	vector<thread> threads;
	for (int i = 0; i<nthreads; i++) threads.push_back(thread(bind(&semaphore_tester1::thread_func, this, 1)));

	coutm.lock(); cout << " Main thread initializing work" << endl; coutm.unlock();
	this_thread::sleep_for(chrono::milliseconds(2000));
	coutm.lock(); cout << " Main thread initializing work finished - posting" << endl; coutm.unlock();
	sem.post();

	for (int i = 0; i<nthreads; i++) threads[i].join();
	cout << " All done ";
}

void semaphore_tester1::go2(int nthreads, int m) 
{
	vector<thread> threads;
	for (int i = 0; i<nthreads; i++) threads.push_back(thread(bind(&semaphore_tester1::thread_func, this, 250)));

	coutm.lock(); cout << " Main thread initializing work" << endl; coutm.unlock();
	this_thread::sleep_for(chrono::milliseconds(1500));
	coutm.lock(); cout << " Main thread initializing work finished - posting" << endl; coutm.unlock();
	sem.post(m);

	for (int i = 0; i<nthreads; i++) threads[i].join();
	cout << " All done ";
}

void semaphore_tester2::go(int nthreads) 
{
	srand(unsigned(time(NULL)));
	vector<thread> threads;
	for (int i = 0; i<nthreads; i++) threads.push_back(thread(bind(&semaphore_tester2::thread_func, this, rand() % 50)));

	for (int i = 0; i<nthreads; i++) threads[i].join();
	cout << " All done ";
}

//	Barrier

void barrier_tester::thread_func(const int waitms) 
{
	cout << "thread id " << this_thread::get_id() << " starting" << endl;
	cout << "thread id " << this_thread::get_id() << " sleeping for " << waitms << " milliseconds" << endl;

	this_thread::sleep_for(chrono::milliseconds(waitms));

	cout << "thread id " << this_thread::get_id() << " at barrier" << endl;
	myBarrier.wait();
	cout << "thread id " << this_thread::get_id() << " past barrier" << endl;
}

void barrier_tester::go() {

	vector<thread> vt;
	srand(unsigned(time(NULL)));

	for (int i = 0; i<5 * myNthread; i++) vt.push_back(thread(&barrier_tester::thread_func, this, rand() % 1000));

	cout << "joining" << endl;
	for (int i = 0; i<5 * myNthread; i++) vt[i].join();
	cout << "joined" << endl;
}

//	Shared mutex and lock

//	Exclusive lock
void shared_mutex::lock() 
{
	unique_lock<mutex> l(myMutex);

	//	Wait if either exclusive or shared lock is already granted, or if the maximum number of exclusive locks has been granted while shared locks are waiting
	if (myActShared > 0 || myActExcl > 0 || myNexclAllowed >= myMaxNexclAllowed) 
	{
		myWaitExcl++;
		while (myActShared > 0 || myActExcl >0 || myNexclAllowed >= myMaxNexclAllowed) myCVexcl.wait(l);
		myWaitExcl--;
	}

	//	Lock is granted

	//	Mark exclusive lock as granted
	myActExcl++;

	//	Reset the number of shared locks granted while exclusive locks are waiting
	myNsharedAllowed = 0;

	//	Increment number of exclusive locks granted if shared locks are waiting
	if (myWaitShared > 0) myNexclAllowed++;
}

//	Shared lock
void shared_mutex::lock_shared() 
{
	unique_lock<mutex> l(myMutex);

	//	Wait if exclusive lock is already granted, or if the maximum number of shared locks has been granted while exclusive locks are waiting
	if (myActExcl > 0 || myNsharedAllowed >= myMaxNsharedAllowed) {
		myWaitShared++;
		while (myActExcl > 0 || myNsharedAllowed >= myMaxNsharedAllowed) myCVshared.wait(l);
		myWaitShared--;
	}

	//	Lock is granted

	//	Increase the number of shared locks currently granted
	myActShared++;

	//	Reset the number of exclusive locks granted while shared locks are waiting
	myNexclAllowed = 0;

	//	Increment number of shared locks granted if exclusive locks are waiting
	if (myWaitExcl > 0) myNsharedAllowed++;
}

//	Unlocks

void shared_mutex::unlock() {

	unique_lock<mutex> l(myMutex);

	//	Decrement exclusive locks
	myActExcl--;

	//	If other exclusive locks are waiting AND the maximum allowed while shared locks are waiting has not been achieved, then give priority to exclusive locks
	if (myWaitExcl > 0 && myNexclAllowed < myMaxNexclAllowed) 
	{
		myCVexcl.notify_one();
	}
	//	Otherwise grant shared locks (if any waiting)
	else if (myWaitShared > 0) 
	{
		myCVshared.notify_all();
	}
}

void shared_mutex::unlock_shared() 
{
	unique_lock<mutex> l(myMutex);

	//	Decrement shared locks
	myActShared--;

	//	If no other shared locks are currently granted, and exclusive locks are waiting, then wake one
	if (myActShared == 0 && myWaitExcl > 0) 
	{
		myCVexcl.notify_one();
	}

	//	Otherwise no need to explicitly wake a waiting shared lock, wince these can only wait on waiting or active exclusives
}

void shared_lock_driver::launch_writer_thread(int len) 
{
	unique_lock<shared_mutex> lock(sm);
	{ lock_guard<mutex> lg(cout_m); cout << "	Writer thread " << this_thread::get_id() << " starting" << endl; }
	this_thread::sleep_for(chrono::milliseconds(len));
	{ lock_guard<mutex> lg(cout_m); cout << "	<-Writer thread " << this_thread::get_id() << " exiting" << endl; }
}

void shared_lock_driver::launch_reader_thread(int len) 
{
	shared_lock<shared_mutex> lock(sm);
	{ lock_guard<mutex> lg(cout_m); cout << "Reader thread " << this_thread::get_id() << " starting" << endl; }
	this_thread::sleep_for(chrono::milliseconds(len));
	{ lock_guard<mutex> lg(cout_m); cout << "<-Reader thread " << this_thread::get_id() << " exiting" << endl; }
}

void shared_lock_driver::go() 
{
	srand(98765);
	for (int i = 0; i<50; i++) {
		this_thread::sleep_for(chrono::milliseconds(rand() % 50));
		int rn = rand() % 100;
		if (i>0 && rn<30) {
			lock_guard<mutex> l(cout_m);
			cout << "....Launching writer thread " << i << " : ";
			threads.push_back(thread(&shared_lock_driver::launch_writer_thread, this, rand() % 100));
			cout << threads.back().get_id() << endl;
		}
		else {
			lock_guard<mutex> l(cout_m);
			cout << "....Launching reader thread " << i << " : ";
			threads.push_back(thread(&shared_lock_driver::launch_reader_thread, this, rand() % 100));
			cout << threads.back().get_id() << endl;
		}
	}
}



