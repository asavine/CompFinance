//	ThreadPool.cpp

#include "threadPool.h"

ThreadPool ThreadPool::myInstance;

thread_local size_t ThreadPool::myTLSNum = 0;