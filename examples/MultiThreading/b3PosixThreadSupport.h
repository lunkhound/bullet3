/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2007 Erwin Coumans  http://bulletphysics.com

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef B3_POSIX_THREAD_SUPPORT_H
#define B3_POSIX_THREAD_SUPPORT_H


#include "Bullet3Common/b3Scalar.h"


#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 //for definition of pthread_barrier_t, see http://pages.cs.wisc.edu/~travitch/pthreads_primer.html
#endif //_XOPEN_SOURCE
#include <pthread.h>
#include <semaphore.h>



#include "Bullet3Common/b3AlignedObjectArray.h"

#include "b3ThreadSupportInterface.h"


// b3PosixThreadSupport helps to initialize/shutdown libspe2, start/stop SPU tasks and communication
class b3PosixThreadSupport : public b3ThreadSupportInterface
{
public:
    typedef enum sStatus {
        STATUS_BUSY,
        STATUS_READY,
        STATUS_FINISHED
    } Status;

	// placeholder, until libspe2 support is there
	struct	b3ThreadStatus
	{
		int	m_taskId;
		int	m_commandId;
		int	m_status;

		ThreadFunc	m_userThreadFunc;
		void*	m_userPtr; //for taskDesc etc
		void*	m_lsMemory; //initialized using PosixLocalStoreMemorySetupFunc

        pthread_t thread;
        //each tread will wait until this signal to start its work
        sem_t* startSemaphore;

        // this is a copy of m_mainSemaphore, 
        //each tread will signal once it is finished with its work
        sem_t* m_mainSemaphore;
        unsigned long threadUsed;
	};
private:

	b3AlignedObjectArray<b3ThreadStatus>	m_activeThreadStatus;

	// m_mainSemaphoresemaphore will signal, if and how many threads are finished with their work
	sem_t* m_mainSemaphore;
	
public:
	///Setup and initialize SPU/CELL/Libspe2

    typedef ConstructionInfo ThreadConstructionInfo;

	b3PosixThreadSupport(const ConstructionInfo& threadConstructionInfo);

///cleanup/shutdown Libspe2
	virtual	~b3PosixThreadSupport();

	void	startThreads(const ConstructionInfo& threadInfo);


	virtual	void runTask(int uiCommand, void* uiArgument0, int uiArgument1);

	virtual	void waitForResponse(int *puiArgument0, int *puiArgument1);


///tell the task scheduler we are done with the SPU tasks
	virtual	void stopThreads();

	virtual void setNumTasks(int numTasks) {}

	virtual int getNumTasks() const
	{
		return m_activeThreadStatus.size();
	}

	///non-blocking test if a task is completed. First implement all versions, and then enable this API
	virtual bool isTaskCompleted(int *puiArgument0, int *puiArgument1, int timeOutInMilliseconds);


	virtual b3Barrier* createBarrier();

	virtual b3CriticalSection* createCriticalSection();

	virtual void deleteBarrier(b3Barrier* barrier);

	virtual void deleteCriticalSection(b3CriticalSection* criticalSection);


	virtual void*	getThreadLocalMemory(int taskId)
	{
		return m_activeThreadStatus[taskId].m_lsMemory;
	}

};


#endif // B3_POSIX_THREAD_SUPPORT_H


