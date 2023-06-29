#include "nxbase_core.h"
#include "nxbase_threads.h"
//#include "nxworkerthread.h"

#if !defined(NX_WINDOWS)
#include <pthread.h>
#endif

/*-----------------------------------------------------------------------------
 *					nxWorkerThreadBoostInstance		2009-5-28*/
/** **/
/*---------------------------------------------------------------------------*/

class nxWorkerThreadBoostInstance
{
	private:
		nxWorkerThreadInstance*		m_maininstance;

	public:
						nxWorkerThreadBoostInstance	( nxWorkerThreadInstance* instance)				{ m_maininstance = instance;}
						nxWorkerThreadBoostInstance	( const nxWorkerThreadBoostInstance&  other )	{ m_maininstance = other.m_maininstance;}
		void			operator()					()												{ m_maininstance->ExecuteThread();
																									  m_maininstance->Manager()->FinishWorkerThread( m_maininstance );}
};


/*-----------------------------------------------------------------------------
 *					nxWorkerThreadManager::GetCurrentThreadIdCode		 2014- 10- 22*/
/** **/
/*---------------------------------------------------------------------------*/

size_t nxWorkerThreadManager::GetCurrentThreadIdCode()
{
	size_t	threadid = 0;
#if defined(NX_WINDOWS)
	threadid = GetCurrentThreadId();
#else
	threadid = (size_t)pthread_self();
#endif
	return threadid;
}

/*-----------------------------------------------------------------------------
 *					nxWorkerThreadManager::nxWorkerThreadManager	2009-5-28*/
/** **/
/*---------------------------------------------------------------------------*/

nxWorkerThreadManager::nxWorkerThreadManager()
{
	m_numinstances = 0;												// The number of instances active
}


/*-----------------------------------------------------------------------------
 *					nxWorkerThreadManager::~nxWorkerThreadManager		2009-5-28*/
/** **/
/*---------------------------------------------------------------------------*/

nxWorkerThreadManager::~nxWorkerThreadManager()
{
	if (m_numinstances != 0)
	{
		nxLog::Record(NXLOG_WARNING, "nxWorkerThreadManager::Destructor, We might be destroying a thread manager while objects are still running. Something is not right (a thread crashed or is still running for example).");
	}
};

/*-----------------------------------------------------------------------------
 *					nxWorkerThreadManager::IncrementInstance		2009-5-28*/
/** Safely increments the number of active threads **/
/*---------------------------------------------------------------------------*/

int nxWorkerThreadManager::IncrementInstance()
{
	int	numinstance;
	std::unique_lock<std::mutex>	lock( m_incmutex );				// lock the increment counter 
	//nxSingleLock	lock(&m_incmutex, TRUE);
	numinstance = ++m_numinstances;										// increment the number of instances
	return numinstance;													// and return the number
}


/*-----------------------------------------------------------------------------
 *					nxWorkerThreadManager::DecrementInstance		2009-5-28*/
/** Safely decrements the number of active threads **/
/*---------------------------------------------------------------------------*/

int nxWorkerThreadManager::DecrementInstance()
{
	int	numinstance;

	std::unique_lock<std::mutex>	lock( m_incmutex );			// lock the increment counter 
//	nxSingleLock	lock(&m_incmutex, TRUE);

	numinstance = --m_numinstances;									// and decrement the number of instances
	return numinstance;												
}

/*-----------------------------------------------------------------------------
 *					ThreadEntryPoint		2008-3-6*/
/** Fires up the thread controller threads.  Starts execution at
 *	each thread's BeginThread member.
 **/
/*---------------------------------------------------------------------------*/

//#if !defined(NX_USEBOOST_THREADLIBS)
//static void ThreadEntryPoint( void * address )
//{
//	nxWorkerThreadInstance*			userthread = (nxWorkerThreadInstance *)address;
//	nxWorkerThreadBoostInstance		instance(userthread);
//	instance();
//}
//#endif



/*-----------------------------------------------------------------------------
 *					nxWorkerThreadManager::StartWorkerThread		2009-5-28*/
/** This starts a worker thread. The paradigm assumes we are running in the
 *	context of the main processing thread, which is lauching a whole bunch of
 *	worker threads (virtualized in class nxWorkerThreadInstance). Typically
 *	once all the threads are launched the main thread will call
 *	WaitForThreadsToFinish ().
 *
 *	\par Notes
 *	We have specifically classed this stuff up so a small special class is actually used
 *	to launch the thread. I like to keep stats about the executing thread available to the
 *	main program. This does not work if the nxWorkerThreadInstance is the instance actually launched by std::thread. 
 *	Boost::Thread makes a new copy instance of the workerthread and passes the new instance to the thread. 
 *	The solution is to use a a small class nxWorkerThreadBoostInstance that contains pointers to our thread instance,
 *	launch the thread and call method ExecuteThread. This also allows our small class to inject thread management code
 *	both before and after the main thread execution. Makes life a little neater and tidier.
 **/
/*---------------------------------------------------------------------------*/

bool nxWorkerThreadManager::StartWorkerThread( nxWorkerThreadInstance* userthread )
{
//	int nthreads;

	IncrementInstance();								// Increment the number of threads running
	userthread->SetManager(this);								// Set the thread manager in the new thread to this object
	userthread->SetIsRunning(true);								// Flag this object that it is running

	std::thread( nxWorkerThreadBoostInstance( userthread ) );	// Launch the actual thread.
//	if (nthreads == 1) m_threadactivecond.ResetEvent();				// reset the thread active condition event
//	_beginthread(  &ThreadEntryPoint, 0, (void *)userthread );		// Now start the thread, clears busy when done 
	return true;	
}


/*-----------------------------------------------------------------------------
 *					nxWorkerThreadManager::FinishWorkerThread		2009-5-28*/
/** This is called by the worker thread. It is exposed as a public member but is
 *	really a private member and should not be called by most users. It executes
 *	within the thread context. It is called in the dying stages of the thread.
 *	Special helper class nxWorkerBoostInstance calls this method after calling
 *	nxWorkerThreadInstance::ExecuteThread and just before the thread terminates
 **/
/*---------------------------------------------------------------------------*/

void nxWorkerThreadManager::FinishWorkerThread( nxWorkerThreadInstance* userthread )
{
	int		n;
	userthread->SetIsRunning(false);							// Flag this object that it is running
	userthread->SetManager(NULL);								// Set the thread manager in the new thread to this object
	n = DecrementInstance();									// Increment the number of threads running
	if (n <= 0)													// If we have no more threads left
	{															// then
		m_threadactivecond.notify_one();						// notify the main thread
//		m_threadactivecond.SetEvent();							// notify the main thread by setting the thread active condition
	}
}

/*-----------------------------------------------------------------------------
 *					nxWorkerThreadManager::WaitForThreadsToFinish		2009-5-28*/
/** This is designed to run in the context of the main thread. It waits for all
 *	of the worker threads to finish. Note it may hit a deadlock if the worker
 *	threads do not terminate properly and unregsiter themselves from the manager.
 *	It should not be called by a worker thread as it will cause a deadlock as
 *	the worker thread will stay alive until it dies!
 **/
/*---------------------------------------------------------------------------*/

bool nxWorkerThreadManager::WaitForThreadsToFinish()
{
	std::unique_lock<std::mutex>	lock(m_incmutex);

	while(m_numinstances > 0)
    {
        m_threadactivecond.wait(lock);
    }
//	WaitForSingleObject( m_threadactivecond, INFINITE);
	return true;
}
