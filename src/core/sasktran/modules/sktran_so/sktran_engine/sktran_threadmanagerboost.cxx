#include "../sasktranv21_internals.h"
#include <float.h>
#include <nxbase_threads.h>
#include <mutex>


class SKTRAN_ThreadInstanceBoost; 

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerBoost		2011-5-31*/
/** The Boos multi-threaded implementation**/
/*---------------------------------------------------------------------------*/

class SKTRAN_ThreadManagerBoost : public SKTRAN_ThreadManager
{


	private:
		SKTRAN_ThreadInstanceBoost*						m_threads;
		size_t										m_numthreads;
		static size_t								NumDefaultProcessingThreads	();


	private:
		size_t										m_nextpoint;				// The index of the next point or ray that will be processed
		std::mutex								m_nextpointmutex;			// A mutex used when child threads call the parent to get the next point to process
		std::condition_variable					m_threadsync;				// A condition variable used to make all of the child threads wait for a notification

	private:
		void										ReleaseResources					();
		bool										ResetPointIndex						( size_t numpoints );


	public:
													SKTRAN_ThreadManagerBoost	();
		virtual									   ~SKTRAN_ThreadManagerBoost	();
		std::condition_variable*					ConditionalVariable			()	{return &m_threadsync;}				//!< Used by child worker threads to wait for requests from manager
		bool										GetNextPointIndex			( size_t* pointindex );					//!< Returns the index of the next point or ray to be processed

	public:
		virtual bool								SetOpticalProps						( const SKTRAN_TableOpticalProperties_V21* optprop);
		virtual bool								CreateThreads						( size_t numthreads, const SKTRAN_SpecsInternal_V21* modelspecifications, const SKTRAN_EngineDiffuseTables* modeltables );
		virtual bool								NotifyWorkerThreadsAndWaitForEnd	( size_t numpoints );
		virtual bool								CloseThreadsAndRelease				( );
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_ThreadInstanceBoost						2010-3-23*/
/** A class for worker threads. The worker threads are "owned" by a
 *	SKTRAN_ThreadManager class. The manager class runs in the main program
 *	thread while the workers run in their own thread. The worker threads
 *	are created once and are re-used throughout the processing. This helps
 *	eliminate inefficiences associated with thread creation etc. I have been
 *	quite careful with this code so please dont just change it unless
 *	you are really (really really) sure there is something wrong.
 *
 *	Each thread once initialized and launched is either BUSY processing
 *	computational requests or is WAITING for new requests. When the thread 
 *	is WAITING the thread's internal mutex is unlocked by the call to 
 *	std::conditional_variable->wait(). The manager coerces the thread out
 *	of the wait state with a conditional_variable->notify_all(). This
 *	automatically locks the threads internal mutex until procesisng is complete.
 *
 *	The manager starts all threads processing with a call to
 *	conditional_variable->notify_all(). The manager then queries the state of
 *	each thread. Because the thread already has the lock the manager is left
 *	in a blocked state until the thread completes the work.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_ThreadInstanceBoost
{
	public:
		enum	ENUM_TASKSTATE		{ STATE_TASK_COMPLETED,      STATE_TASK_PROCESSING,  STATE_NEW_TASK_REQUESTED};
		enum	ENUM_THREADSTATE	{ STATE_THREAD_INITIALIZING, STATE_THREAD_WAITING,   STATE_THREAD_BUSY, STATE_THREAD_CLOSED };

	private:
		volatile ENUM_THREADSTATE			m_threadstate;						// The current overall state of the thread. 
		volatile ENUM_TASKSTATE				m_taskstate;						// The current processing state of the task. This can be set by the manager asynchronously.
		volatile bool						m_keepprocessing;					// If true then continue processing data
		volatile bool						m_actionok;							// Result of the last call to ExecuteAction
		SKTRAN_ThreadManagerBoost*			m_manager;
		SKTRANSO_Quadrature_TLS_V21*			m_quadrature;						// The Sasktran Quadrature object, stores Sasktran thread specific context

		std::mutex						m_mutex;							// A mutex required for the condition_variable locking

	private:
		void								ProcessingLoop			();			// The main processing loop of the thread instance
		bool								CheckThreadIsWaiting	();


	public:
											SKTRAN_ThreadInstanceBoost	();
		virtual							   ~SKTRAN_ThreadInstanceBoost	();
		bool								SetManagerAndQuadrature	( SKTRAN_ThreadManagerBoost* manager, SKTRANSO_Quadrature_TLS_V21* quadrature );
		bool								SetOpticalProps			( const SKTRAN_TableOpticalProperties_V21*	optprop);

		void								SetTerminateFlag		();		//!< Manager calls this during a call to ThreadExecuteAction
		ENUM_TASKSTATE						TaskState				();		//!< Fetch the task state. This will block if the worker thread is working
		bool								SetNewTaskRequested		();		//!< Set the state machine (while the worker thread is asleep in conditional variable) that a new task is ready to be executed.
		bool								StartThread				();		//!< Starts the thread and immediately returns, Only call once during initialization
		bool								WaitUntilThreadIsWaiting();		//!< Waits until the thread is in a waiting state., This will block the current thread if the worker thread is working
		bool								ActionOk				()		{ return m_actionok;}

	friend class SKTRAN_ThreadInstanceBoostLauncher;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoostLauncher		2010-3-24*/
/** A very small class that boost libraries can use to luanch the thread.
 *	In essence boost requires an object that it can copy internally and 
 *	does not require the original instance to be hanging around.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_ThreadInstanceBoostLauncher								// A quick and dirty class
{																// that launches the thread
	private:
		SKTRAN_ThreadInstanceBoost*		m_thread;

	public:
									SKTRAN_ThreadInstanceBoostLauncher(SKTRAN_ThreadInstanceBoost* thread)	{ m_thread = thread;}
		void						operator()					 ()									{ m_thread->ProcessingLoop(); }
};



/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerBoost_CreateNewInstance		2011-5-31*/
/** This is a global function that is used to create new instances of the
 *	BoostThreadManager
 **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadManager* SKTRAN_ThreadManagerBoost_CreateNewInstance()
{
	SKTRAN_ThreadManager* manager;
	
	manager = new SKTRAN_ThreadManagerBoost;
	return manager;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::NumDefaultProcessingThreads		2010-5-26*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_ThreadManagerBoost::NumDefaultProcessingThreads()
{
	size_t	numcpu;
	numcpu = std::thread::hardware_concurrency();		// Get the number of CPU or cores on this machine
	if (numcpu == 0) numcpu = 1;
	return numcpu;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::SKTRAN_ThreadManager		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadManagerBoost::SKTRAN_ThreadManagerBoost()
{
	m_nextpoint = 0;
	m_threads = NULL;
	m_numthreads = 0;


}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::~SKTRAN_ThreadManager		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadManagerBoost::~SKTRAN_ThreadManagerBoost	()
{
	CloseThreadsAndRelease();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::ReleaseResources		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_ThreadManagerBoost::ReleaseResources()
{
	if (m_threads != NULL) delete [] m_threads;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager_Uniform::CloseThreadsDown		2010-3-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerBoost::CloseThreadsAndRelease()
{
	size_t	i;
	bool	ok;

	for (i = 0; i < m_numthreads; i++)
	{
		m_threads[i].SetTerminateFlag();
	}
	ok = NotifyWorkerThreadsAndWaitForEnd( 0 );
	ReleaseResources();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::SetOpticalProps		2010-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerBoost::SetOpticalProps( const SKTRAN_TableOpticalProperties_V21* optprop)
{
	bool	ok = true;
	bool	ok1;
	size_t	i;

	for (i = 0; i < m_numthreads; ++i)
	{
		ok1 = m_threads[i].SetOpticalProps(optprop);
		ok = ok && ok1;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_ThreadManager::SetOpticalProps, Error setting the optical properties for the worker threads");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerBoost::CreateThreads		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerBoost::CreateThreads(   size_t								numthreads,
										         const SKTRAN_SpecsInternal_V21*		modelspecifications,
										         const SKTRAN_EngineDiffuseTables*		modeltables   )
{
	bool									ok;
	bool									ok1;
	size_t									i;
	SKTRANSO_Quadrature_TLS_V21*				quadrature;
	const SKTRAN_Quadrature_Factory_V21*		quadraturefactory;

	if (numthreads == 0)
	{
		numthreads = NumDefaultProcessingThreads();
	}
	CloseThreadsAndRelease();
	NXTRACE(("SKTRAN_ThreadManagerBOost::CreateThreads, creating %u threads \n", (unsigned int)numthreads));

	m_threads = new SKTRAN_ThreadInstanceBoost[numthreads];
	ok = (m_threads != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, " SKTRAN_ThreadManager::CreateThreads, Error allocating memory for %Iu threads", (size_t)numthreads );
	}
	else
	{
		m_numthreads      = numthreads;
		quadraturefactory = modelspecifications->QuadratureSpecs()->QuadratureFactory();
		for (i = 0; i < m_numthreads; ++i)
		{
			ok1  = quadraturefactory->CreateThreadInstance( modelspecifications, modeltables, &quadrature );				// Create a thread specific quadrarture (holds all of the thread sepcific conext)
			ok1  = ok1 && m_threads[i].SetManagerAndQuadrature( this, quadrature );		// Let each thread know who the manager is and give it its own thread specific context
			ok1  = ok1 && m_threads[i].StartThread();									// start the thread
			ok   = ok && ok1;															// make sure everything is ok
			if (quadrature != NULL) quadrature->Release();								// and release our hold on the quadrature object.
		}
		if (ok)																			// We now have all the threads running
		{																				// NOw make sure they are all
			for (i = 0; i < m_numthreads; ++i)											// in the  waiting state
			{																			// before we start sending any commands
				ok1 = m_threads[i].WaitUntilThreadIsWaiting();							// otherwise signals might get lost
				ok = ok && ok1;
			}
		}
		NXTRACE(("SKTRAN_ThreadManager::CreateThreads, finished creating %u threads \n", (unsigned int)numthreads));

		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_ThreadManager::AddThreadInstance, Error adding the thread instance to our list, it will be destroyed");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerBoost::GetNextPointIndex		2010-3-22*/
/** Get the index of the next point in the table to process. Return false
 *	if there no more points to process else return the index of the next point
 *	to process.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerBoost::GetNextPointIndex( size_t* pointindex )
{
	bool	ok;
	std::lock_guard<std::mutex>		lock( m_nextpointmutex );

	ok =       m_nextpoint < NumPoints();
	*pointindex = m_nextpoint++;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerBoost::ResetPointIndex		2010-3-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerBoost::ResetPointIndex( size_t numpoints )
{
	std::lock_guard<std::mutex>		lock( m_nextpointmutex );

	ResetNumPoints(numpoints);
	m_nextpoint      = 0;
	return  true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerBoost::NotifyWorkerThreadsAndWaitForEnd		2010-3-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerBoost::NotifyWorkerThreadsAndWaitForEnd( size_t numpoints )
{
	bool		ok  = true;
	bool		ok1;
	size_t		i;
//	bool		actionok = true;

	ResetPointIndex( numpoints);												// Reset the point/ray counter 
	for (i = 0; i < m_numthreads; i++)											// Now tell each thread
	{																			// that 
		m_threads[i].SetNewTaskRequested();										// there is a new task to process.
	}																			// do all of the threads

//	NXTRACE(("***** MANAGER start waiting for threads\n"));
	m_threadsync.notify_all();													// Wake up all of the child threads waiting for a notification

	for (i = 0; i < m_numthreads; i++)											// Now call each thread
	{																			// although this looks like an inefficient spinlock it actually keeps this thread asleep most of the time
		ok1 = m_threads[i].WaitUntilThreadIsWaiting();							// as it is constantly blocked while it tries to lock the workerthreads internal mutex.
		ok = ok && ok1 &&  m_threads[i].ActionOk();								// Return the status of the overall process including the action executed.
	}
//	NXTRACE(("***** MANAGER end waiting for threads\n"));
//	NXASSERT(( ok ));
	return ok;
}






/*-----------------------------------------------------------------------------

 *					SKTRAN_ThreadInstanceBoost::SKTRAN_ThreadInstanceBoost		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadInstanceBoost::SKTRAN_ThreadInstanceBoost()
					  : m_mutex()
{
	m_taskstate      = STATE_TASK_PROCESSING;								// Flag that we have just complected a task
	m_threadstate    = STATE_THREAD_INITIALIZING;							// flag that 
	m_keepprocessing = true;												// This is the flag that lets us quit the threaded loop
	m_actionok       = true;
	m_manager        = NULL;
	m_quadrature     = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::~SKTRAN_ThreadInstanceBoost		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadInstanceBoost::~SKTRAN_ThreadInstanceBoost()
{
	if (m_quadrature != NULL) m_quadrature->Release();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::SetManagerAndQuadrature		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadInstanceBoost::SetManagerAndQuadrature( SKTRAN_ThreadManagerBoost* manager, SKTRANSO_Quadrature_TLS_V21* quadrature )
{
	quadrature->AddRef();
	m_manager = manager;
	if (m_quadrature != NULL) m_quadrature->Release();
	m_quadrature = quadrature;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::SetOpticalProps		2010-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadInstanceBoost::SetOpticalProps( const SKTRAN_TableOpticalProperties_V21*	optprop)
{
	NXASSERT(( m_quadrature != NULL ));
	return m_quadrature->SetOpticalProps( optprop);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::StartThread		2010-3-24*/
/** Starts the thread. The boos implementation requires that we use 
 *	class with the function operator overload (operator()). The object is
 *	also internally copied by the thread. **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadInstanceBoost::StartThread	()
{
	bool							ok;

	ok = (m_threadstate == STATE_THREAD_INITIALIZING);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_ThreadInstanceBoost::StartThread, This thread has already started. It is an error to start it again");
	}
	else
	{
		ok = (m_manager != NULL);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_ThreadInstanceBoost::StartThread, You must successfully call SetManager before starting the worker thread");
		}
		else
		{
			m_keepprocessing = true;								// Reset the keep processing flag. The flag
			std::thread( SKTRAN_ThreadInstanceBoostLauncher(this) );							// This is the point at which the new thread starts in SKTRAN_ThreadInstanceBoostLauncher::operator()
																	// which immediately calls SKTRAN_ThreadInstanceBoost::ProcessingLoop()
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::TaskState		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadInstanceBoost::ENUM_TASKSTATE	SKTRAN_ThreadInstanceBoost::TaskState()
{
	ENUM_TASKSTATE	state;
	std::lock_guard< std::mutex >			lock(m_mutex);

	state = m_taskstate;
	return state;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::SetNewTaskRequested		2010-3-24*/
/** This is a request from an external thread to start a new task. The actual
 *	must be configured in the derived class via virtual function ExecuteAction
	before calling this guy.
	**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadInstanceBoost::SetNewTaskRequested	()
{
	std::lock_guard< std::mutex >			lock(m_mutex);

	m_taskstate = STATE_NEW_TASK_REQUESTED; 
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::ProcessingLoop		2010-3-24*/
/** This is the main processing loop of the thread. The thread stays in a loop running
 *	couple of state machines. The purpose of the code is to wait for a command
 *	execute the command and wait for next command (and so on).  The 
 *	"wait for command" state is controlled by the thread manager's conditional_variable. 
 *
 *	Note that the loop locks the instance's m_mutex for all of this method apart from the
 *	call to ConditionalVariable()->wait where the mutex is internally released. This generates
 *	the correct behaviour as the worker thread is flagged as inaccessible while it is 
 *	"working" and becomes accessible once it is is "idle".
 **/
/*---------------------------------------------------------------------------*/

void  SKTRAN_ThreadInstanceBoost::ProcessingLoop	()
{
	std::unique_lock< std::mutex >			lock(m_mutex);
	bool										domore;
	size_t										pointindex;
	size_t										numpoints;
	bool										ok;

	NXTRACE(("SKTRAN_ThreadInstanceBoost::ProcessingLoop, Thread for instance 0x%p has started\n", (void *)this));

	#if defined(NX_WINDOWS)
	int threadp = ::GetThreadPriority( ::GetCurrentThread() );					// Change priority. This can be #define away as its purpose is purely cosmetic.
	::SetThreadPriority( ::GetCurrentThread(), THREAD_PRIORITY_LOWEST );		// It helps the response time for other programs if we run this code at low priority.
	#endif

	m_taskstate   = STATE_TASK_COMPLETED;
	while (m_keepprocessing)	
	{
		while (m_taskstate != STATE_NEW_TASK_REQUESTED )			// If the user has not requested a new reques
		{															// then
			NXASSERT(( lock.owns_lock() ));
			m_threadstate = STATE_THREAD_WAITING;					// Flag that this thread is waiting. The manager must wait for this condition using a mutexed lock before issuing any commands
			m_manager->ConditionalVariable()->wait( lock );		// Wait for notify_all from manager, note mutex is unlocked so manager can get state info. Note that condition_variables may exit without proper notification
			NXASSERT(( lock.owns_lock() ));
			m_threadstate = STATE_THREAD_BUSY;								// Flag that we are busy.  Manager cannnot get this info until we unlock the lock
		}																	// Make sure user has requested a new task, sometimes it may not happen
		m_taskstate   = STATE_TASK_PROCESSING;								// A new task is requested so flag we are processing the state
		if (m_keepprocessing)											
		{
//			NXTRACE(("***** Thread for instance 0x%p Starting execute action\n", (void *)this));
			domore       = m_manager->GetNextPointIndex( &pointindex );							// Get the next point in the table to process
			numpoints    = 0;
			m_actionok   = true;
			while (domore)
			{
				ok     = m_manager->ThreadExecuteAction(m_quadrature, pointindex);				// Execute the action provided by the manager ( this method MUST BE THREAD SAFE !!!!)
				domore = m_manager->GetNextPointIndex( &pointindex );
				numpoints++;
				m_actionok = m_actionok && ok;
			}
			m_quadrature->SetNumPointsProcessed( numpoints);
//			NXTRACE(("***** Thread for instance 0x%p Finishing execute action\n", (void *)this));
		}
		m_taskstate   = STATE_TASK_COMPLETED;
	}
	m_taskstate   = STATE_TASK_COMPLETED;
	m_threadstate = STATE_THREAD_CLOSED;					// Flag that this thread is waiting. The manager must wait for this condition using a mutexed lock before issuing any commands
	NXTRACE(("SKTRAN_ThreadInstanceBoost::ProcessingLoop, Thread for instance 0x%p has exited started\n", (void *)this));
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::CheckThreadIsWaiting		2010-5-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadInstanceBoost::CheckThreadIsWaiting()
{
	bool								iswaiting;
	std::lock_guard< std::mutex >	lock(m_mutex);				// Grab the mutex. This will often be delayed until the algorithm has finished do its work

	iswaiting = (m_taskstate == STATE_TASK_COMPLETED);

	return iswaiting;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::WaitUntilThreadIsWaiting		2010-3-24*/
/** This is a main synchronization routine that the manager calls to wait
 *	for each thread to finish its action. **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadInstanceBoost::WaitUntilThreadIsWaiting()
{
	bool	iswaiting;

	do
	{
		iswaiting = CheckThreadIsWaiting();
	} while (!iswaiting);													// keep looping until the worker thread is in the waiting state
	NXASSERT((  (m_threadstate ==STATE_THREAD_WAITING) || ( m_threadstate == STATE_THREAD_CLOSED) ));
	return iswaiting;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::TerminateThread		2010-3-23*/
/** The manager can only exit the threads by invoking ExecuteAction and
 *	then calling SetTerminateFlag from within that function. Upon return from
 *	the ExecuteAction function the thread will terminate.
 **/
/*---------------------------------------------------------------------------*/

void SKTRAN_ThreadInstanceBoost::SetTerminateFlag()
{
	m_keepprocessing = false;										// Force the Thread Instance to exit
}
