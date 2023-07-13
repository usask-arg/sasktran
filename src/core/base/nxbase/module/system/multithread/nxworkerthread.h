// boosttest.cpp : Defines the entry point for the console application.
//
#if !defined(NXWORKERTHREAD_INCLUDEH)
#define NXWORKERTHREAD_INCLUDEH

//#if defined(NX_USEBOOST_THREADLIBS)						// Unix versions must use Boost Multi-threading Libs
//#include <boost/thread.hpp>
#include <mutex>
#include <functional>
//#else
//#include <process.h>
//#endif

class nxWorkerThreadInstance;


/*-----------------------------------------------------------------------------
 *					class nxWorkerThreadManager		2009-5-28*/
/** **/
/*---------------------------------------------------------------------------*/

class nxWorkerThreadManager
{
	private:
		volatile int						m_numinstances;				// Number of threads currently active
		std::mutex							m_incmutex;
		std::condition_variable				m_threadactivecond;
//		nxMutex								m_incmutex;
//		nxEvent								m_threadactivecond;



	private:
		int									IncrementInstance			();
		int									DecrementInstance			();

	public:
		void								FinishWorkerThread			( nxWorkerThreadInstance* userthread );
		static size_t						GetCurrentThreadIdCode		();
	
	public:
											nxWorkerThreadManager		();
										   ~nxWorkerThreadManager		();
		bool								WaitForThreadsToFinish		();
		bool								StartWorkerThread			( nxWorkerThreadInstance* userthread );
};

/*-----------------------------------------------------------------------------
 *					nxThreadStorageMap		 2014- 11- 12*/
/** A templated class to provide simple thread safe storage for applications.
 *	The class avoids user supplied thread indexing by using the unique thread
 *	id code provided by the operating system as an index for thread specific data. This
 *	can be a big programming advantage as it allows the user to drop passing
 *	the current thread index around their code. All the programmer has to do is identify where thread local
 *	storage arrays are required, create an instance of class THREADDATA with suitable default
 *	constructors and destructors and use it to create an instance of nxThreadStorageMap. The user can specify 
 *	a function or method to be called to initialize thread storage after after THREADDATA construction use the std::bind function 
 *	to acquire a functional of the form "bool  f( THREADDATA *)" and passing this to method SetThreadStorageInitializer
 * 
**/
/*---------------------------------------------------------------------------*/

template <class THREADDATA>
class nxThreadStorageMap
{

	private:
		std::mutex										m_threadlock;					// Add a lock to make sure map.insert is thread safe.
		std::map<size_t, THREADDATA>					m_data;							//!< Thread storage. The integer key is the value returned by omp_get_thread_num
		std::function< bool(THREADDATA*) >				m_threadinitializer;			//A function pointer used to initialize thread storage when it is created.

	private:
		bool												LookupUpThreadDataInternal(THREADDATA** data );

	public:
		bool												LookupUpThreadData			(THREADDATA** data ) const;
		bool												Clear						()											{ m_threadlock.lock(); m_data.clear(); m_threadlock.unlock(); return true;}
		void												SetThreadStorageInitializer	( std::function< bool(THREADDATA*) > f)		{ m_threadinitializer = f;}
};


/*-----------------------------------------------------------------------------
 *					nxThreadStorageVector		 2014- 11- 12*/
/** **/
/*---------------------------------------------------------------------------*/

template <class THREADDATA>
class nxThreadStorageVector
{
	private:
		std::vector<THREADDATA>						m_data;							//!< Thread storage. 
		bool										LookupUpThreadDataInternal		(size_t threadindex, THREADDATA** data );

	public:
		bool										Resize							( size_t numentries)	{ m_data.resize(numentries); return m_data.size()== numentries;}
		bool										LookupUpThreadData				( size_t threadindex, THREADDATA** data ) const;
		bool										Clear							()						{ m_data.clear(); return true;}
		THREADDATA&									at								( size_t idx)			{ return m_data.at(idx);}
};


/*-----------------------------------------------------------------------------
 *					nxThreadStorage<THREADDATA>::LookupUpThreadDataInternal	2014- 10- 22*/
/** **/
/*---------------------------------------------------------------------------*/

template <class THREADDATA>
bool nxThreadStorageMap<THREADDATA>::LookupUpThreadDataInternal(THREADDATA** data )
{
	size_t												threadnum;
	typename std::map<size_t, THREADDATA>::iterator				iter;
	bool												ok;

	threadnum = nxWorkerThreadManager::GetCurrentThreadIdCode();
	m_threadlock.lock();													// We must lock the thread to ensure no-one else is using it while we find our data
	iter      = m_data.find(threadnum);												// Find our data
	ok = (iter != m_data.end() );													// See if we succeeded
	if (!ok)
	{
		std::pair< typename std::map<size_t, THREADDATA>::iterator, bool> result;
		result = m_data.insert(  std::pair<size_t, THREADDATA>(threadnum, THREADDATA()) );
		ok     = result.second;
		iter   = result.first;
		if (m_threadinitializer)
		{
			ok = ok && m_threadinitializer( &(*iter).second );
		}
	}
	m_threadlock.unlock();
	if (!ok)
	{
		*data = nullptr;
		nxLog::Record(NXLOG_WARNING,"nxThreadStorage::LookupUpThreadData, error fetching/creating thread local storage for thread id (%d)", (int)threadnum);
	}
	else
	{
		*data =  &(*iter).second;
	}
	return ok;

}

/*-----------------------------------------------------------------------------
 *					nxThreadStorage::LookupUpThreadData		 2014- 11- 12*/
/** **/
/*---------------------------------------------------------------------------*/

template <class THREADDATA>
bool nxThreadStorageMap<THREADDATA>::LookupUpThreadData(THREADDATA** data ) const
{
	nxThreadStorageMap*	ptr;

	ptr = const_cast< nxThreadStorageMap*>(this);
	return ptr->LookupUpThreadDataInternal( data);
}

/*-----------------------------------------------------------------------------
 *					nxThreadStorageVector<THREADDATA>::LookupUpThreadDataIternal		 2014- 11- 12*/
/** **/
/*---------------------------------------------------------------------------*/

template <class THREADDATA>
bool nxThreadStorageVector<THREADDATA>::LookupUpThreadDataInternal(size_t threadindex, THREADDATA** data )
{
	bool			ok;
	THREADDATA*		dataptr = nullptr;

	ok = threadindex < m_data.size();
	if (ok)
	{
		dataptr = &m_data.at(threadindex);
	}
	else
	{
		dataptr = nullptr;
		nxLog::Record(NXLOG_WARNING,"nxThreadStorageVector::LookupUpThreadData, the supplied thread index (%d) is outside the range of teh internal storage (%d). Make sure you have called method Resize successfully.", (int)threadindex, (int)m_data.size());
	}
	*data = dataptr;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxThreadStorageVector<THREADDATA>::LookupUpThreadData		 2014- 11- 12*/
/** **/
/*---------------------------------------------------------------------------*/

template <class THREADDATA>
bool nxThreadStorageVector<THREADDATA>::LookupUpThreadData(size_t threadindex, THREADDATA** data ) const
{
	nxThreadStorageVector<THREADDATA>*	ptr;

	ptr = const_cast< nxThreadStorageVector<THREADDATA>*>(this);
	return ptr->LookupUpThreadDataInternal( threadindex, data);
}

/*-----------------------------------------------------------------------------
 *					class nxWorkerThreadInstance					2009-5-28*/
/** This is a base class for a classes that implement a worker thread. The
 *	derived classes must implement ExecuteThread. The function ExecuteThread will
 *	run in the worker thread using the same instance as that passed into 
 *	nxWortkerThreadManager.StartWorkerThread. This code avoids an issue in
 *	std::thread where std::thread clones the instance to another instance
 *	and uses that instance in the thread. The copied instance is destroyed
 *	when the thread finishes along with any thread specific information.
 *
 *	\par Example
 *	I include an example below. The main program launches three worker threads
 *	The thread manager waits for the three threads to complete before quitting
 *	the program
 *
 *	\code
 *	class mythread : public nxWorkerThreadInstance
 *	{
 *		private:
 *			int				m_id;
 *
 *		public:
 *							mythread		(int id) { m_id = id;}
 *			virtual void	ExecuteThread	();
 *	};
 *
 *
 *	void mythread::ExecuteThread()
 *	{
 *		Do_the_thread_work();
 *	}
 *
 *
 *	void main()
 *	{
 *		nxWorkerThreadManager	threadmanager;
 *		mythread				workera(1);
 *		mythread				workerb(2);
 *		mythread				workerc(3);
 *
 *	threadmanager.StartWorkerThread( &workera );
 *	threadmanager.StartWorkerThread( &workerb );
 *	threadmanager.StartWorkerThread( &workerc );
 *	threadmanager.WaitForThreadsToFinish();
 *	printf("We are done\n");
 * }
 *	\endcode
**/

/*---------------------------------------------------------------------------*/

class nxWorkerThreadInstance
{
	private:
		bool					m_isrunning;
		nxWorkerThreadManager*	m_manager;

	public:
								nxWorkerThreadInstance	()									{m_isrunning = false; m_manager = NULL;}
		virtual				   ~nxWorkerThreadInstance	()									{}
		bool					IsRunning				()									{ return m_isrunning;}
		void					SetIsRunning			( bool isrunning )					{ m_isrunning = isrunning;}
		void					SetManager				( nxWorkerThreadManager* manager)	{ m_manager   = manager;}
		nxWorkerThreadManager*	Manager					()									{ return m_manager;}
		virtual void			ExecuteThread			() = 0;
};


/*-----------------------------------------------------------------------------
 *					nxThread_TLStorage		 2016- 6- 27*/
/** A template for using thread safe storage using a pure boost and std:: implementation
 *	
 **/
/*---------------------------------------------------------------------------*/

template< class STORAGE_OBJECT>
class nxThread_TLStorage
{

	private:
				std::mutex																m_mutex;			// A mutex object used to sync access to m_TLStorage 
		std::map<std::thread::id, std::shared_ptr<STORAGE_OBJECT> >						m_TLStorage;		// The list of regions, one for each thread, indexed by the thread::id 
		typedef typename std::map<std::thread::id, std::shared_ptr<STORAGE_OBJECT> >::value_type	value_type;

	public:

		std::shared_ptr<STORAGE_OBJECT>				Storage		();
		void										Erase		();
		bool										SetStorage	(STORAGE_OBJECT* region	);
		void										clear		();
};



/*-----------------------------------------------------------------------------
 *					nxThread_TLStorage::Storage		 2016- 6- 17*/
/** Get the current storage object for this thread. The object is returned as a shared_ptr
 *	that guarantees the thread storage object wil not be destroyed until the user destroys 
 *	the shared_ptr. The shared_ptr may return a nullptr. The function will return the 
 *	object set by the last call in this thread to SetStorage. Users should not try to delete
 *	the thread storage themselves but should rely upon the shared_ptr reference count to destroy the 
 *	object .
**/
/*---------------------------------------------------------------------------*/

template<class STORAGE_OBJECT>
std::shared_ptr<STORAGE_OBJECT>	nxThread_TLStorage<STORAGE_OBJECT>::Storage()
{
	std::unique_lock<std::mutex>			lock(m_mutex);									// Lock the threads so I can acces sthe TLS storage

	auto x = m_TLStorage.find( std::this_thread::get_id());									// See if this thread local storage already exists for this thread. Get copy of object and increment the lock count if it does. (object lifetime guaranteed until user releases shared_ptr)
	return (x != m_TLStorage.end()) ? (x->second) : std::shared_ptr<STORAGE_OBJECT>(nullptr); 
}


/*-----------------------------------------------------------------------------
 *					nxThread_TLStorage::Erase		 2016- 6- 17*/
/** Removes the thread local storage associated with the current thread. The object
 *	will only be destroyed if all shared_ptr references are also released.
**/
/*---------------------------------------------------------------------------*/

template<class STORAGE_OBJECT>
void nxThread_TLStorage<STORAGE_OBJECT>::Erase()
{	
	std::unique_lock<std::mutex>			lock(m_mutex);						// Lock the threads so I can acces sthe TLS storage

	auto x = m_TLStorage.find( std::this_thread::get_id());				// See if this thread local storage already exists for this thread. Get copy of object and increment the lock count if it does. (object lifetime guaranteed until user releases shared_ptr)
	if (x != m_TLStorage.end()) m_TLStorage.erase(x);			// erase this entry. It will actually destroy when "x" goes out of scope
}

/*-----------------------------------------------------------------------------
 *					nxThread_TLStorage::SetStorage		 2016- 6- 17*/
/** Sets the thread local storage for this thread. The object will be accessible by 
*	subsequent calls to Storage. The STORAGE_OBJECT* pointed to by region  will be 
*	released  using "delete region" by this class so the object region must be allocated with "new"
 **/
/*---------------------------------------------------------------------------*/

template<class STORAGE_OBJECT>
bool nxThread_TLStorage<STORAGE_OBJECT>::SetStorage(STORAGE_OBJECT* region	)
{
	std::unique_lock<std::mutex>		lock(m_mutex);						// Lock the threads so I can access the TLS storage	

	auto status = m_TLStorage.insert( value_type(std::this_thread::get_id(),  std::shared_ptr<STORAGE_OBJECT>(region) ) );
	return status.second;
}


template<class STORAGE_OBJECT>
void nxThread_TLStorage<STORAGE_OBJECT>::clear()
{
	std::unique_lock<std::mutex>		lock(m_mutex);						// Lock the threads so I can access the TLS storage

	m_TLStorage.clear();
}



#endif






