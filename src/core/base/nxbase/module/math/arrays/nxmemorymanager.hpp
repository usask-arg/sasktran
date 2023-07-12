

/*---------------------------------------------------------------------------
 *'					InxMemoryManager<T>::InxMemoryManager		2003-9-16
 *-------------------------------------------------------------------------*/

template <class T >
InxMemoryManager<T>::InxMemoryManager()
{
#if defined(NXDEBUG)
	nxMemoryManagerInitCode();
	g_nxMemoryManager_numinstance++;
//	NXTRACE(("InxMemoryManager CONSTRUCTOR address %08x, # instances created = %d\n", (void *)this, (int)g_nxMemoryManager_numinstance));
#endif


	m_lockcount       = 0;
	m_memorylockcount = 0;
	m_linearsize      = 0;
	m_pointer         = NULL;
	m_endpointer      = NULL;
	m_reservesize     = 0;
}


/*---------------------------------------------------------------------------
 *'					InxMemoryManager<T>::~InxMemoryManager		2003-9-16
 *-------------------------------------------------------------------------*/

template <class T >
InxMemoryManager<T>::~InxMemoryManager	()
{
	if (( m_memorylockcount != 0) && (m_pointer != NULL))
	{
		nxLog::Record(NXLOG_WARNING,"InxMemoryManager<T>::Destructor, there are still reference counts on the memory. It will be destroyed");
		InternalFreemem();
	}
#if defined(NXDEBUG)
	g_nxMemoryManager_numinstance--;
//	NXTRACE(("InxMemoryManager DESTRUCTOR address %08x, # instances left = %d\n", (void *)this, (int)g_nxMemoryManager_numinstance));
#endif

}


//---------------------------------------------------------------------------
//						InxMemoryManager<T>::InternalAllocate
//	Default implementation to  allocate an array of "linearsize" objects.
//	By this stage it is guaranteed that linearsize is not 0 (or less).
//	The code may return NULL pointers if the meory allocation does not work.
//---------------------------------------------------------------------------

template <class T >
T* InxMemoryManager<T>::InternalAllocate( size_t linearsize )
{
	T* array = new T [linearsize];							// if we have some data then allocate
	return array;										// and return the array.
}

/*---------------------------------------------------------------------------
 *'					InxMemoryManager<T>::InternalFree        2002-6-14
 *-------------------------------------------------------------------------*/

template <class T>
void InxMemoryManager<T>::InternalFreemem()
{
	if (m_pointer != NULL)
	{
		delete [] m_pointer;							// release the memory
	}
	m_pointer     = NULL;								// Set the pointer equal to NULL
}


/*---------------------------------------------------------------------------
 *'					InxMemoryManager<T>::Erase		2003-9-16
 *-------------------------------------------------------------------------*/

template <class T>
void InxMemoryManager<T>::Erase()
{
	if (m_pointer != NULL)
	{
		InternalFreemem();
	}
	m_pointer     = NULL;								// Set the pointer equal to NULL
	m_endpointer  = NULL;
	m_reservesize  = 0;
	m_linearsize  = 0;
	m_memorylockcount = 0;
}




/*---------------------------------------------------------------------------
 *'					InxMemoryManager<T>::UnlockMemory		2003-9-16
 *-------------------------------------------------------------------------*/
template <class T>
void InxMemoryManager<T>::UnlockMemory()
{
	if (m_pointer != NULL)
	{
		if (--m_memorylockcount <= 0) Erase();
	}
}


/*---------------------------------------------------------------------------
 *'					InxMemoryManager<T>::Allocate                2002-6-14
 *	Allocate nelements of memory and sets the memory lock count to 1.
 *
 *	1)	The allocation will do nothing if the new size is identical to the
 *		current size. But is will return SUCCESS.
 *
 *	2)	The allocation will use existing memory allocation if it fits within
 *		the bounds of existing memory allocation and there are no excess locks on
 *		memory and the user is requesting that we useexistingallocation and
 *		the memory manager supports memory re-use in this manner.
 *
 *	3)	Finally it will allocate memory using our memory manager.  If succesful
 *		It will place a lock on the memory.
 *
 *	If the function return FAILURE then the caller should call UnlockMemory()
 *	to clear any locks the user may have on the memory.
 *-------------------------------------------------------------------------*/

template  <class T>
nxBOOL InxMemoryManager<T>::AllocateAndLock( const RankSpecification* rankspecs, nxBOOL useexistingallocation, T** ptr )
{
	nxBOOL	ok;
	size_t		linearsize = rankspecs->GetContiguousStorageSize(sizeof(T));

	ok =  (linearsize == m_linearsize);												// Is this array already the correct size
	if (!ok)																		// Nope
	{																				// so lets try re-sizing;
		ok = (m_memorylockcount <= 1);												// and resizing can only done if we have one and only one (we assume the one is held by our caller)
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "InxMemoryManager<T>::Allocate, cannot allocate as there are %d locks on the memory", (int)m_memorylockcount);
		}
		else
		{
			useexistingallocation = useexistingallocation && InternalAllowRealloc();	// Only allow the array to grow if this memory allocator supports the operation
			ok = (useexistingallocation) && ( linearsize <= m_reservesize);				// If we can re-size arrays and current allocation is big enough
			if (ok)																		// then
			{																			// simply
				m_linearsize  = linearsize;												// Set the linear size
				m_endpointer  = m_pointer + linearsize;									// and the endpointer
			}
			else
			{
				Erase();
				m_pointer = InternalAllocate(linearsize);							// allocate the new space
				ok = (m_pointer != NULL );											// Assert that everything is ok.
				if (!ok)															// and that check that its ok
				{
					nxLog::Record(NXLOG_WARNING, "InxMemoryManager<T>::Allocate, Insufficient memory to allocate %d elements", (int)linearsize );
				}
				else
				{
					m_linearsize          = linearsize;				// Set the linear size
					m_reservesize         = linearsize;
					m_endpointer          = m_pointer + linearsize;
					m_memorylockcount     = 1;
				}
			}
		}																// and that is that
	}
	*ptr = m_pointer;
	return ok;
}

