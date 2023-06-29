/*----------------------------------------------------------------------------
 *				class nxSyncObject									2002-9-30*/
/**	\ingroup system_multithread
 *	Base class for thread synchronization objects. Uses
 *	MsgWaitForMultipleObjects to implement the synchronization
 *	It is possible to wait on messages as well
**/
/*------------------------------------------------------------------------*/

class  nxSyncObject
{
	private:
		DWORD					m_wakemask;

	public:
		HANDLE					m_hObject;

	public:
								nxSyncObject();
		virtual				   ~nxSyncObject();
		void					SetWakeMask(DWORD mask) {m_wakemask = mask;}
								operator HANDLE() const {return m_hObject;}

 // Operations
	virtual BOOL					Lock(DWORD dwTimeout = INFINITE);
	virtual BOOL					Unlock() = 0;
	virtual	BOOL					Unlock(LONG /* lCount */, LPLONG /* lpPrevCount=NULL */) { return TRUE; }

 // Implementation
	friend class nxSingleLock;
};


/*----------------------------------------------------------------------------
 *						class nxMutex										*/
/**	\ingroup system_multithread
 *	A class for mutexing threads.
**/
/*--------------------------------------------------------------------------*/

class  nxMutex : public nxSyncObject
{
	public:
						nxMutex( BOOL bInitiallyOwn = FALSE, LPCTSTR lpszName = NULL, LPSECURITY_ATTRIBUTES lpsaAttribute = NULL);

public:
	virtual			   ~nxMutex(){};
	BOOL				Unlock();
};

/*----------------------------------------------------------------------------
 *						class nxEvent										*/
/**	\ingroup system_multithread
 * A class for handling events between threads
**/
/*--------------------------------------------------------------------------*/

class  nxEvent : public nxSyncObject
{
	public:
				nxEvent(BOOL bInitiallyOwn = FALSE, BOOL bManualReset = FALSE, LPCTSTR lpszNAme = NULL, LPSECURITY_ATTRIBUTES lpsaAttribute = NULL);

	public:
		BOOL	SetEvent()	{ return ::SetEvent(m_hObject); }
		BOOL	Unlock();
		BOOL	PulseEvent(){ return ::PulseEvent(m_hObject); }
		BOOL	ResetEvent(){ return ::ResetEvent(m_hObject); }

	public:
		virtual ~nxEvent(){};
};

/*----------------------------------------------------------------------------
 *						class nxSingleLock									*/
/**	\ingroup system_multithread
 * A class for handling single thread locks while executing code.
**/
/*--------------------------------------------------------------------------*/

class  nxSingleLock
{
 // Constructors

	public:
						nxSingleLock(nxSyncObject* pObject, BOOL bInitialLock = FALSE);

 // Operations
	public:
		BOOL			Lock(DWORD dwTimeOut = INFINITE);
		BOOL			Unlock();
		BOOL			Unlock(LONG lCount, LPLONG lPrevCount = NULL);
		BOOL			IsLocked();

 // Implementation
	public:
		~nxSingleLock()	{Unlock();}

	protected:
		nxSyncObject*	m_pObject;
		HANDLE			m_hObject;
		BOOL			m_bAcquired;
};

