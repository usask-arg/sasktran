
#include "nxbase_core.h"
#include "nxbase_threads.h"

static std::mutex g_nxUnknown_refcountlock;				// then add a mutex for AddRef and Release locking purposes.


#if defined (NXDEBUG)
size_t g_nxUnknown_numinstance = 0;
nxBOOL g_nxUnknown_isinitialized = nxFALSE;

static void nxUnknownExitCode()
{
	if ( g_nxUnknown_numinstance != 0)
	{
		NXTRACE(("nxUnknown:: Program exiting with %d reference counts on objects\n", (int)g_nxUnknown_numinstance));
	}
}

static void nxUnknownInitCode()
{
	if (!g_nxUnknown_isinitialized)
	{
		atexit(nxUnknownExitCode);
		g_nxUnknown_numinstance   = 0;
		g_nxUnknown_isinitialized = nxTRUE;
	}
}
#endif


/*-----------------------------------------------------------------------------
 *					nxUnknown::nxUnknown		2007-12-5*/
/** **/
/*---------------------------------------------------------------------------*/

nxUnknown::nxUnknown()
{
#if defined(NXDEBUG)
	nxUnknownInitCode();
	g_nxUnknown_numinstance++;
//	NXTRACE(("nxUnknown CONSTRUCTOR address %08x, # instances created = %d\n", (void *)this, (int)g_nxUnknown_numinstance));
#endif

	m_cref = 0;
	m_isstatic = nxFALSE;
}


/*---------------------------------------------------------------------------
 *'					nxUnknown::SetStatic		2003-12-11
 *		This is used to signify the core object is statically allocated
 *	and must not be deleted.
 *-------------------------------------------------------------------------*/

void nxUnknown::SetStatic()
{
	m_isstatic = nxTRUE;
#if defined(NXDEBUG)
	g_nxUnknown_numinstance--;
#endif
	AddRef();
}


/*-----------------------------------------------------------------------------
 *					nxUnknown::~nxUnknown		2007-12-5*/
/** **/
/*---------------------------------------------------------------------------*/

nxUnknown::~nxUnknown()
{
#if defined(NXDEBUG)
	ULONG	n;

	n = (m_isstatic) ? 1 : 0;
	if (m_cref != (LONG)n)
	{
		NXTRACE(("nxUnknown::Destroying an object with excess reference counts still on it\n"));
	}
	if (!m_isstatic) g_nxUnknown_numinstance--;
//	NXTRACE(("nxUnknown DESTRUCTOR address %08x, # instances left = %d\n", (void *)this, (int)g_nxUnknown_numinstance));
#endif

}


/*-----------------------------------------------------------------------------
 *					nxUnknown::AddRef		2009-5-29*/
/** There is now support for boost. 
 **/
/*---------------------------------------------------------------------------*/

LONG nxUnknown::AddRef()
{
	{
		std::unique_lock<std::mutex>	lock( g_nxUnknown_refcountlock );		// lock the mutex, There must be a better way than this!!!
		++m_cref;																	// increment the reference count
	}
//	::InterlockedIncrement(&m_cref);											//  m_cref++, using Windows multithread safe function
	return m_cref;
}


/*-----------------------------------------------------------------------------
 *					nxUnknown::Release		2009-5-29*/
/** There is now support for boost. The library must be built with
 *	NX_USEBOOST_THREADLIBS
 **/
/*---------------------------------------------------------------------------*/

LONG nxUnknown::Release()
{
	{
		std::unique_lock<std::mutex>	lock( g_nxUnknown_refcountlock );		// lock the mutex, There must be a better way than this!!!
		--m_cref;																	// decrement the reference count
	}
//	::InterlockedDecrement(&m_cref);												//  m_cref--, using Windows multithread safe function
	if ( m_cref > 0) return m_cref;
	if (!m_isstatic ) delete this;
	return 0;
}

/*-----------------------------------------------------------------------------
 *					nxUnknown::AddRef		2009-5-29*/
/** This is a special form where we actually violate the const attribute
 *	as we have to modify the reference count even though the object is
 *	unmodifiable. A little dangerous but not too bad.
 **/
/*---------------------------------------------------------------------------*/


LONG nxUnknown::AddRef() const
{
	nxUnknown*	modifiablethis;
	modifiablethis = (nxUnknown*)(intptr_t)this;
	return modifiablethis->AddRef();
}


/*-----------------------------------------------------------------------------
 *					nxUnknown::Release		2009-5-29*/
/** This is a special form where we actually violate the const attribute
 *	as we have to modify the reference count even though the object is
 *	unmodifiable. A little dangerous but not too bad.
 **/
/*---------------------------------------------------------------------------*/

LONG nxUnknown::Release() const
{
	nxUnknown*	modifiablethis;
	modifiablethis = (nxUnknown*)(intptr_t)this;
	return modifiablethis->Release();
}





