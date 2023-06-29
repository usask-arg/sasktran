/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

#include "nxbase_core.h"
//#include "nxlibCOM_i.c"		// Include the Definition of the COM InxLog interafces


InxLog*		nxLogBase::m_DefaultLogger = NULL;

//---------------------------------------------------------------------------
//						nxLogBase::QueryInterface
//	COM IUnknown support.  Simply delegate to controlling unknown.
//---------------------------------------------------------------------------

STDMETHODIMP nxLogBase::QueryInterface(REFIID riid, LPVOID FAR* ppvObj)
{
	HRESULT status;

	if ( (riid == IID_InxLog) || (riid == IID_IUnknown) )
	{
		*ppvObj = this;
		AddRef();
		status = S_OK;
	}
	else
	{
		status = E_NOINTERFACE;
		*ppvObj = NULL;
	}
	return status;
}

//---------------------------------------------------------------------------
//						nxLogBase::AddRef
//	COM IUnknown support. 
//---------------------------------------------------------------------------

STDMETHODIMP_(ULONG)  nxLogBase::AddRef()
{
	return m_cref++;
}

//---------------------------------------------------------------------------
//						nxLogBase::Release
//	COM IUnknown support. Note that we normally wont do this if the object
//	is created statially or on the stack
//---------------------------------------------------------------------------

STDMETHODIMP_(ULONG)  nxLogBase::Release()
{
	if (--m_cref > 0) return m_cref;
	delete this;
	return 0;
}


//----------------------------------------------------------------------------
//                      nxLogBase::Constructor
//	Normally the object is created on the stack.  For that special reason
//	make a call to AddRef() to keep the object alive. This can cause problems
//  if another object tries to keep it alive longer than it is actually alive.
//	A potential bug waiting to kill me but ...  I do try and track this in the
//	destructor.  If the reference count is greater than 1 then I throw an
//	exception.
//----------------------------------------------------------------------------

nxLogBase::nxLogBase()
{
	m_cref      = 1;
	m_isverbose = nxFALSE;
}

//----------------------------------------------------------------------------
//                      nxLogBase::Destructor
//	The InxLog interface is a "kluge".  Much of my legacy code generates the
//	nxlog object either staically or on the stack.  This is contrary to normal
//	COM objects which are created on the heap.  Consequently if you get an
//	exception that brings you to this destructor then you are deleting the
//	nxLog object (probably off the stack!) while other objects in the program still have a
//	reference count on this logger.  This is really bad news which is why the
//	exception is thrown.  Either get the other objects to Release the logger or
//	keep the logger alive longer.
//----------------------------------------------------------------------------

nxLogBase::~nxLogBase()
{
	if (m_DefaultLogger == this) m_DefaultLogger = NULL;
	//if (m_cref > 1) throw("nxLog Object Still referenced by other objects!!!!");
	/* ****** SEE THE COMMENTS ABOVE IF YOU HAVE THROWN AN EXCEPTION TO THIS POINT **** */
}

void nxLogBase::CheckDefaultLogger()
{
	if (m_DefaultLogger  == NULL) SetAsDefaultLogger();
}

//----------------------------------------------------------------------------
//                      nxLogBase::SetAsDefault
//	Let this instance of nxLog be the default logger.
//----------------------------------------------------------------------------

void nxLogBase::SetAsDefaultLogger()
{
	m_DefaultLogger = (InxLog*)this;
}

//----------------------------------------------------------------------------
//                      nxLogBase::SetAsDefault
//	Let this instance of nxLog be the default logger.
//----------------------------------------------------------------------------

void nxLogBase::SetAsDefaultLogger( InxLog* newlog)
{
	m_DefaultLogger = newlog;
}


// ---------------------------------------------------------------------------
//                       nxLogBase::Record( Format, ... )
// ---------------------------------------------------------------------------

void nxLogBase::Record( nxLogStatus status, const char * file, int line, const char *Format, ... )
{
	nxString	  message;
	if (m_DefaultLogger != NULL)
	{
		va_list ArgPtr;
	    va_start(ArgPtr, Format);
		message.vsprintf( Format, ArgPtr);								// and print the string to the nxString variable.
		m_DefaultLogger->Trace( status, file, line, message);
	    va_end(ArgPtr);
	}
#if defined(NXDEBUG)
	else
	{
 		va_list ArgPtr;
	    va_start(ArgPtr, Format);
		message.vsprintf( Format, ArgPtr);				// and print the string to the nxString variable.
	    va_end(ArgPtr);

		NXTRACE(("nxLogBase::Record, Default logger is undefined\n"));
		NXTRACE(("-------------> %s\n", (const char *)message));
	}
#endif

}

// ---------------------------------------------------------------------------
//                       nxLogBase::Record( Format, ... )
// ---------------------------------------------------------------------------

void nxLogBase::Verbose( nxLogStatus status,  const char * file, int line, const char *Format, ... )
{
	if (m_DefaultLogger != NULL)
	{
		if (m_DefaultLogger->IsVerbose() == S_OK)
		{
			nxString	message;
			va_list		ArgPtr;
			va_start(ArgPtr, Format);
			message.vsprintf( Format, ArgPtr);								// and print the string to the nxString variable.
			m_DefaultLogger->Trace( status, file, line, message);
			va_end(ArgPtr);
		}
	}
}

// ---------------------------------------------------------------------------
//                       nxLogBase::Record( Format, ... )
// ---------------------------------------------------------------------------

void nxLogBase::lrecord( nxLogStatus status,  const char * file, int line, const char *Format, ... )
{
	nxString	message;
	va_list ArgPtr;
	va_start(ArgPtr, Format);
	message.vsprintf( Format, ArgPtr);								// and print the string to the nxString variable.
	m_DefaultLogger->Trace( status, file, line, message );
	va_end(ArgPtr);
}

// ---------------------------------------------------------------------------
//                       nxLogBase::lverbose( Format, ... )
// ---------------------------------------------------------------------------

void nxLogBase::lverbose( nxLogStatus status,  const char * file, int line, const char *Format, ... )
{
	if ( IsVerbose() == S_OK )
	{
		nxString	message;
		va_list ArgPtr;
		va_start(ArgPtr, Format);
		message.vsprintf( Format, ArgPtr);								// and print the string to the nxString variable.
		m_DefaultLogger->Trace( status, file, line, message );
		va_end(ArgPtr);
	}
}


STDMETHODIMP nxLogBase::IsVerbose()
{
	if (m_isverbose) return S_OK;
	return S_FALSE;
} 

