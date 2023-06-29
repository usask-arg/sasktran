/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

#include "nxbase_core.h"
// --------------------------------------------------------------------------
//                               nxLogEntry::SetEntry
// --------------------------------------------------------------------------

void nxLogEntry::Set( int status, const char *mes, const char*file, int linenum )
{
	m_message = mes;
	m_status  = status;
	m_file    = file;
	m_line    = linenum;
	m_mjd.FromSystem();
}

// --------------------------------------------------------------------------
//                               nxLogEntry::StatusString
// --------------------------------------------------------------------------

nxString nxLogEntry::StatusString( nxBOOL padstring ) const
{
	nxString Buffer;

	if (padstring)
	{
		switch (m_status)
		{
			case nxlog_error        : Buffer = "ERROR\t:"; break;
			case nxlog_warning      : Buffer = "WARNING\t:"; break;
			case nxlog_info         : Buffer = "INFO\t:"; break;
			case nxlog_continuation	: Buffer = "\t"; break;
		}
	}
	else
	{
		switch (m_status)
		{
			case nxlog_error        : Buffer = "ERROR";   break;
			case nxlog_warning      : Buffer = "WARNING"; break;
			case nxlog_info         : Buffer = "INFO";    break;
			case nxlog_continuation	: Buffer = ""; break;
		}
	}
	return Buffer;
}

//---------------------------------------------------------------------------
//						nxLog::CreateMutex
//	Create a mutex to enable logging from multiple threads 
//---------------------------------------------------------------------------

void nxLog::CreateMutex()
{
#if defined (NX_WINDOWS)
	m_mutex =  ::CreateMutex( NULL, FALSE, NULL );
#endif
}

void nxLog::CloseMutex()
{
#if defined(NX_WINDOWS)
	CloseHandle( m_mutex);
#endif
}

//---------------------------------------------------------------------------
//						nxLog::Lock
//	Lock the mutex so only current thread can access the log.
//	There is a trick here when we are dealing with loggers drawn in windows.
//	A deadlock problem arises if the thread that owns the logger window is
//	waiting for the mutex when another thread is actively writing another log
//	entry to that window. This is a deadlock as the window thread must be
//	servicing its message queue.  To get around this we use the YieldToSystem
//	to make sure all threads are servicing their message queues while waiting
//	for the mutex.
//---------------------------------------------------------------------------

void nxLog::Lock()
{
#if defined(NX_WINDOWS)
	DWORD status;
	if (m_mutex != NULL)
	{

		for (int i = 0; i < 20; i++)
		{
			status = ::WaitForSingleObject( m_mutex, 100 );			// wait for the mutex
			if (status == WAIT_OBJECT_0) break;						// But I could be locked up because someone is sending a message to me
			YieldToSystem( 20 );									// so if 
		}
		if (status != WAIT_OBJECT_0)
		{
			if (status == WAIT_TIMEOUT)        ::MessageBox( NULL, "nxLog::Lock(), timed out (2000 ms) waiting for m_mutex", "nxLog", MB_OK);
			else if (status == WAIT_ABANDONED) ::MessageBox( NULL, "nxLog::Lock(), wait abandoned (threads finished!)" , "nxLog", MB_OK);
			else                               ::MessageBox( NULL, "nxLog::Lock(), Unidentified error waiting for m_mutex", "nxLog", MB_OK);
		}
	}
#endif

}

//---------------------------------------------------------------------------
//						nxLog::Unlock
//	Release the mutex so the other threads can log to  this guy.
//---------------------------------------------------------------------------

void nxLog::Unlock()
{
#if defined(NX_WINDOWS)
	if (m_mutex != NULL) ::ReleaseMutex(m_mutex);
#endif
}



/*---------------------------------------------------------------------------
 *						nxLog::privateinit
 *	Common constructor code.
 *-------------------------------------------------------------------------*/

void nxLog::privateinit()
{
	CreateMutex();
	m_mode = NXLOG_WINDOW;
	memset(m_logName, 0, sizeof(m_logName) ); 
}

//----------------------------------------------------------------------------
//                      nxLog::Constructor
//	Normally the object is created on the stack.  For that special reason
//	make a call to AddRef() to keep the object alive. This can cause problems
//  if another object tries to keep it alive longer than it is actually alive.
//	A potential bug waiting to kill me but ...  I do try and track this in the
//	destructor.  If the reference count is greater than 1 then I throw an
//	exception.
//----------------------------------------------------------------------------

nxLog::nxLog()
{
	privateinit();
}

//----------------------------------------------------------------------------
//                      nxLog::Constructor
//----------------------------------------------------------------------------

nxLog::nxLog( const char* logname)
{
	privateinit();
	SetLogName( logname );
}

//----------------------------------------------------------------------------
//                      nxLog::Destructor
//	The InxLog interface is a "kluge".  Much of my legacy code generates the
//	nxlog object either staically or on the stack.  This is contrary to normal
//	COM objects which are created on the heap.  Consequently if you get an
//	exception that brings you to this destructor then you are deleting the
//	nxLog object (probably off the stack!) while other objects in the program still have a
//	reference count on this logger.  This is really bad news which is why the
//	exception is thrown.  Either get the other objects to Release the logger or
//	keep the logger alive longer.
//----------------------------------------------------------------------------

nxLog::~nxLog()
{
	CloseMutex();
}

//----------------------------------------------------------------------------
//+
//NAME:
//                      nxLog::setLogName( NameStr )
// Sets the full file name of the log file to be used on subsequent logging calls
//
//-
//----------------------------------------------------------------------------

STDMETHODIMP nxLog::SetLogName(const char* NameStr)
{
	if (NameStr != NULL)
	{
		strncpy( m_logName, NameStr, MAX_PATH );
		m_logName[MAX_PATH] = '\0';
		SetBit( m_mode, NXLOG_FILE);
	}
	else
	{
		strcpy( m_logName, "" );
		ClearBit( m_mode, NXLOG_FILE );
	}
	return S_OK;
}
//---------------------------------------------------------------------------
//						nxLog::GetLogName
//---------------------------------------------------------------------------

STDMETHODIMP nxLog::GetLogName( const char** logname )
{
	*logname = m_logName;
	return S_OK;
}

//----------------------------------------------------------------------------
//                      nxLog::setLogMode( LogMode )
// Sets the log file mode to one of the defined values: NXLOG_FILE,
// NXLOG_WINDOW.
//----------------------------------------------------------------------------

STDMETHODIMP nxLog::SetLogMode(int  LogMode)
{
	SetBit( m_mode, LogMode);
	return S_OK;
}

//----------------------------------------------------------------------------
//                      nxLog::ClearLogMode( LogMode )
// Clears the logging mode using a bitwaise operation.
//----------------------------------------------------------------------------

STDMETHODIMP nxLog::ClearLogMode( int logmode )
{
	ClearBit( m_mode, logmode);
	return S_OK;
}


// ---------------------------------------------------------------------------
//                         nxLog::Trace
// Records the following information to output devices
// ---------------------------------------------------------------------------

STDMETHODIMP nxLog::Trace( int status,  const char * file, int line, const char *message )
{
	nxStringArray Lines;
	int			  nlines;
	nxLogEntry	  entry;

	Lines.Strtok( message, "\n" );						// Look for new lines in the printed string
	nlines = Lines.GetSize();							// Get the total number of lines
	if (Lines[nlines-1].IsEmpty())						// if the last line is empty it means a trailing
	{													// \n charcater in the format string
		nlines--;										// so just dump it as it just wastes space
		if (nlines < 1) nlines = 1;						// but always send at least one line
	}

	Lock();												// Lock this code with the mutex so we dont get multiple thread problems.
	for (int i =0; i < nlines; i++ )					// for all the lines in the string
	{													// output each line
		nxString& aline = Lines[i];						// get the line
		entry.Set( status, aline, file, line );			// Set the entry
		if ((m_mode & NXLOG_FILE)   == NXLOG_FILE)		// if file output enabled
		{												// then
			WriteEntryToFile( entry );					// output to file
		}												//

		if ((m_mode & NXLOG_WINDOW)	== NXLOG_WINDOW)	// If console display is required
		{												// then display it.
			DisplayEntry( entry );						// on the display
		}
		status = nxlog_continuation;					// flag extra lines as a continuation of this line.
	}
	Unlock();
	return S_OK;
}

//----------------------------------------------------------------------------
//			nxLog::WriteEntryToFile()
//----------------------------------------------------------------------------

void nxLog::WriteEntryToFile( const nxLogEntry& entry )
{
	nxString	TheLine;
	FILE*		LOGFile;

	if (strlen(m_logName) > 0)
	{
		TheLine = "[";
		TheLine += entry.Mjd()+"]" + entry.StatusString( nxTRUE ) + entry.Message();

		LOGFile = fopen(m_logName,"at");
		if (LOGFile != NULL)
		{
			fprintf(LOGFile,"%s\n", (const char*)TheLine);
			fclose(LOGFile);										// always close file after use.
		}
#if defined(NXDEBUG)
		else
		{
			NXTRACE(("nxLog::WriteEntryToFile, Error opening file %s", (const char*)m_logName));
		}
#endif
	}
}

