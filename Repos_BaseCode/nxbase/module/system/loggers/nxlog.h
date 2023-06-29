#if !defined(NXBASE_NXLOG_H)
#define NXBASE_NXLOG_H 1

#include "../../system/cominterface/nxunknown.h"
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
*
* HISTORY
* -------
* 21-Dec-1999	NDL
*		Added method privateinit to class nxLog.  Used this so I had common
*		constructor code for the two nxLog constructors.  Follows a report
*		from Harri Auvinen at FMI regarding problems with compiler. Also
*		made changes to nxlog.cxx at this time.
*
****************************************************************************/

#include <stdarg.h>

EXTERN_C const IID IID_InxLog;



/*-----------------------------------------------------------------------------
 *					class InxLog									2004-11-23*/
/** \ingroup system_loggers
 *	This is the pure virtual base class for the logging classes. It used to be
 *	code in MIDL but it is now hard-coded as a class. This interface is
 *	identified by the GUID, IID_InxLog.
**/
/*---------------------------------------------------------------------------*/

class InxLog : public IUnknown
{
    public:
    	virtual							  ~InxLog(){}
        virtual HRESULT STDMETHODCALLTYPE Trace     	( int status, const char *sourcefilename, int linenum, const char *message) = 0;
        virtual HRESULT STDMETHODCALLTYPE IsVerbose 	( void ) = 0;
        virtual HRESULT STDMETHODCALLTYPE SetLogName	( const char *logfilename) = 0;
        virtual HRESULT STDMETHODCALLTYPE SetLogMode	( int mode) = 0;
        virtual HRESULT STDMETHODCALLTYPE ClearLogMode	( int mode) = 0;
        virtual HRESULT STDMETHODCALLTYPE GetLogName	(  const char **logfilename) = 0;

};

//DEFINE_GUID( IID_InxLog, 0xc3c922f0, 0x9a89, 0x11d2, 0xb8,0x98,0x00,0x00,0xc0,0x54,0x85,0x54);

enum nxLogMode    { NXLOG_FILE = 1,  NXLOG_WINDOW = 2};
enum nxLogStatus  { nxlog_error, nxlog_warning, nxlog_info, nxlog_continuation};

#define NXLOG_ERROR		nxlog_error,__FILE__,__LINE__
#define NXLOG_WARNING   nxlog_warning,__FILE__,__LINE__
#define NXLOG_INFO      nxlog_info,__FILE__,__LINE__


/*--------------------------------------------------------------------------
 *                              class nxLogBase								*/
/**	\ingroup system_loggers
 *	A tight implementation that can be included in COM objects for error
 *	reporting back to the main control program.  The program provides a simple
 *	interface that implements nxLog::Record and nxLog::Verbose.  The algorithm
 *	encodes the variable argument list as a string and then calls the DefaultLogger
 *	Trace function.  This internal implementation does nothing.  However if the default
 *	logger is overrided then the messages go to the new DafaultLogger implementation
 *	of Trace.  Base implenations are provided in class  nxLog::
**/
/*--------------------------------------------------------------------------*/

class nxLogBase : public InxLog
{
	protected:
		static InxLog*		m_DefaultLogger;										// The default logger if nothing else defined
		nxBOOL				m_isverbose;											// The logger is verbose (ie nxLog::Verbose doe ssomthing!
		UINT				m_cref;													// The interface reference count for this object

	 public:
							nxLogBase();											// Default constructor
		virtual			   ~nxLogBase();
		void				SetAsDefaultLogger();
		void			    lverbose( nxLogStatus status, const char * file, int line, const char *Format, ... );
		void				lrecord ( nxLogStatus status, const char * file, int line, const char *Format, ... );
		void				SetVerbose( nxBOOL on = nxTRUE) { m_isverbose = on; }
		void				CheckDefaultLogger();

	public:
static InxLog*				DefaultLog() { return m_DefaultLogger;}
static void					SetAsDefaultLogger(  InxLog* newdefault );
static void					Record ( nxLogStatus status, const char * file, int line, const char *Format, ... );
static void					Verbose( nxLogStatus status, const char * file, int line, const char *Format, ... );

	public:
        HRESULT STDMETHODCALLTYPE			QueryInterface( REFIID riid, void ** ppvObject);
        ULONG   STDMETHODCALLTYPE			AddRef();
        ULONG   STDMETHODCALLTYPE			Release();

		STDMETHOD(IsVerbose) (); // {if (m_isverbose) return S_OK; return S_FALSE;}
};

/*-----------------------------------------------------------------------------
 *					nxLogEntry		2004-11-23*/
/** \ingroup system_loggers
 *	This represents a single entry in the Logger. It is used internally only.
**/
/*---------------------------------------------------------------------------*/

class  nxLogEntry
{
	protected:
		nxTimeStamp		m_mjd;			// The current time.
		int				m_status;		// The status of the message (ERROR, WARNING, INFO)
		nxString		m_message;		// The message created by the user
		nxString		m_file;			// The source file wher the problem originated.
		int				m_line;			// The line number in the source file.

	public:
						nxLogEntry() { m_status = nxlog_info;}
		void			Set          ( int status, const char* mes, const char *file, int linenum);
		nxString		StatusString ( nxBOOL padstring = nxTRUE) const;			// returns the status as a string, pads it with white space so all statii are the same width in chars.
		const nxString&	Message()      const { return m_message;}
		const nxString& File()         const { return m_file;}
		int             LineNumber()   const { return m_line;}
		int				Status()       const { return m_status;}
		nxString		Mjd()	       const { return m_mjd.UTCStr();}
		bool			operator ==  (const nxLogEntry& ) const { return nxFALSE;}
		bool			operator !=  (const nxLogEntry& ) const { return nxFALSE;}
		bool			operator <   (const nxLogEntry& ) const { return nxFALSE;}
		bool			operator >   (const nxLogEntry& ) const { return nxFALSE;}
};

/*----------------------------------------------------------------------------
 *						class nxLog											*/
/**	\ingroup system_loggers
 * Logger devived form nxLogBase that implements thread safe behaviour and
 * enables streaming of log data to files.
**/
/*------------------------------------------------------------------------*/

class nxLog : public nxLogBase
{
	private:
		char						m_logName[MAX_PATH+1];							// The name of the log file, Dont use nxString for thread safety reasons.
		int							m_mode;											// The current log mode.
#if defined (NX_WINDOWS)															// on Windows systems
		HANDLE						m_mutex;										// provide mutex support for multiple threads
#endif																				// and that is that
	private:
		void						privateinit();
		void						WriteEntryToFile( const nxLogEntry& entry );		// write the entry to file.
		void						CreateMutex();										// Creates the mutex for locking purposes on a thread by thread basis
		void						Lock();												// Lock this thread so it uses the mutex
		void						Unlock();											// Unlock this thread
		void						CloseMutex();										// Close the mutex

	protected:
		virtual void   				DisplayEntry    ( const nxLogEntry& entry) = 0;		// display the log entry on a display device

	public:
									nxLog();
									nxLog( const char * logname );
								   ~nxLog();

	public:
		STDMETHOD(Trace)	    ( int status, const char * file, int line, const char *message);
		STDMETHOD(SetLogName)   ( const char * logname );
		STDMETHOD(GetLogName)	( const char** logname );
		STDMETHOD(SetLogMode)   ( int mode );
		STDMETHOD(ClearLogMode) ( int mode );
};


/*--------------------------------------------------------------------------
 *						class nxLogConsole									*/
/**	\ingroup system_loggers
 *	Implements streaming of logger messages to the console STDERR and/or STDOUT.
 *	This is the logger of choice for console based programs.
**/
/*------------------------------------------------------------------------*/

class  nxLogConsole : public nxLog
{
	private:
		nxBOOL				m_showUT;

	protected:
		virtual void		DisplayEntry( const nxLogEntry& entry);

	public:
							nxLogConsole();
		void				ShowUT( nxBOOL show );

};

#endif



