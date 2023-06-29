/*--------------------------------------------------------------------------
 *						class nxFileSessionName								*/
/**	\ingroup system_fileio
 *	This is class used for generating filenames when archiving experimental
 *	data.   The basic idea is to automatically generate a unique name as the
 *	experiment continues. *	The basic paradigm is store "sessions" to a file.  The basic form
 *	of the file is..
 *
 *	[base dirspec][maincategory][Data][session Number].[extension]
 *	All files are written to the [base dir spec]
**/
/*--------------------------------------------------------------------------*/

class nxFileSessionName
{
	private:
		nxString			m_outputdir;						// Output directory for the file
		nxString			m_experimentname;					// The experiment name typically 0,1, 2 or 3 characters
		nxString			m_extension;						// the file extension (includes leasding '.')
		nxString			m_fullname;							// The current full filename
		int					m_fileindex;						// the current full fileindex
		nxTimeStamp			m_filedate;							// the current file date
		nxBOOL				m_onefileperday;					// Flags that we should only use one file per day (aka no session number)
		nxBOOL				m_synclog;							// Flags that we should sync the log with this exeperiment

	private:
		void				ResetSessionIndex();
		void				GetBaseName( nxString* basename );
		void				GetCoreName( nxString* corename );

	public:
							nxFileSessionName( );
		void				Initialize( const char * basedir, const char * experimentname, const char *extension, nxBOOL onefileperday, nxBOOL synchronizelogs);
		nxBOOL				CheckSessionName( nxBOOL forcenewsession = nxFALSE );
		void				UpdateSessionNumber()            { CheckSessionName(nxTRUE);}
		nxString			GetCurrentFilename()		     { return m_fullname;}
		int					GetCurrentSessionNumber()        { return m_fileindex;}
};


