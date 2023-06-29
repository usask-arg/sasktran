#include "nxbase_core.h"
//#include "nxfilesessionname.h"

//---------------------------------------------------------------------------
//						nxFileSesionName::nxFileSessionName
//---------------------------------------------------------------------------

nxFileSessionName::nxFileSessionName()
{
	m_filedate      = -1.0;
	m_fileindex     = -1;
	m_synclog       = nxFALSE;
	m_onefileperday = nxFALSE;
}

//---------------------------------------------------------------------------
//						nxFileSessionName::Filename
//---------------------------------------------------------------------------

nxBOOL nxFileSessionName::CheckSessionName( nxBOOL forcenewsession )
{
	nxTimeStamp		ut;
	nxBOOL			newfile;
	nxString		logname;

	forcenewsession = forcenewsession && !m_onefileperday;		// Cant force a new session if only one file per day

	ut.FromSystem();											// Get the current time from the system
	ut = ut.ZeroUT();											// and get the day number
	newfile = (ut != m_filedate);								// see if we have changed date
	if (newfile)												// if we have a new day the we must definite
	{															// then
		m_filedate = ut;										// Update the new date
		ResetSessionIndex();									// scan the directory for the last session on this date.
	}
	newfile |= forcenewsession;									// Do we have to force a new file
	if (newfile)												// we do
	{															// then
		nxString name;											// a local variable to the name of the file

		GetBaseName( &name );									// Get the base name of the next file
		if (m_onefileperday)									// if we are just using one file per day
		{														// then
			m_fullname = name;									// use just the basename
		}														// and that is that
		else													// otherwise we are using sessions
		{														// so
			m_fileindex++;										// point to the next session and generate filename without extension
			m_fullname.sprintf( "%s%04d", (const char *)name, (int)m_fileindex);
		}														// and that is that

		if (m_synclog)											// if we are syncing logs
		{														// then
			InxLog* log = nxLogBase::DefaultLog();				// Get the current default log
			if (log != NULL)									// if the log is not NULL
			{													// then
				logname     = m_fullname + ".log";				// Generate the associated log name
				log->SetLogName( logname );					// and set the log name
			}
		}
		m_fullname += m_extension;
	}
	return newfile;
}


//---------------------------------------------------------------------------
//						nxFileSessionName::CoreName
//---------------------------------------------------------------------------

void nxFileSessionName::GetCoreName( nxString* corename )
{
	int				year, month,day, hour, mins, secs;
	double			tick;

	m_filedate.GetUTC( &day, &month, &year, &hour, &mins, &secs, &tick );
	year = year%100;
	corename->sprintf( "%s%02d%02d%02d", (const char *)m_experimentname, (int)year, (int)month, (int)day );
}


//---------------------------------------------------------------------------
//						nxFileSessionName::GetBaseName
//---------------------------------------------------------------------------

void nxFileSessionName::GetBaseName( nxString* basename )
{
	nxString	corename;
	GetCoreName(&corename );
	basename->sprintf( "%s%s", (const char *)m_outputdir, (const char *)corename );
}

//---------------------------------------------------------------------------
//						nxFileSessionName::ResetSessionIndex
//---------------------------------------------------------------------------

void nxFileSessionName::ResetSessionIndex()
{
	nxDirectory		dir;
	nxString		wildcard;
	nxString		name;
	nxString		corename;
	nxFileSpec		spec;
	int				index;
	int				idx;
	int				clen;

	m_fileindex = -1;													// reset the file index
	if (!m_onefileperday)												// If we are alloweed multiple files per day
	{																	// then
		GetCoreName( &corename );										// get the corename
		clen = corename.GetLength();
		GetBaseName( &wildcard );										// Get the base name
		wildcard += "*" + m_extension;									// add on the wildcard
		dir.ScanDirectory( wildcard );									// now scan the directory for this session
		index       = -1;

		nxStringArray&	list = dir.List();								// Now get the list of object
		for (int i =0; i < list.GetSize(); i++)							// for each object
		{																// then
			spec = list[i];												// get the filename and get the file specs
			name = spec.Name();											// get just the name
			if (name.GetLength() >= (clen+4))							// and if its reasonable
			{															// then
				idx = name.Find(corename);								// See if we can see the corename
				if (idx == 0)											// and if its at the beginning
				{														// then
					index = atoi( name.Right(name.GetLength()- clen ));	// get the characters to right of the corename as an index;
					if (index > m_fileindex) m_fileindex = index;		// and update the index
				}
			}
		}
	}
}

//---------------------------------------------------------------------------
//						nxFileSessionName::Initialize
//---------------------------------------------------------------------------

void nxFileSessionName::Initialize( const char * basedir, const char* experimentname, const char *extension, nxBOOL onefileperday, nxBOOL synchronizelogs)
{
	m_filedate       = -1.0;															// force an update at the next request
	m_onefileperday  = onefileperday;													// save if we want just one file per day
	m_synclog        = synchronizelogs;
	m_experimentname = experimentname;
	m_outputdir      = basedir;															// Copy over the base directory
	m_fileindex      = -1;

	if (!m_outputdir.IsEmpty())															// if we have a non-empty dir char
	{																					// then
		if ( m_outputdir[ m_outputdir.GetLength()-1 ] == DIRECTORY_CHAR)				// if the last char is the not a directory char
		{																				// then
			m_outputdir = m_outputdir.Left( m_outputdir.GetLength()-1);					// add one
		}																				// otherwise we'll get bad filenames later on
		nxDirectory::CreateADirectory( m_outputdir );
		m_outputdir += DIRECTORY_CHAR;
	}

	m_extension = extension;															// copy over the extension
	if (!m_extension.IsEmpty())															// if we have an extension string
	{																					// and the first char is not a dot
		if ( m_extension[0] != '.') m_extension = '.' + m_extension;					// the add a dot at the beginning of the char.
	}

}



