/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"
#include <boost/thread.hpp>


#if defined(NX_WINDOWS) || defined(_WIN32)
	#if defined(__BORLANDC__)
	  #include <dos.h>
	  #define MKDIR(x)              mkdir(x)
	#else
		#include <direct.h>
		#include <errno.h>
		#define MKDIR(x)              _mkdir(x)
	#endif

#else
//	#include <sys/mode.h>
	#include <sys/stat.h>
	#include <errno.h>									// added 2-NOV-1999
	#define MKDIR(x)              mkdir(x,0776)
	//#error nxdir.cxx need to implement UNIX VERSION of nxDirectory::CreateDirectory
#endif
//---------------------------------------------------------------------------
//						CreateDirectory
//---------------------------------------------------------------------------

nxBOOL nxDirectory::CreateADirectory( const char * dirname )
{
	int			status;
	nxBOOL		ok;
	nxString	pathname(dirname);
	static		boost::recursive_mutex	threadlock;

	if (pathname.GetLength() < 1)
	{
		ok = nxTRUE;
	}
	else
	{
		if (pathname[ pathname.GetLength()-1 ] != DIRECTORY_CHAR) pathname += DIRECTORY_CHAR;

		nxFileSpec		spec(pathname);
		nxStringArray	tokens;
		nxString		path;
		nxString		s;

		tokens.Strtok( spec.Directory(), "\\/" );					// This is threadsafe in our later code.

		path = spec.Drive();										// On windows systems include the drive
		if (!path.IsEmpty()) path += DIRECTORY_CHAR;				// and the first directory separator

		if (!spec.Directory().IsEmpty() && (spec.Directory().GetAt(0) == DIRECTORY_CHAR)) path += DIRECTORY_CHAR;	// On Unix include the first direcortu symbol if required
		ok = nxTRUE;
		boost::lock_guard< boost::recursive_mutex> lock(threadlock);	// Make sure only thread at a time is trying to create a directory structure
		for (int i=0; i < tokens.GetSize(); i++)
		{
			s = tokens[i];
			path += s;
			if (!nxDirectory::FileExists(path))
			{
				if ( !(s==".") && !(s==".."))
				{
					status = MKDIR( (const char*)path );
					ok = (status == 0);
					if (!ok)
					{
						int errcode = errno;
						ok = (errcode == EACCES) || (errcode == EEXIST);
					}
					if (!ok)
					{
						nxLog::Record( NXLOG_WARNING, "CreateDirectory, Error creating directory <%s>", (const char *)path );
						break;
					}
				}
			}
			path += DIRECTORY_CHAR;

		}
	}
	return ok;
}

//---------------------------------------------------------------------------------
//						nxDirectory::ScanDirectory
//	Scans the directory for the given filespec, which might include
//	wildcards.
//---------------------------------------------------------------------------------

int nxDirectory::ScanDirectory( const char * cfilespec, bool includedirectoriesasfiles, const char * directorytoscan )
{

	m_filespec = cfilespec;													// get the current filespec
	m_list.RemoveAll();														// clear the current list
	DirectoryScan( directorytoscan, *this, nxFALSE, m_filespec.DangerousTypecast(), includedirectoriesasfiles);		// scan the current directory and dont recurse
	return m_list.GetSize();												// return the number of points in the list.
}

void nxDirectory::operator() ( const char *fullfilename, nxBOOL isAdirectory )
{
	if (!isAdirectory) m_list.Add( fullfilename );				// add the name to the list of files
}

//---------------------------------------------------------------------------------
//						nxDirectory::FileExists
//	Static member that returns true if a file exists.
//---------------------------------------------------------------------------------


nxBOOL nxDirectory::FileExists( const char *pathname )
{
	struct stat buffer;

	return (stat( pathname, &buffer ) == 0);
}


/*---------------------------------------------------------------------------
 *'					nxFileLocator::AutoSearchDrivePartition        2001-12-19
 *	Scans through the directories looking for sub-directories of the form
 *	xxxxxxx\disk1   xxxxxx\disk2.  We use these to partion 
 *-------------------------------------------------------------------------*/

nxBOOL nxFileLocator::AutoSearchDrivePartition( const char *drivepartition, int startindex )
{
	STL(list)<nxString>							alist;
	typedef STL(list)<nxString>::iterator		iterator;

	iterator	iter;
	nxString	dirname;
	nxString	fullname;
	int			index;
	nxBOOL		ok;
	nxBOOL		newdirs = nxFALSE;
	int			i,n;

	n = m_paths.GetSize();													// get the number of directories in the path
	for (i = 0; i < n; i++ )												// Now check each path
	{																		// and for each path
		nxString& dirstr = m_paths.GetAt(i);								// get a reference to the path 
		alist.push_back(dirstr);												// still search this directory for the file
		index = startindex;
		do
		{
			dirname.sprintf("%s%1d", (const char*)drivepartition, (int)index);
			fullname = dirstr + dirname;
			ok = nxDirectory::FileExists(fullname);
			if (ok)
			{
				fullname += DIRECTORY_CHAR;
				alist.push_back( fullname );
				newdirs = nxTRUE;
				index++;
			}
		} while (ok);
	}
	if (newdirs)
	{
		m_paths.RemoveAll();
		for (iter = alist.begin(); !(iter == alist.end()); iter++)
		{
			m_paths.Add( *iter );
		}
	}
	return newdirs;
}

//---------------------------------------------------------------------------
//						nxFileLocator::FromString
//---------------------------------------------------------------------------

void nxFileLocator::FromString( const char *str )
{
	char	delimiter[2] = { DIRECTORY_SEARCH_PATH_DELIMITER, '\0'};
	int     n;
	int		l;
	int		i;

	m_paths.RemoveAll();														// remove all of the entries in the string
	if (str != NULL )															// if the string is valid
	{																			// then
		m_paths.Strtok(str, delimiter);											// find all of the paths by searching for the delimiters
		n = m_paths.GetSize();													// get the number of directories in the path
		for (i = 0; i < n; i++ )												// Now check each path
		{																		// and for each path
			nxString& dirstr = m_paths.GetAt(i);								// get a reference to the path 
			l = dirstr.GetLength();												// get the length of the string
			if (dirstr[l-1] != DIRECTORY_CHAR) dirstr += DIRECTORY_CHAR;		// make sure we have a directory char appended to the end of the string
		}																		// do all of the elements.
	}
}


//---------------------------------------------------------------------------
//						nxFileLocator::FromEnvironmentVar
//---------------------------------------------------------------------------

void nxFileLocator::FromEnvironmentVar( const char *envar )
{
	char*	str;

	if (envar != NULL) str = getenv(envar);
	else               str = NULL;
	FromString( str );
}

//---------------------------------------------------------------------------
//						nxFileLocator::FindFile
//---------------------------------------------------------------------------

nxBOOL nxFileLocator::FindFile( const char*filename, nxString* finalfullname )
{
	int	n;
	int i;
	nxString	basename(filename);
	nxString	fullname;
	nxBOOL		ok = nxFALSE;

	n = m_paths.GetSize();
	for (i=0; (i < n) && (!ok); i++)
	{
		fullname = m_paths[i] + basename;
		ok = nxDirectory::FileExists( fullname );
	}
	if (ok) *finalfullname = fullname;
	else    *finalfullname = "";
	return ok;
}

																	 
