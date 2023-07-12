
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/



#if defined(NX_WINDOWS) 


#if defined(__BORLANDC__)
	#include <dos.h>
	#define CHDIR(x)              chdir(x)
	#define FINDDATA	 			ffblk
	#define FINDHDL    			intptr_t
	#define FINDFIRST( path, blk)	findfirst(path, blk, (FA_RDONLY | FA_HIDDEN | FA_SYSTEM | FA_DIREC | FA_ARCH))
	#define FINDNEXT( hdl, blk  ) findnext( blk )
	#define FINDCLOSE(hdl)
	#define ATTRIB_TYPE			unsigned long
	#define ATTRIB_FIELD(x)  		x.ff_attrib
	#define NAME_FIELD(x)			x.ff_name
	#define ATTRIB_HIDDEN			FA_HIDDEN
	#define ATTRIB_SYSTEM 		FA_SYSTEM
	#define ATTRIB_DIR    		FA_DIREC
#else
	#if defined(_WIN32)
		#define CHDIR(x)				_chdir(x)
		#define FINDDATA				_finddata_t
		#define FINDHDL    				intptr_t
		#define FINDFIRST( path, blk)	_findfirst(path, blk )
		#define FINDNEXT(hdl, blk)      _findnext( hdl, blk)
		#define FINDCLOSE(hdl)          _findclose(hdl)
		#define ATTRIB_TYPE				unsigned
		#define ATTRIB_FIELD(x)  		x.attrib
		#define NAME_FIELD(x)			x.name
		#define ATTRIB_HIDDEN			_A_HIDDEN
		#define ATTRIB_SYSTEM 			_A_SYSTEM
		#define ATTRIB_DIR    			_A_SUBDIR
	#elif defined (_WIN64)
		#define CHDIR(x)				_chdir(x)
		#define FINDDATA				_finddata64_t
		#define FINDHDL    				intptr_t
		#define FINDFIRST( path, blk)	_findfirst64(path, blk )
		#define FINDNEXT(hdl, blk)      _findnext64( hdl, blk)
		#define FINDCLOSE(hdl)          _findclose(hdl)
		#define ATTRIB_TYPE				unsigned
		#define ATTRIB_FIELD(x)  		x.attrib
		#define NAME_FIELD(x)			x.name
		#define ATTRIB_HIDDEN			_A_HIDDEN
		#define ATTRIB_SYSTEM 			_A_SYSTEM
		#define ATTRIB_DIR    			_A_SUBDIR
#endif

#endif


/*--------------------------------------------------------------------------
 *					DirectoryScan<ANALYSEFILEOBJECT>						*/
/**	\ingroup system_fileio
 *	Scans the directory tree starting at this directory and calls the function
 *	object to execute a specific action on each file. Typical usage is to
 *	quickly scan through the directory tree. The subroutine does not look at
 *	HIDDEN or SYSTEM files.
 *
 *	\param fulldirname
 *		The full path name of the directory to scan
 *
 *	\param action
 *		The function object to call for each object of each file ANALYSEFILEOBJECT&
 *
 *	\param recurse
 *		If nxTRUE then recurse through the sub-directories. The action
 *		object will still receive a call for each directory entry at this level
 *
 *`	\param wildcard
 *		The wild card specification to be used for selecting files.
 *
 * \par ANALYSEFILEOBJECT;
 *	This is a templated A function object. A function of the form:
 *	\code
 *	void ANALYSEFILEOBJECT(const char* fullfilename, nxBOOL isAdirectory)
 *	\endcode
 *	or a class that has the operator
 *	\code
 *	void	ANALYSEFILEOBJECT::operator() (const char* fullfilename, nxBOOL isAdirectory )
 *	\endcode
 *
 *	\par HISTORY
 *	2002-9-30\
**/
/*------------------------------------------------------------------------*/

template <class ANALYSEFILEOBJECT>
void DirectoryScan(
				   const char *			fulldirname,
				   ANALYSEFILEOBJECT&	action,
				   nxBOOL				recurse,
				   char*				wildcard,
				   bool					adddirectories = false)

{
	FINDDATA			DirStruct;
	FINDHDL				handle;
	int					result;
	nxString			name;
	ATTRIB_TYPE			attrib;
	nxBOOL				ignore;
	char				olddir[1100];
	STL(list<nxString>)	alist;

	int thisdrive = _getdrive();												// Get the current drive
	_getdcwd( thisdrive, olddir, sizeof(olddir) );								// Get the current working directory on this drive.

	nxString	basename = fulldirname;
	int nc = basename.GetLength();
	if (nc > 0)
	{
		if (basename[nc-1] != '\\')
		{
 			CHDIR( basename );													// Change to the required directory
			basename += "\\";
		}
		else
		{
			nxString dummy;
			if (basename[nc-2] != ':')
			{
				dummy = basename.Left(nc-1);
			}
			else
			{
				dummy = basename;
			}
			CHDIR( dummy );
		}
	}

	handle = FINDFIRST( wildcard, &DirStruct);									// get the first file in the list
	if (handle != -1)															// if we found nothing then dont bother
	{																			// otherwise
		result = 0;																// set the result equal to zero
		while (result == 0 )													// while we have good files
		{																		// then
			attrib = ATTRIB_FIELD(DirStruct);									// check out the file attributes
			ignore =     ((attrib & ATTRIB_HIDDEN) != 0)						// skip over hidden files.
					  || ((attrib & ATTRIB_SYSTEM) != 0);						// and system files
			if (!ignore)														// shall we ignore this file
			{																	// nope
				name = NAME_FIELD(DirStruct);									// Get the name of this file (without path info)
				if ((attrib & ATTRIB_DIR) != 0)									// is this a sub directory
				{																// yes
					ignore = (name == ".") || (name == "..");					// see if its special dir info
					if (!ignore)												// it isnt
					{															// so
						nxString	newdir = basename + name;
						alist.push_back( newdir );								// and put it in our list
						if (adddirectories) action(newdir, nxFALSE);
					}
				}																// otherwise
				else															// file is not a directory
				{																// so
					char scandir[256];
					int  scandrive = _getdrive();								// Get the current scan drive
					_getdcwd( scandrive, scandir, sizeof(scandir) );			// Get the current scan directory on this drive.
					CHDIR( olddir );											// change back to the users call directory
					nxString newname = basename + name;
					action(newname, nxFALSE);
					CHDIR( scandir );											// change back to the scan directory
				}
			}
			result = FINDNEXT( handle, &DirStruct );							// find the next file in the list
		}																		// Do all the files
		FINDCLOSE(handle);														// now close the file handle
	}

	CHDIR( olddir );																			// change back to the old directory
	for ( STL(list<nxString>::iterator) ptr = alist.begin(); !(ptr == alist.end()); ptr++)		// Now for scan through the list of sub-directories
	{																							// and
		action( (const char *)(*ptr), nxTRUE);
		if (recurse) DirectoryScan( (const char *)(*ptr), action,  recurse, wildcard );
	}
}



#else

//---------------------------------------------------------------------------
//					Unix Version
//---------------------------------------------------------------------------

template <class ANALYSEFILEOBJECT>
void DirectoryScan( const char * fulldirname,  ANALYSEFILEOBJECT& action,  nxBOOL recurse, char *wildcard, bool adddirectories = false)
{
	nxWildcard			wspec( wildcard );
	STL(list<nxString>)	alist;
	DIR*				handle;
	struct dirent*		Dirstruct;
	nxString			name;
	nxBOOL				ignore;
	nxString			fullfilename;
 	struct stat			attrib;
	nxString			basename = fulldirname;							// basename hold directory name with trailing "/"
	nxString			dirspec;										//dirspec hold directory name with no trailing "/"

	int nc = basename.GetLength();										// get the length of the directory entry
	if (nc > 0)															// if we have a string
	{																	// then
		if (basename[nc-1] == '/')										// if it is terminated with a final slash
		{																// then remove the
			dirspec = basename.Left(nc-1);								// slash for the directort spec
		}																// and that is that
		else															// otherwise there is no terminating slash
		{																// so
			dirspec = basename;											// save this as the dirspec variable
			basename += "/";											// and add the slash to the directory spec
		}
	}
	else																// if basename is empty
	{																	// then the
		dirspec = ".";													// directory spec is the current one
	}

	handle = opendir( dirspec );										// open the specified directory specification.
	if (handle == NULL) return;											// if it failed then simply return
	do																	// now scan through the files in this directory
	{																	// so
		Dirstruct = readdir(handle);									// get an entry from the directory
		if (Dirstruct != NULL)											// see if we succeeded
		{																// then
			name = Dirstruct->d_name;									// Get the name of this file (without path info)
			ignore = !wspec.Match(name) || (name == ".") || (name == "..");				// see if its special dir info'
			if (!ignore)												// if it isnt one of those two
			{															// then
				fullfilename = basename + name;							// get the fullname of this entry
				ignore = stat(fullfilename, &attrib) != 0;				// and get the "stats" on this file
				if (!ignore)											// if we got the stats ok
				{														// then
					if (S_ISDIR(attrib.st_mode))						// if this is a directory
					{													// then
						alist.push_back( fullfilename);					// push it onn the list for later use.
						if (adddirectories) action(fullfilename, nxFALSE);

					}													// and that finished directory entries.
					else												// no do all the other possibilities
					{
						if (S_ISREG(attrib.st_mode))					// so if this is a regular file
						{												// then
							action(fullfilename, nxFALSE);				// perform the action on this file
						}
					}
				}
			}
		}																// Do all the files
	} while (Dirstruct != NULL);										// unless the user requests an abort
	closedir(handle);													// now close the file handle

	for ( STL(list<nxString>::iterator) ptr = alist.begin(); !(ptr == alist.end()); ptr++)		// Now for scan through the list
	{																							// and
		action( (const char *)(*ptr), nxTRUE);
		if (recurse) DirectoryScan( (const char *)(*ptr), action,  recurse, wildcard );
	}
}
#endif


