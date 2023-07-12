
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/


/*--------------------------------------------------------------------------
 *					class nxDirectory										*/
/**	\ingroup system_fileio
 *	Perform directory operations.  Typically used to get of files in a
 *	directory.
**/
/*-------------------------------------------------------------------------*/

class  nxDirectory
{
	private:
		nxStringArray	m_list;				// a current list of files.
		nxString		m_filespec;


	public:
		static nxBOOL	FileExists( const char *pathname );
		static nxBOOL	CreateADirectory( const char * dirname );
		int				ScanDirectory( const char * filespec,  bool includedirectoriesasfiles= false, const char * directorytoscan="" );
						nxDirectory(){}
						nxDirectory( const char *filespec ){ ScanDirectory(filespec);}
		nxStringArray&	List() { return m_list;}
		void			operator()( const char* fullfilename, nxBOOL isAdirectory );

};

/*-------------------------------------------------------------------------------
 *					class nxFileSpec											*/
/**	\ingroup system_fileio
 *	Class used to parse filenames into constituent pieces.
**/
/*------------------------------------------------------------------------------*/

class  nxFileSpec
{
	private:
		nxString	m_extension;
		nxString	m_name;
		nxString	m_directory;
		nxString	m_drive;

    public:
					nxFileSpec(){};
					nxFileSpec( const char *filnam);
		const char*	operator= ( const char *filnam);
		nxString	Extension()                   {return m_extension;}
		nxString    Name()                        {return m_name;}
		nxString	Directory()                   {return m_directory;}
		nxString	Drive()						  {return m_drive;}
		nxString	FullDirSpec()				  {return (m_drive + m_directory);}
		nxString	FullPathName()				  {return (m_drive + m_directory + m_name + m_extension);}
		void		SetExtension(const char *ext) { m_extension = ext;}
		void		SetName     (const char *nme) { m_name      = nme;}
		void		SetDirectory(const char *dir) { m_directory = dir;}
		void		SetDrive    (const char *drv) { m_drive     = drv;}
};

/*----------------------------------------------------------------------------
 *						class nxFileLocator									*/
/**	\ingroup system_fileio
 *	Used to locate files relative to an array of search path directories.
**/
/*---------------------------------------------------------------------------*/

class nxFileLocator
{
	private:
		nxStringArray	m_paths;


	public:
						nxFileLocator(){};
						nxFileLocator( const char* environmentvar ) { FromEnvironmentVar(environmentvar);}
		nxBOOL			AutoSearchDrivePartition( const char *drivepartition, int startindex );
		void			FromString( const char * str );
		void			FromEnvironmentVar( const char* environmentvar);
		nxBOOL			FindFile( const char* filename, nxString* fullname );
};


