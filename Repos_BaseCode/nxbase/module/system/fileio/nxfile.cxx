
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"

//---------------------------------------------------------------------------
//	Note that strict ANSI does not include the followig functions
//
//	stat
//	fileno
//	fstat
//---------------------------------------------------------------------------

#if !defined(SUPPORT_ANSI_FUNCTIONS_ONLY)

nxDWORD nxFile::FileSize( const char *path )
{
	struct stat	buffer;
	nxDWORD			size;

	if ( stat( path, &buffer ) == 0 )
	{
		size = (nxDWORD)buffer.st_size;
	}
	else size = 0;
	return size;
}

nxBOOL nxFile::Exists ( const char *path )
{
	struct stat	buffer;
	return stat( path, &buffer ) == 0;
}

int nxFile::Handle()
{
	int hdl;

	if (IsOpen())
	{
#if defined(NX_WINDOWS)
		hdl =  _fileno(m_file);
#else
		hdl =   fileno(m_file);
#endif
	}
	else
	{
		hdl = -1;
	}
	return hdl;
}

nxDWORD nxFile::FileSize()
{
	struct stat	buffer;
	nxDWORD			size;

	
	if ( IsOpen() && fstat( Handle(), &buffer ) == 0 )
	{
		size = (nxDWORD)buffer.st_size;
	}
	else size = 0;
	return size;
}

#endif		// End of SUPPORT_ANSI_FUNCTIONS ONLY

//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------

nxFile::nxFile( FILE* file )
{
	m_file = file;
	m_logerrors = nxTRUE;
	m_filename = "<Unknown>";
}

//---------------------------------------------------------------------------
//+
//NAME:
//					nxFile::Constructor
//-
//---------------------------------------------------------------------------

nxFile::nxFile()
{   
	m_file = NULL;
	m_logerrors = nxTRUE;
	m_filename.Empty(nxTRUE);
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxFile::Destructor
//-
//---------------------------------------------------------------------------


nxFile::~nxFile()
{
	Close();
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxFile::Close
//-
//---------------------------------------------------------------------------


nxBOOL nxFile::Close()
{
	nxBOOL	idx;
	
	if (m_file != NULL)
	{
		idx  = (fclose(m_file)==0);
		m_file = NULL;
		m_filename.Empty(nxTRUE);
	}
	else idx = nxTRUE;
	return idx;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxFile::Open
//-
//---------------------------------------------------------------------------


nxBOOL	nxFile::Open( const char * filename, const char * mode )
{
	if (m_file != NULL)
	{
		if (m_logerrors)
		{
			nxLog::Record( NXLOG_WARNING,"nxFile::Open, closing file %s before opening %s ", (const char *)m_filename, filename);
		}
		Close();
	}
	m_filename = filename;
	if (!m_filename.IsEmpty() )		// Empty filenames cause a crash in Visual C++
	{
		m_file = fopen( filename, mode );
	}
	else
	{
		m_file = NULL;
	}
	return (m_file != NULL);
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxFile::ReadString
//-
//---------------------------------------------------------------------------

nxBOOL nxFile::ReadString( nxString & str)
{
	nxBOOL	status;
	char	buffer[1027];
	int		n;
	
	if (m_file == NULL)
	{
		nxLog::Record( NXLOG_ERROR,"nxFile::ReadString, attempting to read a string from a NULL file handle");
		str.Empty(nxFALSE);
		status = nxFALSE;
	}
	else
	{
		status = (fgets( buffer, sizeof(buffer), m_file ) != NULL);
		if (status)
		{
			n = (int)strlen(buffer);
			if (n > 0 && buffer[n-1] == '\n') buffer[n-1] = '\0';
			str = buffer;
		}
		else
		{
			str.Empty(nxFALSE);
		}
	}
	return status;
}
//---------------------------------------------------------------------------
//+
//NAME:
//						nxFile::WriteString
//-
//---------------------------------------------------------------------------

nxBOOL nxFile::WriteString( const char *str)
{
	nxBOOL	status;
	
	if (m_file == NULL)
	{
		nxLog::Record( NXLOG_WARNING, "nxFile::WriteString, attempting to write a string to a NULL file handle [%s]", (const char *)m_filename);
		status = nxFALSE;
	}
	else
	{
		if (str == NULL) 
		{
			status = nxFALSE;
		}
		else
		{
			status = (fprintf(m_file, "%s", str) >= 0);
		}
	}
	return status;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxFile::operator >>
//-
//---------------------------------------------------------------------------


//nxFile&		nxFile::operator >> ( nxString& str )
//{
//	ReadString(str);
//	return *this;
//}
//---------------------------------------------------------------------------
//+
//NAME:
//						nxFile::operator <<
//-
//---------------------------------------------------------------------------


//nxFile&		nxFile::operator << ( nxString& str )
//{
//	WriteString(str);
//	return *this;
//}

//nxFile& nxFile::operator << (const char * str )
//{
//	WriteString(str);
//	return *this;
//}

nxBOOL nxFile::eof() 
{
	nxBOOL	eof_flag;

	if (m_file != NULL)
	{
		eof_flag =  feof(m_file) != 0;
	}
	else
	{
		eof_flag = nxTRUE;
	}
	return eof_flag;
}

nxString nxFile::ReadALine()
{
	nxString	theanswer;
	char*	    buffer;
	nxBOOL	    ok;
	int         n = 1024;
	bool		more;

	if (m_file == NULL)
	{
		nxLog::Record( NXLOG_ERROR, "nxFile::ReadAline, trying to read a line from file [%s] with a NULL file handle", (const char *)m_filename);
		return "";
	}

	buffer = new char[n+1];
	if (buffer == NULL)NXTHROW(("Memory allocation in nxFile::ReadAline\n"));
	ok = !eof();
	more = true;
	while (ok && more)
	{
		ok = (fgets(buffer, n, m_file ) != NULL);
		if (ok)
		{
			buffer[n] = '\0';
			n = (int)strlen(buffer);
			more = (buffer[n-1] != '\n');
			if (!more) buffer[n-1] = '\0';
			theanswer += buffer;
		}
	}
	delete [] buffer;
	return theanswer;
}

  




