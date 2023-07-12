#if !defined(NXBASE_NXFILE_H)
#define NXBASE_NXFILE_H 1

#include "../../system/strings/nxstring.h"

/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/


/*----------------------------------------------------------------------------
 *					class nxFile											*/
/**	\ingroup system_fileio
 *	A Wrapper class for the basic FILE* handler.  Also see nxBinaryFile
**/
/*------------------------------------------------------------------------*/

class  nxFile
{
	protected:
		nxString	m_filename;
		nxBOOL		m_logerrors;
		FILE*		m_file;


	public:
					nxFile();
					nxFile( FILE* file );
		virtual	   ~nxFile();

		const nxString&	Filename() const {return m_filename;}
		nxString		ReadALine();
		nxBOOL			ReadString( nxString& str);
		nxBOOL			WriteString( const char * str );
		inline			operator FILE* () { return m_file;}
//		nxFile&			operator >> ( nxString& str );
//		nxFile&			operator << ( nxString& str );
//		nxFile&			operator << ( const char * str );


#if !defined( SUPPORT_ANSI_FUNCTIONS_ONLY )
static  nxDWORD		FileSize( const char *path );
static 	nxBOOL		Exists( const char *path);
		virtual	 nxDWORD		FileSize();
		int			Handle();
#endif

	public:
		virtual nxBOOL		Open( const char * filename, const char * mode );
		virtual nxBOOL	    IsOpen ()const {return m_file != NULL;}
		virtual nxBOOL		Close();
		virtual nxBOOL		eof();
		virtual size_t		Read  ( void *buffer, size_t size, size_t count )       {return ::fread ( buffer, size, count, m_file);}
		virtual size_t		Write ( const void *buffer, size_t size, size_t count)	{return ::fwrite( buffer, size, count, m_file);}
		virtual int			Seek  ( size_t offset, int origin )    					{return ::fseek ( m_file, (long)offset, origin);}
		virtual size_t		Tell  ()												{return ::ftell ( m_file );}
};



/*----------------------------------------------------------------------------
 *					class nxBinaryFile										*/
/**	\ingroup system_fileio
 *	class used for binary I/O operations. In particular it overloads I/O for
 *	many of the basic data types.
**/
/*------------------------------------------------------------------------*/

class  nxBinaryFile : public nxFile
{
	private:
		nxBOOL			m_ok;

	public:
		nxBinaryFile( FILE* file) : nxFile( file ) { m_ok = nxTRUE; }
						nxBinaryFile() { m_ok = nxTRUE;}
						operator FILE* () { return m_file;}
		nxBinaryFile&	operator<<( const char * str );
		nxBinaryFile&	operator<<( const double val );
		nxBinaryFile&	operator<<( const float  val );
		nxBinaryFile&	operator<<( const nxBYTE val );
		nxBinaryFile&	operator<<( const nxWORD val );
		nxBinaryFile&	operator<<( const nxDWORD val );
		nxBinaryFile&	operator<<( const int val );
		nxBinaryFile&	operator>>( nxString& str );
		nxBinaryFile&	operator>>( double&  val );
		nxBinaryFile&	operator>>( nxBYTE&  val );
		nxBinaryFile&	operator>>( nxWORD&  val );
		nxBinaryFile&	operator>>( nxDWORD& val );
		nxBinaryFile&	operator>>( int& val );
		nxBinaryFile&	operator>>( float& val );
		nxBinaryFile&	operator<<( const nxBOOL bval );
		nxBinaryFile&	operator>>( nxBOOL& bval );


		nxBOOL			IsOk() const {return m_ok;}
};

#endif

