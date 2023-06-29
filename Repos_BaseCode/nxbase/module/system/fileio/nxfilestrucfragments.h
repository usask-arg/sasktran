
/*--------------------------------------------------------------------------
 *						class nxFileStrucFragments							*/
/**	\ingroup system_fileio
 *	This is a class that allows a user to read specific sections structures
 *	from a file of fixed length structures.  The idea is to just skip over
 *	the undesired sections and make the other sections look like a contiguous
 *	stream.  I use it to decode the OSIRIS housekkeping words from the ODIN SHK files
**/
/*--------------------------------------------------------------------------*/

class nxFileStrucFragments : public nxBinaryFile
{

	private:
		size_t				m_totalrecordsize;						// Size of 1 record as stored in the file
		size_t				m_cacherecordsize;						// Size of the users portion
		nxBYTE*				m_totalrecord;							// The total record holds one total record from the user.
		nxBYTE*				m_cache;								// The cache holds one "user record" of size m_cacherecordsize (ie unwanted fields have been stripped)
		nxBYTE*				m_cachepointer;							// Pointer to the current locationin the cache
		nxBYTE*				m_cacheendpointer;						// Pointer to the end of the cache.

	private:
		nxBOOL				UpdateCache();							// Loads a new record into the cache
		void				ReleaseCache();
		nxBYTE*				CacheCurrentPointer();					// Current location in cache.
		void				UpdateCachePointer( nxBYTE*	newptr );	// Update the cache pointer
		size_t				CacheNumBytesLeft();					// NUmber of bytes left in cache
		void				EmptyCacheContents();

	protected:
	virtual	nxBOOL			CopyRecordFieldsToCache( nxBYTE* rawrecord, size_t rawsize, nxBYTE* cacherecord, size_t usersize ) = 0;

	public:
							nxFileStrucFragments();
		virtual			   ~nxFileStrucFragments();
		nxBOOL				SetRecordSizes( size_t rawsize, size_t usersize );

	public:
		virtual	nxDWORD		FileSize();
		virtual nxBOOL		Open( const char * filename, const char * mode );
		virtual nxBOOL		IsOpen() const;
		virtual nxBOOL		Close( );
		virtual nxBOOL		eof();
		virtual size_t		Read  ( void *userbuffer, size_t tsize, size_t count);
		virtual size_t		Write ( const void *userbuffer, size_t size, size_t count);
		virtual int			Seek  ( size_t offset, int origin );
		virtual size_t		Tell  ();

};
