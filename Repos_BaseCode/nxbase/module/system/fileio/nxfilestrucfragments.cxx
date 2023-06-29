#include "nxbase_core.h"

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::nxFileStrucFragments
 *-------------------------------------------------------------------------*/

nxFileStrucFragments::nxFileStrucFragments()
{
	m_totalrecordsize = 0;						// Size of 1 record as stored in the file
	m_cacherecordsize = 0;						// Size of the users portion
	m_totalrecord     = NULL;					// The total record holds one total record from the user.		
	m_cache           = NULL;					// The cache holds one "user record" of size m_cacherecordsize (ie unwanted fields have been stripped)
	m_cachepointer    = NULL;					// Pointer to the current locationin the cache
	m_cacheendpointer = NULL;					// Pointer to the end of the cache.
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::ReleaseCache
 *-------------------------------------------------------------------------*/

void nxFileStrucFragments::ReleaseCache()
{
	if (m_totalrecord != NULL) delete [] m_totalrecord;
	if (m_cache       != NULL) delete [] m_cache;
	m_totalrecordsize = 0;						// Size of 1 record as stored in the file
	m_cacherecordsize = 0;						// Size of the users portion
	m_totalrecord     = NULL;					// The total record holds one total record from the user.		
	m_cache           = NULL;					// The cache holds one "user record" of size m_cacherecordsize (ie unwanted fields have been stripped)
	m_cachepointer    = NULL;					// Pointer to the current locationin the cache
	m_cacheendpointer = NULL;					// Pointer to the end of the cache.
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::~nxFileStrucFragments
 *-------------------------------------------------------------------------*/

nxFileStrucFragments::~nxFileStrucFragments()
{
	ReleaseCache();
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::SetRecordSizes
 *-------------------------------------------------------------------------*/

nxBOOL nxFileStrucFragments::SetRecordSizes( size_t rawsize, size_t usersize )
{
	nxBOOL ok;

	ReleaseCache();													// Release any existing cache.
	ok = (rawsize > 0) && (usersize > 0);							// make sure we hae sensible buffer sizes
	if (!ok)														// If it doesnt then tell the user
	{																// about it
		nxLog::Record(NXLOG_WARNING, "nxFileStrucFragments::SetRecordSizes, record sizes must be greater than 0 ");
	}																// otherwise
	else															// params are ok
	{																// so do the allocation
		m_totalrecordsize = rawsize;								// copy over the size of the raw record stored on file
		m_cacherecordsize = usersize;								// copy over the size of the user buffer
		m_totalrecord     = new nxBYTE[rawsize];					// allocate space for the raw record as stored on file
		m_cache           = new nxBYTE[usersize];					// allocate memory for the user buffer
		ok = (m_totalrecord != NULL)  && (m_cache != NULL);			// Make sure the memory allocation worked
		if (!ok)													// if it did not then
		{															// tell the user about it
			nxLog::Record(NXLOG_WARNING, "nxFileStrucFragments::SetRecordSizes, Error allocating memory");
			ReleaseCache();											// and release the resources.
		}															// and that is that
		else														// otehrwise everything is OK
		{															// so
			m_cacheendpointer = m_cache + usersize;					// get the end of the user cache
			m_cachepointer    = m_cacheendpointer;					// and setup cache so it is currently empty
		}
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::UpdateCache
 *-------------------------------------------------------------------------*/

nxBOOL nxFileStrucFragments::UpdateCache()
{
	nxBOOL ok;

	ok = (m_cache != NULL) && (	m_totalrecord != NULL);							// Make sure we have memory allocated
	if (!ok)																	// if we dont then
	{																			// log a message
		nxLog::Record(NXLOG_WARNING, "nxFileStrucFragments::UpdateCache, The cahce is empty, there is nothing to do");
	}																			// and that is that
	else																		// otherwise
	{																			// Read in the next record
		ok = (nxBinaryFile::Read( m_totalrecord, m_totalrecordsize, 1) == 1);	// Fetch the record and make sure it is ok.
		if (ok)																	// if it was okay
		{																		// then ask derived classes
			ok = CopyRecordFieldsToCache( m_totalrecord,						// to process this record
										  m_totalrecordsize,
										  m_cache,
										  m_cacherecordsize );
			m_cachepointer = m_cache;										// and updatethe cache pointer	
		}																		// and that 
	}																			// is that
	if (!ok) m_cachepointer = m_cacheendpointer;								// if there was an error then empty the cache
	return ok;																	// and return status to the caller.
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::CacheCurrentPointer
 *-------------------------------------------------------------------------*/

nxBYTE* nxFileStrucFragments::CacheCurrentPointer()
{
	return m_cachepointer;
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::UpdateCachePointer
 *-------------------------------------------------------------------------*/

void nxFileStrucFragments::UpdateCachePointer( nxBYTE* ptr)
{
	m_cachepointer = ptr;
	NXASSERT( m_cachepointer <= m_cacheendpointer);
	NXASSERT( m_cachepointer >= m_cache);
}


/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::CacheNumBytesLeft
 *-------------------------------------------------------------------------*/

size_t nxFileStrucFragments::CacheNumBytesLeft()					// NUmber of bytes left in cache
{
	return (m_cacheendpointer - m_cachepointer);
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::FileSize
 *	Returns the size of the file
 *-------------------------------------------------------------------------*/

nxDWORD nxFileStrucFragments::FileSize()
{
	nxDWORD nbytes  = nxBinaryFile::FileSize();
	return (nxDWORD)((nbytes/m_totalrecordsize)*m_cacherecordsize);
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::Open
 *-------------------------------------------------------------------------*/

nxBOOL nxFileStrucFragments::Open( const char * filename, const char * /*mode*/ )
{
	nxBOOL ok;
	nxString	rb = "rb";

	rb.MakeLower();											// Make the string lower case
	ok = (rb.Find('r') >= 0);								// Look for the "r" for read only
	if (!ok)												// if we did not find it
	{
		nxLog::Record( NXLOG_WARNING, "nxFileStrucFragments::Open, Currently we only support reading binary files. Please use a readonly mode (letter 'r' in the access mode) ");
		ok = nxFALSE;
	}
	else
	{
		ok = nxBinaryFile::Open( filename, "rb");			// Open the file, only open binary files
	}														// and that is that
	EmptyCacheContents();
	return ok;
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::EmptyCacheContents
 *-------------------------------------------------------------------------*/

void nxFileStrucFragments::EmptyCacheContents()
{
	m_cachepointer = m_cacheendpointer;						// resync the cachepoint so it is empty
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::eof
 *-------------------------------------------------------------------------*/

nxBOOL nxFileStrucFragments::eof()
{
	return (nxBinaryFile::eof() && (m_cachepointer >= m_cacheendpointer));
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::Seek
 *-------------------------------------------------------------------------*/

int nxFileStrucFragments::Seek( size_t offset, int origin )
{
	size_t	newpos = 0;
	size_t	nrec;
	size_t	ofs;
	nxBOOL	ok = nxTRUE;
	int		iok;
	int		status;

	switch (origin)
	{
		case SEEK_CUR : newpos = Tell() +offset;       break;				// Seek relative to current location
		case SEEK_END : newpos = FileSize() + offset;  break;				// Seek relative to end of disk
		case SEEK_SET : newpos = offset;			   break;				// Seek relative to beginning of file
		default       : ok     = nxFALSE;			   break;				// Otherwise seek current location
	};

	if (ok)
	{
		nrec = newpos/m_cacherecordsize;										// Get number of records into file
		ofs	 = newpos%m_cacherecordsize;										// get offset into file
		status = nxBinaryFile::Seek(nrec*m_totalrecordsize, SEEK_SET ); 
		ok   = ( status == 0);													// Now reposition to start of last "partial record"
		if ((ofs == 0) || !ok)													// if we dont need any of the partial record
		{																		// Then
			EmptyCacheContents();												// Empty the cache contents
		}																		// and that is that
		else																	// otherwise we need a partial record
		{																		// so
			ok = UpdateCache();													// update the cache
			if (ok) UpdateCachePointer( m_cache + ofs );						// and reposition te pointer
		}
	}
	if (ok) iok = 1;
	else    iok = 0;
	return  iok;
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::Tell
 *-------------------------------------------------------------------------*/

size_t nxFileStrucFragments::Tell()
{
	intptr_t filepos = nxBinaryFile::Tell();									// Get the current file position
	intptr_t nrecs   = filepos/m_totalrecordsize;								// Get the number of records read in

	filepos  = (nrecs-1)*m_cacherecordsize + (m_cachepointer - m_cache);	// Get the fileposition (account for stuff in buffer)
	return (size_t)filepos;
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::IsOpen
 *-------------------------------------------------------------------------*/

nxBOOL nxFileStrucFragments::IsOpen() const
{
	return (nxBinaryFile::IsOpen() && (m_cache != NULL));
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::Close
 *-------------------------------------------------------------------------*/

nxBOOL nxFileStrucFragments::Close( )
{
	nxBOOL ok;

	ok = nxBinaryFile::Close();
	EmptyCacheContents();
	return ok;
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::Read
 *-------------------------------------------------------------------------*/

size_t nxFileStrucFragments::Read( void *userbuffer, size_t tsize, size_t count)
{
	size_t	numleft = tsize*count;
	size_t	nb;
	size_t	numincache;
	nxBYTE*	userptr = (nxBYTE *)userbuffer;
	nxBYTE* ptr;
	nxBYTE*	endptr;

	while ( (numleft > 0 ) && !eof() )					// While we need more data and its not the end of the file
	{													// then
		numincache = CacheNumBytesLeft();				// Get the number of bytes left in the cache
		if (numincache == 0)							// If we have no data in the cache
		{												// then 
			UpdateCache();								// update our cache
			numincache = CacheNumBytesLeft();			// Get the number of bytes left in the cache
		}												// and that is that
		nb         = nxmin( numleft, numincache);		// Get the number of bytes required
		ptr        = CacheCurrentPointer();				// Point to current location in cache
		endptr = ptr + nb;								// Point to the end of this location
		for (; ptr < endptr; ) *userptr++ = *ptr++;		// copy the bytes over one by one.
		UpdateCachePointer(ptr);
		numleft -= nb;									// decrement the number of bytes left to read.
	}
	return (tsize*count - numleft)/tsize;				
}

/*---------------------------------------------------------------------------
 *						nxFileStrucFragments::Write
 *-------------------------------------------------------------------------*/

size_t nxFileStrucFragments::Write( const void * /*userbuffer*/, size_t /*size*/, size_t /*count*/)
{
	nxLog::Record(NXLOG_WARNING, "nxFileStrucFragments::Write, Function not supported");
	return 0;
}
