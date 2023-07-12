
#define NXMEMORYFIFO_MAX_READER	32

/*--------------------------------------------------------------------------
 *						class nxMemoryFifo									*/
/**	\ingroup system_fileio
 *	A class for streaming data using shared memory.  It is very similar to pipes
 *	except it allows the reader to keep a certain amount of history.  This is
 *	vital when decoding telemetry and we have to backtrack after losing sync.
 *	This code is only available for Windows NT.
**/
/*--------------------------------------------------------------------------*/

class nxMemoryFifo : public nxBinaryFile
{
	private:
		struct	memorylayout										// The structure of the memory
		{															// shared between processes
			nxDWORD m_RecordAtStartOfHistory;						// record number at beginning of the memory
			nxDWORD m_RecordAtEndOfHistory;							// record number at end of the memory
			nxDWORD m_lowestrecord;									// The lowest reacord (being read)
			nxDWORD	m_issynchronous;
			nxDWORD	m_EofFlag;
			nxDWORD	m_numreaders;
			nxDWORD m_size;											// the number of bytes in the cached buffer
			nxDWORD	m_synced_lookaheadsize;							// the number of bytes when looking ahead (default 1/4 of m_size)
			nxDWORD	m_readerrecord[NXMEMORYFIFO_MAX_READER];
			nxBYTE	m_data[1];										// the actual data array.
		};

	private:
		memorylayout*	m_mem;										// pointer to the shared memory
		HANDLE			m_mapobject;								// WIN32 handle to the shared memory
		nxMutex*		m_iosync;									// Mutex for exclusive access to shared memory
		nxEvent*		m_evntnewdata;								// Event indicating more data are available.
		nxEvent*		m_evntnewspace;								// Event indicating more data are available.
		nxEvent*		m_evntexists;								// Event indicating shared memory exists
		nxBOOL			m_readonly;
		nxBYTE*			m_bufferstart;
		nxBYTE*			m_bufferend;
		int				m_readerid;
		nxDWORD			m_recordid;
		DWORD			m_timeout;
		nxString		m_name;

	private:
		nxBOOL			AllocateSyncObjects( const char * name );
		nxBOOL			WaitForNewData();
		nxBOOL			WaitForSpaceInBuffer();
		nxBOOL			Create( const char *name, nxDWORD NumElements);
		nxBOOL			Attach( const char *name );
		void			UpdateLowestreadRecord();

	public:
						nxMemoryFifo();
					   ~nxMemoryFifo();
		void			SetReadTimeout( DWORD millisecs ) {m_timeout = millisecs;}
		nxBOOL			FifoIsEmpty();
		void			SynchronizeFifo( nxBOOL syncisenabled );
		void			MarkEndOfFile();
		nxBOOL			ChangeSyncedSize( nxDWORD lookahead);


	public:
		virtual	nxDWORD	FileSize() { if (m_mem != NULL) return (m_mem->m_RecordAtEndOfHistory); else return 0;}
		virtual nxBOOL	Open( const char * filename, const char * mode );
		virtual nxBOOL	IsOpen() const;
		virtual nxBOOL	Close( );
		virtual nxBOOL	eof();
		virtual size_t	Read  ( void *userbuffer, size_t tsize, size_t count);
		virtual size_t	Write ( const void *userbuffer, size_t size, size_t count);
		virtual int		Seek  ( size_t offset, int origin );
		virtual size_t	Tell  ();
};

