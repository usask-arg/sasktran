


/*---------------------------------------------------------------------------
 *					nxArrayLinear<T>::UnlockMemoryIfLocked		2003-9-17 */
/**	Unlock any memory in the manager if we have any locks in place. A lock
 *	is identified by the prescence of a valid manager and a valid memory
 *	pointer.
 */
/*-------------------------------------------------------------------------*/

template <class T>
void nxArrayLinear<T>::UnlockMemoryIfLocked( bool force)
{
	if ((m_manager != NULL) && (m_baseptr != NULL) && (force || !m_allowmemoryreuse))
	{
		m_manager->UnlockMemory();
		m_baseptr = NULL;
		m_endptr  = NULL;
	}
}

/*---------------------------------------------------------------------------
 *					nxArrayLinear<T>::ReleaseMemoryManager           2003-6-4*/
/**	Release the memory manager (if we have one) and unlock the memory if it
 *	is locked
 */
 /*-------------------------------------------------------------------------*/

template <class T>
void nxArrayLinear<T>::ReleaseMemoryManager	()
{
	UnlockMemoryIfLocked(true);							// Unlock any memory if we have any locked
	if (m_manager != NULL) m_manager->Release();	// and release the manager if we have one
}

/*---------------------------------------------------------------------------
 *					nxArrayLinear<T>::Detach						2003-9-15*/
/**	Erase the current Detach from the memory, which is basically the same as deleting all
 *	of the resources associated with this array.
 */
/*-------------------------------------------------------------------------*/

template <class T>
void nxArrayLinear<T>::Detach()
{
	bool	allow = m_allowmemoryreuse;				// Cache the allowedtogrow value

	m_allowmemoryreuse = false;						// Now set it so we are not allowed to grow
	erase();										// erase the memory, this will happen
	ReleaseMemoryManager();
	m_allowmemoryreuse = allow;						// restore the memory reuse flag
}

/*---------------------------------------------------------------------------
 *					nxArrayLinear<T>::erase							2002-9-10*/
/**	Erase the memory and layout associated with this array but keep the
 *	current memory manager alive.
 */
/*-------------------------------------------------------------------------*/

template <class T>
void nxArrayLinear<T>::erase()
{
	UnlockMemoryIfLocked(false);
	m_rankspecs.Erase();
	m_arrayisattached = false;
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::SetTo							2002-9-10*/
/**	Set each element of the array equal to the value constant.
 */
/*-------------------------------------------------------------------------*/

template <class T>
void nxArrayLinear<T>::SetTo( const T& constant )
{
	nxArrayIter<T>	iter = begin();
	nxArrayIter<T>	aend = end();

	while( iter != aend ) {*iter = constant; ++iter;};
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::ShallowCopy         2002-9-9*/
 /**	Shallow copy simply transfers the array from another array to this array
 *	and the array is shared.  It provides a fast method of moving arrays
 *	around. The alternative is #DeepCopy which makes a physical copy of the entire
 *	array which is going to be slower as arrays have to be allocated
 *	and physically copied.
 *
 *	ShallowCopy only works if the destination array is compatible with the
 *	memory organization of the source array.  For example it is not
 *	to translate a 3-D array with variable strides to a 2-D array. ShallowCopy
 *	is used extensively in the copy constructor as it allows us to
 *	pass an array by value without too much overhead.
 *
 *	ShallowCopy performs the following steps,
 *		- ShallowCopy will return successfully, with no action taken if it is copying to itself
 *		- ShallowCopy will erase itself if the other array is empty.
 *		- ShallowCopy will fail if it cannnot reshape this array dimensions to be compatible with the other arrays dimensions and strides
 *		- ShallowCopy will increment the memory lock of the other arrays memory manager and add a
 *		  reference onto the memory manager itself.
 *		- ShallowCopy will attach this array if the other array is also attached
 *
 *	Note that the other array cannot change the size of the shared memory until this array
 *  unlocks its hold on that memory.
 */
/*-------------------------------------------------------------------------*/

template <class T>
bool nxArrayLinear<T>::ShallowCopy(const nxArrayLinear<T>& other)
{
	bool					ok;
	bool					oldislocked = (m_baseptr != NULL);
	InxMemoryManager<T>*	oldmanager  = m_manager;			// Save the pointer to our memory manager

	if (this == &other)											// Shallow copying this object to itself
	{															// does nothing
		ok = true;
	}
	else
	{
		ok = !m_arrayisattached;								// See if this array is attached to memory somewhere else
		if (!ok)												// IF it is then we have no business
		{														// trying to assign this pointer 
			nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::ShallowCopy, Cannot shallowCopy to an array that is attached to memory. Use Detach or erase before invoking ShallowCopy");
		}
		else
		{
			if (other.IsEmpty())
			{
				erase();
				ok = true;
			}
			else
			{
				ok = m_rankspecs.ReshapeToMandatoryRank( true, MandatoryRankValue(), sizeof(T), *other.ArrayRankSpecs());
				if (!ok)													// as we cannot copy 1-D arrays into 3-D arrays etc.
				{															// if there is a problem log the message
					nxLog::Record( NXLOG_WARNING, "nxArrayLinear<T>::ShallowCopy, Cannot perform copy as there is a problem reshaping the memory layout to the desired dimensionality");
					erase();												// erase the current array, but leave the memory manager alone.
				}															// and that is that
				else														// otherwise we are okay
				{															// so
					m_manager          = other.m_manager;					// Copy over the manager
					m_arrayisattached  = other.m_arrayisattached;			// If the other array is attached to memory then so is this shallow copy.
					if (m_manager != NULL)									// and if we have a real manager
					{														// then
						m_manager->AddRef();								// Add a reference to the memory manager
						m_manager->LockMemory();							// And increase the lock count on the managed memory
					}
					ok = ConfigureMemoryLayout(other.m_baseptr);			// and configure the
					if (!ok)
					{
						nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::ShallowCopy, Error configuring memory layout.  Something is dysfunctional");
					}
					if (oldmanager != NULL)
					{
						if (oldislocked) oldmanager->UnlockMemory();
						oldmanager->Release();
					}
				}
			}
		}
	}
	if (!ok) Detach();
	return ok;
}


/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::DeepCopy						2002-9-9*/
 /** The purpose is to physically copy the source array to a new array.
 *   Physically copying executes the following rules:
 *
 *	 -  DeepCopy will not copy an array to itself. It will test and immediately
 *		return.
 *
 *	 -  If the source array has the temporary flag set then it is assumed this
 *		array is part of a sequence of binary operators (as only they should create
 *		nxArrayTemporary) and DeepCopy invokes ShallowCopy.
 *
 *	 -	DeepCopy will not change the memory layout (rank, stride and location) of this array
 *		if it has exactly the same number of elements as the source.  This is one method that
 *		allows you to copy the contents from any dimensional array to any other dimensional array
 *
 *	 -	DeepCopy will not change the destination memory layout if it is attached to memory.
 *	 -	DeepCopy is free to change the entire memory layout (rank, stride and location) of the destination
 *		array.  This includes emptying the destination if the source is empty.
 *   .
 */
/*-------------------------------------------------------------------------*/

template <class T>
bool nxArrayLinear<T>::DeepCopy(const nxArrayLinear<T>& other, bool copycontents)
{
	bool				ok;

	if (&other == this) return true;										// Dont copy an array onto itself
	if (other.IsTemporary() && (m_baseptr == NULL) )						// IF we have no memory allocated
	{																		// and the other guy is temporary
		ok = ShallowCopy(other);											// then do a quick shallow copy (i.e. transfer)
	}																		// and that is that
	else																	// otherwise
	{
		if (m_arrayisattached)												// If this array is attached to memory							
		{																	// then it must have identical
			ok = (size() == other.size());									// number of elements
			if (!ok)														// as anything else requires a change of size
			{
				nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::DeepCopy, Cannot DeepCopy to an attached array that does not have the same number of elements. This (attached) array has (%u) elements while the other array has (%u) elements.", (unsigned int)size(), (unsigned int)other.size() );
			}
		}
		else
		{
			ok = SetSize( (size_t)other.ArrayRankSpecs()->Rank(), other.ArrayRankSpecs()->Dims() );		// Always call set size to reconfigure the rank.
			if (!ok)
			{
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::DeepCopy, failed to due to memory allocation error");
				}
			}
		}
		if (ok && copycontents)												// Generally user will want contents copied but some of the operator macros dont as we are going to overwrite immediately
		{
			nxArrayIter<T>	otherb = other.begin();
			nxArrayIter<T>  iter   = begin();
			nxArrayIter<T>	iend   = end();
			while (iter != iend) {*iter = *otherb; ++iter; ++otherb;}
		}
	}
	if (!ok) Detach();
	return ok;
}


/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>:CheckRankIsCompatible        2002-9-11 */
 /** Check that the anticpated rank of this array is compatible with the
 *	 mandatory rank of this array (as defined by the derived base class
 */
/*-------------------------------------------------------------------------*/

/*
template <class T>
bool nxArrayLinear<T>::CheckRankIsCompatible()
{
	bool ok;
	int		mandrank = MandatoryRankValue();

	ok = (mandrank < 1) || ( mandrank == m_rankspecs.Rank());
	if (!ok)																	// Generally these ranks are not compatible
	{																			// however there is one special case that is worthy
		if (mandrank == 1)														// if this array must be 1-D
		{																		// then we can covert a multi-dimensional into this
			ok = m_rankspecs.IsFixedStride();									// as long as its strides are fixed
			if (ok)																// if they are fixed
			{																	// then
				size_t			n      = m_rankspecs.N_Elements();				// get the number of elements
				size_t	stride = m_rankspecs.Stride();						// get the stride

				m_rankspecs.Configure( 1, &n, sizeof(T), &stride );				// and reconfigure the rankspec object
			}																	// and that is that
			else																// but if the source is not fixed stride
			{																	// then it cannot be compatible
				nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::CheckRankIsCompatible, The source array cannot be converted to a 1-D array because the source array's strides are not fixed");
			}																	// and that is that
		}																		// and that is the end of the special case
		else																	// if its not the special case
		{																		// then it just cannot be done. so log a message
			nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::CheckRankIsCompatible, An array of rank %d cannot be converted to an array of rank %d", (int)m_rankspecs.Rank(), (int)mandrank);
		}
	}
	return ok;
}
*/
/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::LastAddressableLocation		2003-9-14 */
 /**	Gets the address of the last location of this array. This is unlike
 *	m_endptr which is the first address after the last location
 */
/*-------------------------------------------------------------------------*/

/* ************ NEEDS WORK. m_baseptr and m_endptr need proper offset control. */

template< class T>
T* nxArrayLinear<T>::LastAddressableLocation() const
{
	const size_t*	dims		= m_rankspecs.Dims();
	const size_t*	strides		= m_rankspecs.Strides();
	size_t			nrank		= m_rankspecs.Rank();
	size_t			offset;
	T*				ptr;

	if (nrank < 1)
	{
		ptr = NULL;
	}
	else
	{
		offset    = 0;
		for (size_t i = 0; i < nrank; i++) offset += (dims[i]-1)*strides[i];
		ptr  = (T*)((nxBYTE *)m_baseptr + offset);
	}
	return ptr;
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::ConfigureMemoryLayout		2003-9-14
 *-------------------------------------------------------------------------*/

/* ************ NEEDS WORK. m_baseptr and m_endptr need proper offset control. */

template< class T>
bool	nxArrayLinear<T>::ConfigureMemoryLayout( T* ptr )
{
	const size_t*		dims		= m_rankspecs.Dims();
	const size_t*	strides		    = m_rankspecs.Strides();
	size_t				rr          = m_rankspecs.Rank();
	size_t			    r		    = m_rankspecs.Rank()-1;
	bool				ok			= (rr >= 1);

	if (ok)
	{
		m_baseptr = ptr;
		m_endptr  = (T*)((nxBYTE *)m_baseptr + dims[r]*strides[r]);
		ConfigureIndexingOperator();
	}
	else
	{
		m_baseptr = ptr;
		m_endptr  = ptr;
		m_indextopointer = &nxArrayLinear<T>::IndexToPointer_EmptyArray;
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::SetMemoryManager               2003-6-7
 *-------------------------------------------------------------------------*/

template <class T>
bool nxArrayLinear<T>::SetMemoryManager( InxMemoryManager<T>* manager, T* baseptr )
{
	if (manager != NULL) manager->AddRef();			// Add a reference to the new memory manager
	ReleaseMemoryManager();							// release the old memory manager.
	m_manager = manager;
	m_baseptr = baseptr;
	if (m_manager != NULL) m_manager->LockMemory();
	return true;
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::InternalAttach		2003-9-15
 *-------------------------------------------------------------------------*/

template< class T>
bool	nxArrayLinear<T>::InternalAttach( size_t nrank, const size_t* dims, T* ptr, const size_t* strides, InxMemoryManager<T>* manager  )
{
	bool	ok1;
	bool	ok2;
	bool	ok3;
	bool	ok;
	RankSpecification sourcerank;

	sourcerank.Configure( nrank, dims, sizeof(T), strides );

	ok1 = m_rankspecs.ReshapeToMandatoryRank( true, MandatoryRankValue(), sizeof(T), sourcerank );	// Configure the rank specifications,
//	ok1 = m_rankspecs.Configure( nrank, dims, sizeof(T), strides );		// Configure the dimensions and strides organization
	ok2 = SetMemoryManager     ( manager,  ptr );						// Configure the memory manager (i.e. get rid of existing)
	ok3 = ConfigureMemoryLayout( ptr );									// Configure
	ok  = ok1 && ok2 && ok3;
	if (ok && (m_manager != NULL))
	{
		ok = manager->CheckMemoryRange( m_baseptr, LastAddressableLocation() );
		if (!ok)
		{
			nxLog::Record(NXLOG_ERROR, "nxArrayLinear<T>::InternalAttach, the requested array is not within the bounds of the memory manager");
			NXASSERT(ok);
		}
	}
	m_arrayisattached = ok;
	if (!ok) Detach();
	return (ok);
}

/*---------------------------------------------------------------------------
 *					nxArrayLinear<T>::Attach					2002-9-12*/
/**	Attach to a chunk of memory. If successful then this object will use the
  * rank specifications and memory specified in the parameters.
  *
  *	\param nrank
  *	The rank of the new array.
  *
  *	\param dims
  *	 An array of at least nrank elements. Each element specifies the size of each dimension.
  *	 Must not be NULL.
  *
  *	\param ptr
  *	Points to the first element of the array.
  *
  *	\param manager
  *	The memory manager used to control the memory used by the array.  If it is not NULL then the
  *	lockcount of the memory manager is incremented.  This parameter can be NULL in which case the
  *	array accesses memory without a memory manager and the user is responsible for ensuring the accesses memory is valid in memory.  This parameter has a default
  *	value of NULL.  If the memory manager is not NULL then the instance will continue using this memory
  *	manager until #Detach() is called or the object is destroyed.
  *
  *	\param strides
  *	An array of at least nrank elements (or NULL if not used). Species the stride in bytes of each
  *	dimension of the array.  If this parameter is NULL then the array strides are configured to be contiguous.
 */
/*-------------------------------------------------------------------------*/

template< class T>
bool nxArrayLinear<T>::Attach( size_t nrank, const size_t* dims, T* ptr, InxMemoryManager<T>* manager,  const size_t* strides )
{
	return InternalAttach( nrank, dims, ptr, strides, manager  );
}

/*---------------------------------------------------------------------------
 *					nxArrayLinear<T>::SetSize					2002-9-12*/
/**	This function sets the size of the array.  It is treated as a request
 *	to discard the current memory and request new memory using the current
 *	memory manager. The new memory will be contiguous. Note that calls to
 *	SetSize on arrays that are attached will fail unless the two arrays
 *	have compatible memory layouts.
 *
 */
/*-------------------------------------------------------------------------*/


template< class T>
bool nxArrayLinear<T>::SetSize( size_t nrank, const size_t* dims, const size_t* strides)
{
	bool				ok;
	bool				willbeempty;
	T*					ptr;
	int					mandatoryrank;

	ok = m_rankspecs.IsSameLayout( nrank, dims, strides );		// Are these two layouts the same
	if (!ok)													// Nope
	{															// So something must change
		ok = !m_arrayisattached;								// that means we must not be Attached to memory
		if (!ok)												// If we are
		{														// then issue an error
			nxLog::Record(NXLOG_WARNING,"nxArrayLinear<T>::SetSize, you cannot change the size of an array that is Attached, you must call Detach or erase first");
			NXASSERT(false);
		}
		else										
		{
			if (( m_manager != NULL) && m_manager->MemoryIsShared()) Detach();		// If we are sharing the memory with others then we must detach, from the memory manager as we cannot change its size
			
			mandatoryrank = MandatoryRankValue();									// Get the mandatory rank of this array 
			ok = m_rankspecs.Configure( nrank, dims, sizeof(T), strides);			// Configure the dimensions and strides passed in 
			if ((nrank != (size_t)mandatoryrank) && (mandatoryrank != 0))					// If the user has requested the wrong size rank
			{																									// If the user has requested a rank that is not this dimension
				ok = ok && m_rankspecs.ReshapeToMandatoryRank( false, mandatoryrank, sizeof(T), m_rankspecs);	// Reshape newrank to the mandatory rank
			}
			if (!ok)																// If it failed then let the user know
			{																		// because this should not fail unless things are quite dysfunctional
				nxLog::Record( NXLOG_WARNING, "nxArrayLinear<T>::SetSize, Error configuring the rankspecification object for rank (%d)", (int)nrank);
			}
			else															// we have a default rank specification
			{																// now we can go about re-allocating the memory
				willbeempty = (m_rankspecs.N_Elements() < 1);				// See if the new array will be empty
				if (willbeempty)											// It will be empty
				{															// so
					erase();												// erase the array
				}
				else														// otherwise we will be allocating some memory
				{															// so
					if (m_manager == NULL)									// If we have no memory manager
					{														// Attached to this object
						m_manager = new InxMemoryManager<T>;				// then create a default implementation (uses C++ new)
						ok = (m_manager != NULL);							// and make sure it worked
						if (!ok) nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::SetSize, Error allocating InxMemoryManager<T>");
						else     m_manager->AddRef();						// lock the manager
						m_baseptr = NULL;									// and make sure our memory pointer is not defined to something else
						m_endptr = NULL;
					}
					if (ok)													// so we now have a valid memory manger and we have a requirement to allocate some memory
					{														// so
						m_baseptr = NULL;																// make sure our memory pointer is not defined
						m_endptr  = NULL;																// and clear the end pointer while we are at it.
						ok = m_manager->AllocateAndLock( &m_rankspecs, m_allowmemoryreuse, &ptr);		// Allocate the number of elements
						if (!ok)
						{
							nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::SetSize, Error allocating memory");
						}
						else
						{
							ok  = ConfigureMemoryLayout( ptr );
							ok  = ok && m_manager->CheckMemoryRange( m_baseptr, LastAddressableLocation() );
							if (!ok)
							{
								nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::SetSize, Error checking rank and dimensions of array");
							}
						}
					}
				}
			}
		}
	}
	if (!ok) Detach();
	return ok;
}



/*-----------------------------------------------------------------------------
 *					nxArrayLinear<T>::CheckSlicingDims		2005-1-7*/
/** **/
/*---------------------------------------------------------------------------*/

template< class T>
bool nxArrayLinear<T>::CheckSlicingDims( const size_t indexlo[], const size_t indexhi[], size_t rank, size_t* numsubsetdims ) const
{
	bool				ok;
	bool				dimok;
	const size_t*		dims;
	size_t				lowindex;
	size_t				hihindex;
	size_t				i;
	size_t				nelem;

	*numsubsetdims = 0;														// reset the number of subset dimensions
	ok = (rank == m_rankspecs.Rank());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::CheckSlicingDims, Rank mismatch. This array has rank (%d) while the passed indices have rank (%d)", (int)rank, (int)m_rankspecs.Rank() );
	}
	else
	{
		dims          = m_rankspecs.Dims();										// get the dimensions array of the slice array
		for (i=0; i < rank; i++)											// First scan through the array and check that everything is kosher.
		{																		// check that the indexing is sensible
			lowindex = indexlo[i];												// get the lower index
			hihindex = indexhi[i];												// get the upper index
			if (lowindex == NXARRAY_STARSELECT )lowindex = 0;										// If its negative then extend to the lower bound
			if (hihindex == NXARRAY_STARSELECT )hihindex = (dims[i]-1);						// If its negative then extend to the upper bound
			dimok = (hihindex < dims[i]) && (lowindex <=hihindex);
			if (dimok)
			{
				nelem = (hihindex - lowindex) + 1;									// Get the number of elements sliced out of this dimension
				if (nelem > 1) ++(*numsubsetdims);				// If we have more than 1 elements (or forced dimension from -ve number)
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::CheckSlicingDims, Slicing indices (%d .. %d) of dimension(%d) not valid for dimension range(0..%d)", (int)lowindex, (int)hihindex, (int)(dims[i]-1));
				ok = false;
			}
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::CheckSlicingDims, There was an error with the slicing specifications.");
			*numsubsetdims = 0;
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxArrayLinear<T>::Slice		2005-1-7*/
/** Slices a sub-section of an array.
 *	\par Slice indexing
 *	Slice indexing is based upon the 0 based indexing normally used in this
 *	array class but there are two special instances. Normally a slice where
 *	the start and finish index are the same (i.e. a length of 1 element)
 *	indicates the "colum" or "row" to be extracted and the code will normally
 *	collapse this dimension.
 *
 *	\param indexlo
 *		An array of integers, indexlo[rank]. Each entry is the 0 based
 **/
/*---------------------------------------------------------------------------*/
/* ************ NEEDS WORK. GOT TO CHECK THE WHOLE DARN THING . */

template< class T>
bool nxArrayLinear<T>::Slice( const size_t indexlo[], const size_t indexhi[], size_t rank, nxArrayLinear<T>* column) const
{
	bool				ok;
	T*					baseptr;
	size_t			baseoffset;
	size_t*				slicedimsptr;
	size_t*				slicestrideptr;
	const size_t*		dims;
	const size_t*	strides;
	size_t					numsubsetdims;
	size_t					i;
	size_t					subidx;
	size_t				lowindex;
	size_t				hihindex;
	size_t				nelem;
	bool				isscalar;

	ok = CheckSlicingDims( indexlo, indexhi, rank, &numsubsetdims );
	if (ok)
	{
		isscalar =  (numsubsetdims == 0);												// Rank 0 if all dimension extract just one element
		if (isscalar) numsubsetdims++;													// If that is the case then its goeing to be a rank one array

		ok       = column->m_rankspecs.AllocateUninitializedSpace( numsubsetdims );		// Allocate the space for the rank of this slice
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::Slice, Error allocating rank storage area for sliced array");
		}
		else
		{
			slicedimsptr    = column->m_rankspecs.UnsafeDims();					// get the dimensions array of the slice array
			slicestrideptr  = column->m_rankspecs.UnsafeStrides();				// get the strides array of the slice array
			dims            = m_rankspecs.Dims();								// get the dimensions array of the slice array
			strides         = m_rankspecs.Strides();							// get the strides array of the slice array
			baseoffset = 0;
			subidx     = 0;
			for (i = 0; (i < rank); i++)										// Scan through all of the indices in this array
			{																	// and for each dimension
				lowindex = indexlo[i];											// get the lower index
				hihindex = indexhi[i];											// get the upper index
				if (lowindex == NXARRAY_STARSELECT )lowindex = 0;									// Handle the special STAR value
				if (hihindex == NXARRAY_STARSELECT )hihindex = (dims[i]-1);				// for both lower and upper bounds
				nelem = (hihindex - lowindex) + 1;								// Get the number of elements in this dimension
				if (nelem > 1)													// and see if we do need to copy over this dimension
				{																// if we do
					slicedimsptr   [subidx]   = nelem;							// Then copy over the number of elements
					slicestrideptr [subidx]   = strides[i];						// and the strides
					subidx++;													// point to the then next valid output dimensions
				}																// otherwise ignore this dimension
				baseoffset += lowindex*strides[i];								// get the offset to the first element of this dimension
			}
			baseptr = (T*)( (nxBYTE*)m_baseptr + baseoffset);					// Point to the memory

			if (isscalar)														// If we have a scalar stride
			{																	// then patch it up
				slicedimsptr   [0]   = 1;										// so its a
				slicestrideptr [0]   = strides[0];
			}

			ok = column->InternalAttach( numsubsetdims, slicedimsptr, baseptr, slicestrideptr, m_manager  );
			if (!ok) column->Detach();
		}
	}
	if (!ok)
	{
		column->erase();
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::SetTemporaryArrayParams    2003-6-4
 *-------------------------------------------------------------------------*/

template< class T>
bool nxArrayLinear<T>::SetTemporaryArrayParams( InxMemoryManager<T>* manager, const RankSpecification* rank, T* baseptr )
{
	bool ok = (m_istemporary);
	if (ok)
	{
		NXASSERT( (manager != NULL) && (baseptr != NULL) );
		m_rankspecs = *rank;
		m_manager   = manager;
		m_manager->AddRef();
		m_manager->LockMemory();
		ConfigureMemoryLayout( baseptr );
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					nxArrayLinear<T>::TrimSize		2004-12-13*/
/** Trims the size of the last dimension in the array**/
/*---------------------------------------------------------------------------*/

template< class T>
bool nxArrayLinear<T>::TrimSize( size_t nx)
{
	bool			ok;
	RankSpecification	copy(m_rankspecs);
	size_t*			dims;
	size_t			ndims;

	ndims = copy.Rank();
	dims  = copy.UnsafeDims();
	ok = (ndims > 0 );
	if (ok)
	{
		ok = (nx <= dims[ndims-1] );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "nxLinearArray<T>::TrimSize, you cannot trim an array to a BIGGER size without destroying old array");
		}
		else
		{
			dims[ndims-1] = nx;
			SetReuseMemory( true );
			SetSize(ndims, dims);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxArrayLinear<T>::At		2004-12-13*/
/** **/
/*---------------------------------------------------------------------------*/

template< class T>
bool nxArrayLinear<T>::At( nxArrayLinear<size_t>& logicalindex, nxArrayLinear<T>* answer)
{
	size_t	diff;
	size_t		lastindex;
	size_t		thisindex;
	size_t			n = logicalindex.size();
	nxArrayIter<T>	iter;
	bool			ok = false;

	if (n < 1)
	{
		ok = true;
		answer->erase();
	}
	else
	{
		ok = answer->SetSize(1,&n, NULL);
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "nxArrayLinear<T>::At(logicalindex), failed when setting size of index");
		}
		else
		{
			diff = 0;
			lastindex = 0;
			iter      = logicalindex.begin();
			for (size_t i=0; i < n; i++)
			{
				thisindex = *iter;
				diff      =  thisindex - lastindex;
				iter     += diff;
				lastindex = thisindex;
				++iter;
			}
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					nx1dArray<T>::STLVector							2007-8-8
 *-------------------------------------------------------------------------*/

/* ************ NEEDS WORK. ADAM JUST STUCK THIS IN WITHOUT TALKING TO NICK. */

template< class T>
std::vector<T> nx1dArray<T>::STLVector( ) const
{
	std::vector<T> v;
	v.reserve( this->size() );

	nxArrayIter<T>	iter   = this->begin();
	nxArrayIter<T>	finish = this->end();

	while (!(iter == finish))
	{
		v.push_back( *iter );
		++iter;
	}

	return( v );
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::CheckBounds		2003-9-13
 *	Check the bounds of the indexes provided by the user.
 *-------------------------------------------------------------------------*/

template <class T>
bool nxArrayLinear<T>::CheckBounds( const size_t* dimsindex) const
{
	const size_t*	rankdims = m_rankspecs.Dims();
	bool		ok       = m_rankspecs.Rank() > 0;

	for (size_t i= 0; i < m_rankspecs.Rank(); i++)
	{
		ok = ok && ((dimsindex[i] >= 0) && (dimsindex[i] < rankdims[i]));
	}
	if (!ok)									// Ok so the users indices are out of bounds
	{
		nxString	s;
		nxString	users;

		s.sprintf    ("[%1d", (int)rankdims[0]);
		users.sprintf("[%1d", (int)dimsindex[0]);

		for (size_t i= 1; (i < m_rankspecs.Rank()); i++)
		{
			s.sprintf    ("%s,%1d", (const char*)s,     (int)rankdims[i]);
			users.sprintf("%s,%1d", (const char*)users, (int)dimsindex[i]);
		}
		s.sprintf    ("%s]", (const char*)s);
		users.sprintf("%s]", (const char*)users);
		nxLog::Record(NXLOG_ERROR, "nxArrayLinear<T>::CheckBounds, Users indices %s are out of bounds of %s", (const char*)users, (const char *)s );
		NXTRACE(("nxArrayLinear<T>::CheckBounds, Users indices %s are out of bounds of %s\n", (const char*)users, (const char *)s));
	}
	return ok;
}



/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::nxArrayLinear		2003-9-13
 *-------------------------------------------------------------------------*/

template <class T>
nxArrayLinear<T>::nxArrayLinear()
{
	ConstructorInit();
}


/*-----------------------------------------------------------------------------
 *					nxArrayLinear<T>::ConstructorInit		2005-1-7*/
/** **/
/*---------------------------------------------------------------------------*/

template <class T>
void nxArrayLinear<T>::ConstructorInit()
{
	m_manager         = NULL;
	m_istemporary     = false;
	m_baseptr         = NULL;
	m_endptr          = NULL;
	m_allowmemoryreuse = false;
	m_arrayisattached  = false;
	m_indextopointer  = &nxArrayLinear<T>::IndexToPointer_EmptyArray;
#if defined(NXDEBUG)
	m_checkbounds     = true;
#else
	m_checkbounds     = false;
#endif
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::nxArrayLinear		2003-9-17
 *	Copy Constructor
 *-------------------------------------------------------------------------*/

template <class T>
nxArrayLinear<T>::nxArrayLinear(const nxArrayLinear<T>& other)
{
	ConstructorInit();
	m_istemporary     = other.m_istemporary;
	m_checkbounds     = other.m_checkbounds;
	DeepCopy(other);
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::IndexToLinearGeneral		2003-9-13
 *	Calculates a pointer to the element specified by the index in dimsindex.
 *	All indexing is assumed to be zero based.
 *
 *	Note that this option does not do any bounds checking.
 *-------------------------------------------------------------------------*/

/* ************ NEEDS WORK. GOT TO CHECK OFFSETS . */

template <class T>
T* nxArrayLinear<T>::IndexToPointer_General( const size_t* dimsindex) const
{
	size_t			offset  = 0;
	const size_t*	strides = m_rankspecs.Strides();

	if (m_checkbounds) CheckBounds( dimsindex);
	for (size_t i= 0; i < m_rankspecs.Rank(); i++)
	{
		offset	+= (*dimsindex++)*(*strides++);
	}
	return (T*)( (nxBYTE*)m_baseptr + offset);
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::IndexToPointer_EmptyArray		2003-9-13
 *	Indexing for empty arrays.  Returns a pointer to a piece of memory
 *	(not thread safe) but at least not an invald pointer.
 *-------------------------------------------------------------------------*/

template <class T>
T* nxArrayLinear<T>::IndexToPointer_EmptyArray(  const size_t* /*dimsindex*/ ) const
{
	static T	dummy;

	nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::IndexToPointer_EmptyArray, You are indexing an empty array");
	return &dummy;
}

template <class T>
T* nxArrayLinear<T>::IndexToPointer_1D_Contiguous(  const size_t* dimsindex ) const
{
	if (m_checkbounds) CheckBounds( dimsindex);
	return m_baseptr + *dimsindex;
}
/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::IndexToPointer_1D_Fixed		2003-9-13
 *	Fast indexing for 1-D contiguous arrays
 *-------------------------------------------------------------------------*/

template <class T>
T* nxArrayLinear<T>::IndexToPointer_1D_Fixed(  const size_t* dimsindex ) const
{
	if (m_checkbounds) CheckBounds( dimsindex);
	return (T*)( (nxBYTE*)m_baseptr + (*dimsindex)*m_rankspecs.Stride());
}

/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::IndexToPointer_2D_Contiguous		2003-9-13
 *-------------------------------------------------------------------------*/

template <class T>
T* nxArrayLinear<T>::IndexToPointer_2D_Contiguous(  const size_t* dimsindex ) const
{
	if (m_checkbounds) CheckBounds( dimsindex);
	return m_baseptr + ((m_rankspecs.Dims()[0])*dimsindex[1] + dimsindex[0]);
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::IndexToPointer_2D_Fixed		2003-9-13
 *	Fast indexing for 2-D contiguous arrays
 *-------------------------------------------------------------------------*/

template <class T>
T* nxArrayLinear<T>::IndexToPointer_2D_Fixed(  const size_t* dimsindex ) const
{
	const size_t*	stride = m_rankspecs.Strides();

	if (m_checkbounds) CheckBounds( dimsindex);
	return (T*)( (nxBYTE*)m_baseptr + (stride[0]*dimsindex[0] + stride[1]*dimsindex[1]));
}


/*---------------------------------------------------------------------------
 *'					nxArrayLinear<T>::ConfigureIndexingOperator		2003-9-14
 *-------------------------------------------------------------------------*/

template <class T>
void nxArrayLinear<T>::ConfigureIndexingOperator()
{
	if (m_rankspecs.N_Elements() < 1)
	{
		m_indextopointer = &nxArrayLinear<T>::IndexToPointer_EmptyArray;
	}
	else if (m_rankspecs.IsContiguous())
	{
		switch ( m_rankspecs.Rank())
		{
		case 1:		m_indextopointer = &nxArrayLinear<T>::IndexToPointer_1D_Contiguous;
					break;

		case 2:		m_indextopointer = &nxArrayLinear<T>::IndexToPointer_2D_Contiguous;
					break;

		default:	m_indextopointer = &nxArrayLinear<T>::IndexToPointer_General;
					break;
		};
	}
	else if (m_rankspecs.IsFixedStride())
	{
		switch ( m_rankspecs.Rank())
		{
		case 1:		m_indextopointer = &nxArrayLinear<T>::IndexToPointer_1D_Fixed;
					break;

		case 2:		m_indextopointer = &nxArrayLinear<T>::IndexToPointer_2D_Fixed;
					break;

		default:	m_indextopointer = &nxArrayLinear<T>::IndexToPointer_General;
					break;
		};
	}
	else
	{
		m_indextopointer = &nxArrayLinear<T>::IndexToPointer_General;
	}
}





/*-----------------------------------------------------------------------------
 *					InputXYStructExText		2004-11-24*/
/** Reads 2-D column major input from a stream until it reaches the end of the
 *	file or has read in the users requirements. The code assumes all columns of
 *	text are of the same type and are read in
 *	using the streaming operator >>.  The code will only succeed if this instance
 *	of nxLinearArray supports creation of a 2-D array. This probably means it is either an
 *	instance of nx2dArray or nxLinearArray. The first, fastest changing, index
 *	of the array is derived from one column of text "down" the file. The second
 *  slowly changing index is each column across the line of text. The code also assumes that
 *	each line of numerical text contains exactly the same number of columns as all other
 *	lines.
 *
 *	\par Text Layout and Comments
 *	I have written the code to so simple comments may be added to the
 *	text files. I find this helps keep track of what the text files
 *	contain. Comments are special characters that act like "end of line": anything
 *	after the comment char is ignored until the next new line. Comments are
 *	identified by one of the 3 characters ";%#". Blank lines or lines with only whitespace
 *	are ignored.  Individual numbers should be separated by whitespace. Commas
 *	may be used between numbers but they are interpreted as white space and not as
 *	text delimiters: a double comma construct such as "45, ,68" is interpreted as
 *	the two numbers 45 and 68 and not the 3 numbers 45, 0, 68. The number formats must
 *	be consistent with the streaming operator >> and the array template type (T).
 *
 *	\par Details
 *	The code scans once through the stream to count/check the number of rows and
 *	a second time to read in the data. The code will parse the file, if required, to
 *	count the number of columns and rows in the file.
 *
 *	\param instr
 *	The text input stream.  This a stream for lines of text. Each line is read
 *	in one line at a time and parsed. The stream must support the ability to
 *	reposition the stream pointers using seekg.  Upon exit the stream will point beyond
 *	the last line processed. This is usually the end-of-file if the user does not specify
 *	the number of columns and number of lines to read in.
 *
 *	\param numcolumns
 *	The number of columns across the line. If this value is set to zero, which is the default,
 *	then the code will scan through the file and count the number of columns.  It is
 *	recommended the user set this value if they know the number of columns that
 *	should be in the file.
 *
 *	\param numdatalines
 *	The number of valid data lines in the file. This is the same as the number of lines
 *	in the stream after all commented and blank lines are removed. If this value is set
 *	to zero, which is the default, then the code will scan through the file and count
 *	the number of valid data lines.  It is recommended that the user set this value
 *	if they know, a-priori, the number of rows that should be in the file.
 *
 *	\return
 *	Returns true if the routine managed to read in the requested data. Note it is
 *	possible for this routine to return "SUCCESS" even after encountering badly
 *	formatted lines
**/
/*---------------------------------------------------------------------------*/

template <class T>
bool nxArrayLinear<T>::InputColumnMajorText(  std::istream& instr, size_t numcolumns, size_t numdatalines )
{
	size_t							lineno;
	nxString						oneline;
	nxArrayIter<T>					iter;
	size_t							numx;
	size_t							numy;
	size_t							xctr;
	size_t							dims[2];
	std::istream::pos_type			startpos;
	bool							ok = true;
	size_t							iy;
	size_t							stride;

	// ---- Scan through the file to get the size
	// ---- unless the caller has already specified these values.

	numx  = numcolumns;												// number of columns of valid, uncommented, text
	if ((numdatalines ==0) || (numx == 0))							// If the user chooses not to specifyu both params
	{																// then scan through the stream
		numy  = 0;													// reset the number of lines of valid text
		startpos = instr.tellg();									// record the position of the stream
		lineno   = 0;												// track the actual line number
		while ( ok && !instr.eof() && ((numdatalines == 0) || (numy < numdatalines)) )		// Scan until end of file, or until we have all the lines the user wants
		{															// so
			lineno++;												// increment the line number
			ok = oneline.InputTextLine( instr);						// input a complete line of text
			if (!ok)												// if it did not work
			{														// properly then something is messed up
				nxLog::Record(NXLOG_WARNING, "nxLinearArray::InputColumnMajorText, error parsing a line of text from the stream. This should not happen with normal text files." );
			}														// otherwise
			else													// we have a line of text
			{														// so
				oneline.TrimToComments(";%#");						// Trim comments and check that its not all white space
				if (!oneline.IsEmpty())								// If we have anything left
				{													// then
					oneline.Replace( ',', ' ');						// eliminate any commas on the line
					xctr = oneline.CountNonWhiteFieldsOnLine();		// count the number of fields on the line
					if (numx == 0) numx = xctr;						// set the number of fields if the user has not set it already
					ok = (xctr == numx);							// make sure the columns across the line are correct
					if (!ok)										// If they are not
					{												// then complain
						nxLog::Record(NXLOG_WARNING, "nxLinearArray::InputColumnMajorText, error parsing file, inconsistent number of fields on line %d, expected %d fields but detected %d", (int)lineno, (int)numx, (int)xctr );
					}												// and that is that
					else numy++;									// otherwise this is a valid data line, so count it in our running total
				}													// finished processing this line
			}														// finsihed processing o fthis line of text
		}															// scan through stream until done
		instr.clear();												// once done clear the error bits
		instr.seekg(startpos);										// and reposition stream
	}																// and that is the end of the file scan
	else															// otherwise if the user knows the number of columns and number of data lines
	{																// then
		numx = numcolumns;											// set the number of columns
		numy = numdatalines;										// and the number of datalines
	}

	// ----- We have parsed the file for columns and rows
	// ----- Now allocate memory and load it in.

	if (ok)
	{
		ok = ( (numy > 0) && (numx > 0) );								// Check that we have at least
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "nxLinearArray::InputColumnMajorText, cannot read in array as either no columns or no valid lines wee found in this stream");
		}
		else
		{
			dims[0]  = numy;											// The fastest index in the new array is the column down the file
			dims[1]  = numx;											// The slowest index in the new array are the fields across the line of text
			ok = SetSize( 2, dims );									// Create the new array
			if (!ok)													// if it failed
			{															// then report a message
				nxLog::Record(NXLOG_WARNING, "nxLinearArray::InputColumnMajorText, This instance could not create a %d x %d array. The requested size is propbably too big or this instance does not support 2D arrays", (int)numx, (int)numy );
			}															// and that is that
			else														// otherwise we have memory allocated
			{															// so
				iy = 0;													// reset the numbe rof lines read in
				lineno = 0;
				while (ok && !instr.eof() && (iy < numy) )
				{
					lineno++;
					ok = oneline.InputTextLine( instr);
					if (!ok)												// if it did not work
					{														// properly then something is messed up
						nxLog::Record(NXLOG_WARNING, "nxLinearArray::InputColumnMajorText, error reading a line of text from the stream. This should not happen with normal text files." );
					}														// otherwise
					else													// we have a line of text
					{														// so
						oneline.TrimToComments(";%#");						// Trim Comments and check that its not all white space
						if (!oneline.IsEmpty())								// If we have any string left
						{													// then
							xctr = oneline.CountNonWhiteFieldsOnLine();		// count the number of fields on the line
							ok = (numx == xctr);							// see if we got the correct number of fields
							if (!ok)										// if we did not
							{												// then tell the user, especially the ones who hardwired the number of columns
								nxLog::Record(NXLOG_WARNING, "nxLinearArray::InputColumnMajorText, error reading stream, inconsistent number of fields on line %d, expected %d fields but detected %d", (int)lineno, (int)numx, (int)xctr );
							}												// and that is that
							else											// otherwise
							{												// we have the proper number of columns
								iter = begin();
								iter += iy;									// so describe whic row we are going to parse
								if (m_rankspecs.Rank() == 1) stride = 1;
								else                         stride = m_rankspecs.Dims()[0];

								NXTRACE_ONCEONLY(firsttime,("***** WARNING **** nxArrayLinear<T>::InputColumnMajorText has changed, we need to walk through this code to make sure it works ****** WARNING ****\n"));
								std::string	thisline;
								thisline.assign( oneline.DangerousTypecast(), oneline.GetLength() );

								std::stringstream	istr( thisline, std::ios_base::in );	// attach the current line to a string stream
								for (size_t ix = 0; ix < numx; ix++)		// for each column on this line
								{											// simply
									istr >> (*iter);						// and use the >> stream operator to read in the element
									iter += stride;
								}											// do all columns on this row
								iy++;										// point to the next row
							}
						}
					}
				}
			}
		}
	}
	if (ok && (numdatalines != 0))
	{
		ok = (numy == numdatalines);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "nxLinearArray::InputColumnMajorText, The user requested %d lines of text but I have read in %d lines from the file. The returned array is a (%d x %d) array.", (int)numdatalines, (int)numy, (int)numy, (int)numx);
		}
	}

	if (!ok) erase();
	return ok;
}

template <class T>
bool nxArrayLinear<T>::InputColumnMajorText( const char * filename, size_t numcolumns, size_t numdatalines )
{
	std::ifstream		strm;
	bool				ok;

	strm.open(filename, std::ios_base::in );
	ok = strm.is_open();
	ok = ok && InputColumnMajorText( strm, numcolumns, numdatalines );
	strm.close();
	return ok;

}

/*-----------------------------------------------------------------------------
 *					input stream operator <<		2009-6-18*/
/** Output the array to a stream. Reads the array in one element at a time
 *	using the iterator to set the values **/
/*---------------------------------------------------------------------------*/

template <class T>
std::ostream& operator << ( std::ostream& outputstream, nxArrayLinear<T>& userdata )
{
	nxArrayIter<T>		iter;
	nxArrayIter<T>		finish;
	T					value;

	finish = userdata.end();
	for( iter = userdata.begin(); !(iter == finish); ++iter)
	{
		value = *iter;
		outputstream << value;
	}
	return outputstream;
}

/*-----------------------------------------------------------------------------
 *					input stream operator >>		2009-6-18*/
/** INput the array from a stream. Reads the array in one element at a time
 *	using the iterator to set the values **/
/*---------------------------------------------------------------------------*/

template <class T>
std::istream& operator >> ( std::istream& inputstream, nxArrayLinear<T>& userdata )
{
	nxArrayIter<T>		iter;
	nxArrayIter<T>		finish;
	T					value;

	finish = userdata.end();
	for( iter = userdata.begin(); !(iter == finish); ++iter)
	{
		inputstream >> value;
		*iter = value;
	}
	return inputstream;
}


/*---------------------------------------------------------------------------
 *'					nx1dArray<T>::Indgen							2003-9-3
 *-------------------------------------------------------------------------*/

/* ************ NEEDS WORK. THIS METHOD SHOULD BE OUTSIDE OF THE CLASS AS A STANDALONE FUNCTION. */

template< class T>
nx1dArray<T>& nx1dArray<T>::Indgen( size_t nx)
{
	SetSize(nx);
	T value = 0;
	nxArrayIter<T>	iter   = this->begin();
	nxArrayIter<T>	finish = this->end();

	while (!(iter == finish))
	{
		*iter = value++;
		++iter;
	}
	return *this;
}


template <class T>
bool	nx1dArray<T>::Slice( size_t lo, size_t hi, nxArrayLinear<T>* column) const
{
	return nxArrayLinear<T>::Slice( &lo, &hi, 1, column );
}

template <class T>
nx1dArray<T> nx1dArray<T>::Slice( size_t lo, size_t hi ) const
{
	nx1dArray<T>	slice;

	Slice( lo, hi, &slice);
	return slice;
}

template <class T>
bool nx2dArray<T>::Slice( size_t lo0, size_t hi0, size_t lo1, size_t hi1, nxArrayLinear<T>* column) const
{
	size_t	lo[2] = {lo0, lo1};
	size_t	hi[2] = {hi0, hi1};
	return nxArrayLinear<T>::Slice( lo, hi, 2, column );
}

template <class T>
nx2dArray<T> nx2dArray<T>::Slice(size_t lo0, size_t hi0, size_t lo1, size_t hi1 ) const
{
	nx2dArray<T>	slice;

	Slice( lo0, hi0, lo1, hi1, &slice);
	return slice;
}

template <class T>
bool nx3dArray<T>::Slice( size_t lo0, size_t hi0, size_t lo1, size_t hi1, size_t lo2, size_t hi2, nxArrayLinear<T>* column) const
{
	size_t	lo[3] = {lo0, lo1, lo2};
	size_t	hi[3] = {hi0, hi1, hi2};
	return nxArrayLinear<T>::Slice( lo, hi, 3, column );
}

template <class T>
nx3dArray<T> nx3dArray<T>::Slice(size_t lo0, size_t hi0, size_t lo1, size_t hi1, size_t lo2, size_t hi2 )const
{
	nx3dArray<T>	slice;

	Slice( lo0, hi0, lo1, hi1, lo2, hi2, &slice);
	return slice;
}

template <class T>
bool nx1dArray<T>::WriteToTextFile( const char* filename, bool columnmajor, const char* formatstr )
{
	nxFile	f;
	bool	ok;
	size_t	ix;
	T		value;

	f.Open( filename, "wt");
	ok = f.IsOpen();
	if (ok)
	{
		if (columnmajor)
		{
			for (ix = 0; ix < this->size(); ix++)
			{
				value = At( ix );
				fprintf( f, formatstr, value );
				fprintf( f, " " );
			}
		}
		else
		{
			for (ix = 0; ix < this->size(); ix++)
			{
				value = At( ix );
				fprintf( f, formatstr, value );
				fprintf( f, "\n" );
			}
		}
		f.Close();
	}
	return ok;
}


template <class T>
bool nx2dArray<T>::WriteToTextFile( const char* filename, bool columnmajor, const char* formatstr )
{
	nxFile	f;
	bool	ok;
	size_t	ix;
	size_t	iy;
	T		value;

	f.Open( filename, "wt");
	ok = f.IsOpen();
	if (ok)
	{
		if (columnmajor)
		{
			for (ix = 0; ix < XSize(); ix++)
			{
				for (iy = 0; iy < YSize(); iy++)
				{
					value = At( ix,iy);
					fprintf( f, formatstr, value );
					fprintf( f, " " );
				}
				fprintf( f, "\n" );
			}
		}
		else
		{
			for (iy = 0; iy < YSize(); iy++)
			{
				for (ix = 0; ix < XSize(); ix++)
				{
					value = At( ix,iy);
					fprintf( f, formatstr, value );
					fprintf( f, " " );
				}
				fprintf( f, "\n" );
			}
		}
		f.Close();
	}
	return ok;
}


