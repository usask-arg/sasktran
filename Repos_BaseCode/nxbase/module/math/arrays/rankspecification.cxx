#include "nxbase_math.h"

#if defined (NXDEBUG)
size_t g_nxMemoryManager_numinstance = 0;
nxBOOL g_nxMemoryManager_isinitialized = false;

static void nxMemoryManagerExitCode()
{
#if defined(NXDEBUG)
	if (g_nxMemoryManager_numinstance == 0)
	{
		NXTRACE( ("g_nxMemoryManager Numinstances is not zero\n") );
	}
#endif
}

void nxMemoryManagerInitCode()
{
	if (!g_nxMemoryManager_isinitialized)
	{
		atexit(nxMemoryManagerExitCode);
		g_nxMemoryManager_numinstance   = 0;
		g_nxMemoryManager_isinitialized = true;
	}
}
#endif


/*---------------------------------------------------------------------------
 *'					RankSpecification::RankSpecification             2002-9-1*/
/**	Default constructor */
/*-------------------------------------------------------------------------*/

RankSpecification::RankSpecification()
{
	init();
}


/*-----------------------------------------------------------------------------
 *					RankSpecification::init		2006-1-18*/
/** **/
/*---------------------------------------------------------------------------*/

void RankSpecification::init()
{
	m_flags         = 0;
	m_rank       = 0;
	m_nelements  = 0;
	m_dims       = &m_staticdims[0];
	m_strides    = &m_staticstrides[0];
	*m_dims      = 0;
	*m_strides   = 0;
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::ReleaseResources              2002-9-1*/
/**	Release any memory allocated within this class. */
/*-------------------------------------------------------------------------*/

void RankSpecification::ReleaseResources()
{
	if (m_rank > 0)
	{
		if (m_dims != &m_staticdims[0] )
		{
			if( m_dims != NULL) delete [] m_dims;
			m_dims  = &m_staticdims[0];
		}

		if ( m_strides != &m_staticstrides[0] )
		{
			if (m_staticstrides != NULL) delete [] m_strides;
			m_strides = &m_staticstrides[0];
		}
		m_flags = 0;
		m_rank     = 0;
		*m_dims    = 0;
		*m_strides = 0;
	}
	m_nelements = 0;
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::ReleaseResources              2002-9-1*/
/**	Erase this array layout descriptor. The described array is empty at
 *	the end of this procedure
 */
/*-------------------------------------------------------------------------*/

void RankSpecification::Erase()
{
	m_rank         = 0;
	*m_dims        = 0;
	*m_strides     = 0;
	m_nelements    = 0;
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::~RankSpecification            2002-9-1*/
/** Default Destructor */
/*-------------------------------------------------------------------------*/

RankSpecification::~RankSpecification()
{
	ReleaseResources();
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::AllocateSpace                 2002-9-1*/
/**	Allocate space to hold the dimensions of this array layout descriptor.
 *	If the rank of the array is less than #NX_RANKSPEC_MAXDIM then it will
 *	use internal class storage while if it is greate r it will dynamically
 *	allocate space on the heap.
 *-------------------------------------------------------------------------*/

nxBOOL RankSpecification::AllocateSpace()
{
	nxBOOL	ok = true;

	if (m_rank > NX_RANKSPEC_MAXDIM)
	{
		m_dims    = new size_t   [m_rank];
		m_strides = new size_t[m_rank];
		ok = (m_dims != NULL) && ( m_strides != NULL);
		if (!ok)
		{
			nxLog::Record(NXLOG_ERROR, "RankSpecification::AllocateSpace, Error allocating space for rank (%d) layout", (int)m_rank);
			m_dims    = &m_staticdims[0];
			m_strides = &m_staticstrides[0];
			m_rank    = 0;
		}
	}
	else
	{
		m_dims    = &m_staticdims[0];
		m_strides = &m_staticstrides[0];
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::GetContiguousStorageSize                  2003-6-5*/
/**	Returns the amount of memory in bytes to be allocated for this array as the
 *	number of base elements.  This is normally only called from InxMemoryManager
 *	and that is only called to create contiguous chunks.
*/
/*-------------------------------------------------------------------------*/

size_t	RankSpecification::GetContiguousStorageSize(size_t elementsize) const
{
	size_t	nelem = 0;

	if (m_rank > 0 && elementsize > 0)
	{
		nelem = (m_strides[m_rank-1]*m_dims[m_rank-1] + elementsize -1)/elementsize;
	}
	return nelem;
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::CopyOther		2003-9-15 */
/**	Part of the copy constructor.  resources have already been released
 *	by the caller, so just copy the other guy over to here.
*/
/*-------------------------------------------------------------------------*/

bool RankSpecification::CopyOther( const RankSpecification& other )
{
	nxBOOL		ok;

	ok = (&other  == this);
	if (!ok)
	{
		ReleaseResources();
		m_rank          = other.m_rank;
		m_nelements     = other.m_nelements;
		m_flags         = other.m_flags;
		ok = AllocateSpace();
		if (!ok)
		{
			ReleaseResources();
		}
		else
		{
			for (size_t i =0; i < m_rank; i++)
			{
				m_dims[i]     = other.m_dims[i];
				m_strides[i]  = other.m_strides[i];
			}
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::RankSpecification             2002-9-1*/
/** default copy cponstructor */
/*-------------------------------------------------------------------------*/

RankSpecification::RankSpecification( const RankSpecification& other )
{
	init();
	CopyOther(other);
}


/*---------------------------------------------------------------------------
 *'					RankSpecification::operator=                     2002-9-1*/
/** The assignement operator. Copy other layout descriptor to this instance. */
/*-------------------------------------------------------------------------*/

RankSpecification&	RankSpecification::operator= (const RankSpecification& other)
{
	CopyOther( other );
	return *this;
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::AllocateUninitializedSpace    2003-6-7*/
/**	This routine is created specifically for the slicing and sub-setting
 *	operations that want to create a rank specification for a smaller dimension
 *	and its easier to get this code to preallocate the arrays rather than try
 *	to do do it  from the caller.
 *
 *	The user must fill out the Dims and Strides structures after calling this
 *	guy and then call UpdateStorageSpecs
*/
/*-------------------------------------------------------------------------*/

nxBOOL RankSpecification::AllocateUninitializedSpace ( size_t rank )
{
	nxBOOL ok;

	m_flags      = 0;								// Clear the allocation flags in this
	ok = (rank == m_rank);							// If the rank has changed
	if (!ok)
	{												// then
		if (rank > 2000000000)						// make sure its positive
		{
			nxLog::Record(NXLOG_WARNING, "RankSpecification::Configure, User has requested an amazingly large rank (%d), it has been set to 0", rank);
			rank = 0;
		}											// and
		ReleaseResources();							// release resources from previous rank
		m_rank = rank;								// set up the rank of this object
		ok = AllocateSpace();						// and allocate internal m_dim and m_strides arrays if necessary
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::Configure                     2002-9-1*/
/** Configure the RankSpecification according to the specified dimensions.
 *	This is the only way to define the memory layout of a RankSpecification
*/
/*-------------------------------------------------------------------------*/

nxBOOL RankSpecification::Configure ( size_t rank, const size_t* dims, size_t elementsize, const size_t* strides  )
{
	nxBOOL		iscontiguous;
	nxBOOL		isfixedstride;
	size_t	defaultstride;
	size_t	fixedstride;
	nxBOOL	ok = true;

	ok = AllocateUninitializedSpace( rank);						// Allocate uninitialized data storage areas
	if (ok)														// and if that worked then continue
	{															// so
		if (rank < 1)
		{
			Erase();
		}
		else
		{

			m_nelements   = 1;										// reset the number of elements
			defaultstride = elementsize;							// IF the data are contiguous then the 1st stride equals the lement size.
			if (strides == NULL )
			{
				for (size_t i =0; i < m_rank; i++)							// Work out the strides for each
				{														// dimension of the array
					m_dims[i]      = dims[i];							// copy over the dimension size
					m_strides[i]   = defaultstride;						// If we have no strides passed from user
					defaultstride *= dims[i];							// and get size of next contiguous stride in bytes
					m_nelements   *= m_dims[i];							// get the total number of elements
				}
				SetBit( m_flags, (RS_ISCONTIG | RS_ISFIXED)	);				// Flag that the array is contiguous and fixed stride
			}
			else														// User has passed in strides
			{															// so now we have to check
				iscontiguous      = true;
				isfixedstride     = true;
				fixedstride       = *strides;
				for (size_t i =0; i < m_rank; i++)							// Work out the strides for each
				{														// dimension of the array
					m_dims[i]      = dims[i];							// copy over the dimension size
					m_strides[i]   = strides[i];						// If we have no strides passed from user
					isfixedstride  = isfixedstride && (fixedstride   == m_strides[i] );
					iscontiguous   = iscontiguous  && (defaultstride == m_strides[i] );
					defaultstride  *= m_dims[i];						// and get size of dimension assuming points are uniformly spaced
					fixedstride    *= m_dims[i];						// and get size of dimension assuming points are uniformly spaced
					m_nelements    *= m_dims[i];						// get the total number of elements
				}
				if (iscontiguous  ) SetBit( m_flags, RS_ISCONTIG	);
				if (isfixedstride ) SetBit( m_flags, RS_ISFIXED	);
			}
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *					RankSpecification::LogicalToPointer		2003-9-13 */
/**	Converts a logical (ie linear elemental position) into a pointer to the
 *	location.
*/
/*-------------------------------------------------------------------------*/

nxBYTE*	RankSpecification::LogicalToPointer( size_t logicalposition, nxBYTE* baseptr) const
{
	size_t		offset = 0;
	size_t		x;
	size_t		j;
	size_t		nelem;

	NXASSERT( ( logicalposition >= 0) && (logicalposition <= (size_t)m_nelements) );
	if (m_rank > 0)
	{
		j     = m_rank-1;
		nelem = m_nelements/m_dims[j];							// Get the total number of elements that make up one element of last dimension
		while (j > 0)											// While the rank is above 0
		{														// then
			x                = logicalposition / nelem;			// get the number of digits in this row
			offset          += x*m_strides[j];					// Add this number of strides to the offset
			logicalposition -= x*nelem;							// remove this number of positions from the logical position
			j--;												// step down to the next dimesnion
			nelem           /= m_dims[j];						// get the number of elements in this dimension
		}														// repeat until we done all except lower dimension
	}
	NXASSERT(( logicalposition < m_dims[0] ));
	offset += logicalposition*m_strides[0];					// now do the lowest dimension.
	return (baseptr + offset);								// return the pointer to this element
}


/*---------------------------------------------------------------------------
 *'					RankSpecification::operator==                  2002-10-22*/
/**	returns TRUE if the dimensions of this array descriptor equal the
 *	dimensions of the other array.  Note that no consideration is given to
 *	strides of each dimension.
*/
/*-------------------------------------------------------------------------*/

bool RankSpecification::operator==( const RankSpecification& other ) const
{
	nxBOOL ok;

	ok =    (m_nelements == other.m_nelements)
		 && (m_rank      == other.m_rank);
	if (ok)
	{
		for (size_t i =0; i < m_rank; i++) ok = ok && (m_dims[i] == other.m_dims[i]);
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					RankSpecification::IsSameLayout            2002-10-22 */
/** Return TRUE if this array has the same layout as the requested number of elements
 *	the input rank and dimensions.
 *-------------------------------------------------------------------------*/

bool RankSpecification::IsSameLayout( size_t nrank, const size_t* dims, const size_t* strides  ) const
{
	nxBOOL	ok = true;
	size_t		lowestrank;
	size_t		i;

	lowestrank = (nrank < m_rank) ? nrank : m_rank;				// Get the lowest rank o fthe two elements
	ok = (dims != NULL) && (nrank > 0) && (N_Elements() > 0);	// Make sure we have something to check
	if (strides == NULL) ok = ok && IsContiguous();				// If the user is not passing in strides then this spec must be contiguous.
	if (ok)														// If we have passed the basic tests
	{															// then check each dimension
		for (i =0; ok && (i < lowestrank); i++)					// Check the ranks which are common to both sets
		{														// For each common dimension
			ok = ok && (m_dims[i] == dims[i]);					// Check that they are the same size
			if (strides != NULL)								// and if the user has given us strides
			{													// then
				ok = ok && (m_strides[i] == strides[i]);		// check that the strides are good
			}													// thats that
		}														// do all of the common dimensions

		for (i = lowestrank; i < m_rank; i++)					// Now check the higher dimensions (if any) of this sepc
		{														// all of these
			ok = ok && (m_dims[i] == 1);						// dimensions must be 1 to be compatible
		}

		for (i = lowestrank; i < nrank; i++)					// Now check the higher dimensions (if any) of the other spec
		{														// all of these
			ok = ok && (dims[i] == 1);							// dimesnions must be 1 for the two arrays to have the same layout
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					RankSpecification::CopyMismatchedDimensions		2005-8-5*/
/** **/
/*---------------------------------------------------------------------------*/

nxBOOL RankSpecification::CopyMismatchedDimensions( size_t mandatoryrank, const RankSpecification& source)
{
	size_t			m1;
	size_t		ns;
	size_t			srank;
	size_t			i;

	srank = source.Rank();
	m1     = mandatoryrank-1;
	if (mandatoryrank < srank)							// If we are making an array
	{													// smaller than the source array
		for (i=0;  i < m1;  i++)						// Then copy all but the last
		{												// dimension and stride
			m_dims[i]     = source.m_dims[i];			// as there is no reason
		}												// do all of the lower dimensions

		ns=1;
		for (i=m1; i < srank; i++)						// last dimension of this array is the sum of dimensions from source
		{												// so add up
			ns  *= source.m_dims[i];					// all of the extended dimensions
		}
		m_dims[m1]    = ns;
	}
	else														// otherwise we are making an array
	{															// that is bigger than the source
		for (i=0;  i < srank;  i++)								// so simply
		{														// copy
			m_dims[i]     = source.m_dims[i];					// all of the dimension
		}														// from the source.

		for (i=srank; i < mandatoryrank; i++)					// last dimension of this array is the sum of dimensions from source
		{														// so add up
			m_dims[i]    = 1;									// all of the extended dimensions
		}
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					RankSpecification::CopyMismatchedStrides		2005-8-8*/
/** **/
/*---------------------------------------------------------------------------*/

nxBOOL RankSpecification::CopyMismatchedStrides( size_t mandatoryrank, const RankSpecification& source)
{
	size_t			m1;
	size_t			srank;
	size_t			i;
//	intptr_t	str;
	size_t		str;
	nxBOOL		ok = true;

	srank   = source.Rank();
	m_flags = source.m_flags;							// Copy the layout flags over. Our new ShallowCopy array cannot have different layout flags.
	m1      = mandatoryrank-1;								// so
	if (mandatoryrank < srank)							// If we are making an array
	{													// smaller than the source array
		for (i=0;  i < m1;  i++)						// Then copy all but the last
		{												// dimension and stride
			m_strides[i]  = source.m_strides[i];		// these cannot be the same
		}												// do all of the lower dimensions
		m_strides[m1] = source.m_strides[m1];
		str           = m_strides[m1];
		for (i=m1; i < srank; i++)
		{
			ok = ok && (str == source.m_strides[i]);	// Make sure the higher dimension strides are all compatible with one big "fixed or contiguous" stride.
			str *= source.m_dims[i];
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "RankSpecification::ReshapeToMandatoryRank, Cannot reshape array to lower rank as upper dimensions do not have proper fixed stride");
		}
	}
	else														// otherwise we are making an array
	{															// that is bigger than the source
		for (i=0;  i < srank;  i++)								// so simply
		{														// copy
			m_strides[i]  = source.m_strides[i];				// and strides
		}														// from the source.
		str = source.m_strides[srank-1]*source.m_dims[srank-1];	// get the last stride of the "next" dimension from the source
		for (i=srank; i < mandatoryrank; i++)					// last dimension of this array is the sum of dimensions from source
		{														// so add up
			m_strides[i] = str;									// and there is no problem with strides "not matching"
		}
	}
	if (!ok) Erase();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					RankSpecification::ReshapeToMandatoryRank		2005-8-5*/
/** This makes a new rank object for a destination array based upon the
 *	current contents of the source rank.  This function
 *	is normally used by ShallowCopy. It is trivial when the source and destination
 *	have the same rank (its just a complete copy of dimensions and strides) but is
 *	more complicated when the source and destination have different dimensionality.
 *	The shallow copy operation has to consider two aspects to consider.
 *	First the dimensions have to be made compatible. This is always possible by either
 *	combining the higher dimensions into 1 big dimension or by adding "unit"
 *	dimensions onto the end of an array. Second the shallow copy must ensure that the strides of the two arrays are
 *	still compatible.  This is not always possible
 */
/*---------------------------------------------------------------------------*/

nxBOOL RankSpecification::ReshapeToMandatoryRank( nxBOOL doshallowcopy, int mandatoryrank, size_t elementsize, const RankSpecification& source )
{
	nxBOOL				ok;
	const size_t*		strideptr = NULL;

	ok =    (mandatoryrank < 1 )							// Trivial if no restrictions on destination rank
		 || (source.Rank() == (size_t)mandatoryrank)				// Trivial if ranks are equal
		 || (source.Rank() == 0);							// Trivial if source is empty.
	if (ok)													// If we have a trivial case
	{														// then
		ok = CopyOther( source );							// and copy the source rank (dims, strides and flags)
	}														// and that is all we have to do for the trivial cases
	else														// otherwise
	{															// we have non-matching ranks
		RankSpecification rankbuffer(*this);

		ok = rankbuffer.CopyMismatchedDimensions( mandatoryrank, source);
		if (doshallowcopy)
		{
			ok = ok && rankbuffer.CopyMismatchedStrides( mandatoryrank, source);	// so copy over the flags and the strides
			strideptr = rankbuffer.Strides();
		}
		ok = ok && Configure( mandatoryrank, rankbuffer.Dims(), elementsize, strideptr);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					RankSpecification::Reshape		2005-8-8*/
/** Reshapes the current RankSpecification to the new **/
/*---------------------------------------------------------------------------*/

nxBOOL RankSpecification::Reshape( size_t nrank, const size_t* dims )
{
	nxBOOL      ok;
	size_t	    npts;
	intptr_t	stride;
	size_t			i;

	ok = (IsContiguous() || IsFixedStride()) && (nrank > 0);
	if (ok)
	{
		npts = 1;
		for (i=0; i < nrank; i++) npts *= dims[i];
		ok = (npts == N_Elements());
		if (ok)
		{
			stride = Stride();
			ok = AllocateUninitializedSpace(nrank);
			if (ok)
			{
				ok = Configure( nrank, dims, (int)stride, NULL );
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "RankSpecification::Reshape, abnormal memory re-shaping error, erasing whole rank-specification. This may cause knock on problems");
				Erase();
			}
		}
	}
	return ok;
}
