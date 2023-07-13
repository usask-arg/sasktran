#if !defined(NXBASE_NXARRAYLIB_H)
#define NXBASE_NXARRAYLIB_H

#include "../../system/debug/nxtrace.h"
#include "../../system/nxbitmanip.h"
//#include "nxbase_nxlog.h"

/** \defgroup nxLinearArray_Internals Internal Classes
 *  \ingroup nxLinearArray
 */

/** \ingroup nxLinearArray
	\page Page1 Outline
	\par Introduction
    This is a C++ template library for handling and processing of multi-dimensional arrays. The
    purpose of the library is to provide programmers with array handling capabilities seen
    in other languages like IDL, Matlab and NumPy. I have strived to make the code reasonably efficient
	in both memory and speed while imposing no constraints on the size, rank and layout of the
	array. I have defined an array as having N dimensions defined by the rank.

	\par What is an Array?
	Arrays mean different things to different people. The arrays supported in this code are
	grids of numerical objects, where each individual element supports normal mathematical
	operations. So \b float, \b double, \b int and \b complex are sensible types for array while
	general purpose structures and classes are unlikely to work. Arrays of arbitrary
	dimension can be defined through nxArrayLinear while convenient wrappers for 1, 2 and 3
	dimensional arrays have been provided in nx1dArray, nx2dArray and nx3dArray. Each dimension
	of an array consists of N elements separated in memory by a stride distance. The stride distance
	for each dimension can be chosen (and usually is) so that all of the elements from all dimensions are
	contiguous in memory. However, it is useful on occassion to have a fixed stride distance that is
	not consistent with a contiguous array, for example a vertucal slice through an image has
	elements separated by one row of elements rather than one element. In addition, arrays derived from
	arrays of structures have non-standard strides.  This implementation provides support for three
	arranagements of memory
	  - Elements of the entire array are contiguous in memory
	  - Elements of the entire array are a fixed stride.
	  - Each dimension has its own unique stride between successive elements of that dimension.

	\par Arrays and STL
	The arrays are compatible with the C++ Standard Template Library and implement
	intelligent iterators that are compatible with the STL random access iterators.
	The iterator is from class nxArrayIter<T> and provides a programmatic way to step through
	all the elements of an array.  The nxArrayIter class has been written so it can take advantage
	of the memory layout of the underlying array: indexing of contiguous arrays is easier than
	indexing fixed strides and both are easier than indexing variable stride arrays.  The important point
	is that the nxArrayIter class can be used with any STL algorithm that accepts random access iterators (eg std::sort)

	\par Logical Indexing

	general purpose base class #nxArrayLinear which can be used for arrays of arbitrary rank while classes #nx1dArray,
	#nx2dArray and #nx3dArray are provided as convenient 1-D, 2-D and 3-D wrappers of the base
	class.
	The library predefines classes #nx1dArray, #nx2dArray and #nx3dArray  which most
	programmers will use for their everyday needs.  Rather than explain the classes a
	simple example is given below.

	\par Example

	\code
	bool myprocedure( nx2dArray<double>* image)
	{
		nx1dArray<int>			tau;
		int						i,j;
		bool					ok;

		ok =       tau.SetSize( 100);					// Set the size of tau to 100 elements
		ok = ok && image->SetSize( 100, 45 );			// Set the size of umage to 100 x 45 (the 100 elements are contiguous in memory, not the 45 elements)
		if (ok)											// if it all worked ok
		{												// then
			image->SetTo(0.0);							// Set each element of image to 0 (just for good measure)
			for (i =0; i < 100; i++)					// now loop throug the lements
			{											// first
				tau.At(i) = i;							// set tau to this loop index
				for (j = 0; j < 45; j++)				// now loop over the second index
				{										// so
					image->At(i,j) = tau.At(i) + j*j;	// set each element o fthe image
				}										// do all of the images second index
			}											// do all of the images first index
		}												// if not ok then return immediately
		return ok;										// and return the status to the user.
	}
	\endcode

	The templated nature of the class is readily apparent as we have an array of \c int
	in variable \c tau and an array of \c double in variable \c image. The array classes manage their
	own internal storage and look after memory disposal when the array goes out of scope. The user is not burdened with
	remebering to deallocate the internal array storage. I have implemented the array indexing using the \c At method.
	The astute reader will note that the At method for the 2D case specifies the X then the Y element which is opposite to the
	way one would use the C [] operator. This is intentional as we feel it is much more natural to specify the indexing in the same
	order as the dimensions. The [] operator is still supported by #nx1dArray, #nx2dArray and #nx3dArray and behaves the same way
	as normal C.

	For the purposes of this class, a multidimensional array is an area of memory that has an
	organized layout specified by the #RankSpecification class.  An array
	consists of N dimensions (rank N). Each dimension has a fixed number of elements, (the dimension).
	The spacing between consecutive elements of each dimension is a fixed number of bytes (the stride).
	The first element of each dimension may be offset by a number of bytes from the start of the memory
	area (the offset). With this definition we can easily define normal 1, 2 and 3 dimensional arrays
	that are contiguous in memory.  In addition we can define arrays which are not contiguous in memory.
	An example is an array generated from a single field of an array of structures while another is
	slicing a 2-D image out of a 3-D structure.

	\par Memory Management
	The array classes generally look after their own memory management. The class will normally allocate
	and deallocate memory as required and will ensure and allocated memory is deallocated during
	its destruction.  Memory allocation is performed by the #InxMemoryManager class: the default
	implementation simply uses the C++ new and delete operators.  However the user is free to use other memory
	managers (eg. IDL, MATLAB, Python, Basic SAFEARRAY) and hook these memory managers into a given class.
	Since the memory managers only affect how memory is allocated they have no impact on all of the other
	methods developed for analysis and processing of arrays.

	\par Array Rank
	A key element of the library is the concept of array rank. We use rank to mean the number of dimensions
	The rank of an array, within the terms of this library,
	is implemented as class #RankSpecification which describes how the array is organized in memory.
	This allows the library through the use of strides and offsets for each dimension to handle arrays
	which are not contiguous in memory. It is expected the main use for this capability	is to extract a
	single field from an array of structures or to slice multi-dimensional images into smaller sub-sets.

	\par Array Indexing
	The classes have been written so they use the At method to index various elements.  The At method always
	indexes the array elements starting with the first dimension (most rapidly varying) and finishing with the
	last index.  This is the opposite order to classic "C" [] operator which specifies the last dimension first.
	So b.At(x,y) is equivalent to b[y][x].  The #nx1dArray, #nx2dArray and nx3dArray classes support the []
	operator in the same sense as the C operator but it is quite inefficient.  The At operator is much more
	efficient with code paths optimized for contiguous and fixed stride layouts.

	\par Array Iteration
	Users who wish to process arrays as a contiguous series of elements and have little concern about the
	array dimensionality can use the array iterator #nxArrayIterator.  These iterators are compatible with
	the Standard Template Library random access iterator and can be used to move to any element within
	the array.  The processing order starts with the first dimension  and its way
	through all of the other dimensions in ascending order: this is true regardless of the physical
	layout (strides and offsets) of the array.  The iterators can be used in STL algorithms for searching
	and sorting although they do not work for insertion and deletion.

	\par Operator Overloading
	The nxArray classes have had most of the mathematical operators overloaded so that entire array arithmetic
	can be performed in single C++ statements.  The operators *=, +=, -= and /= are very efficiently
	implemented and occur almost no extra overhead. The *,+,-+ operators are not as optimal as they require
	the creation of an intermediate temporary array (however they are very convenient and quite often worth the
	price). One implication of this is that the template objects used in the array should be an object that
	support all of the normal mathematical operations.

	\par Template Types:
 	- \b T.  The array data type. Usually a simple mathematical data type like \c int or \c double.
	.

	\par Not Yet Implemented
	The biggest omission from the classes is that there are no I/O routines .

*/

#include <iterator>
#include <deque>

//using namespace std;


/*-----------------------------------------------------------------------------
 *					ContainerIsAscendingOrder		 2014- 4- 30*/
/** **/
/*---------------------------------------------------------------------------*/

template< class ITERTYPE>
bool ContainerIsAscendingOrder( ITERTYPE first, ITERTYPE last)
{
	bool		ok = true;
	ITERTYPE	next(first);

	++next;
	while (ok && !(next == last))
	{
		ok = ok && (*next >= *first );
		++first;
		++next;
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *
 *-------------------------------------------------------------------------*/


#define NX_RANKSPEC_MAXDIM	3											// The number of dimensions statically cached in the RankSpecification object

/*--------------------------------------------------------------------------
 *					class RankSpecification                        2002-9-1*/
/**	\ingroup nxLinearArray_Internals
 *	Class used to describe the shape and physical layoutof an array. The class
 *	holds the rank and dimension specification details. The class
 *	also manages calculating indexing of zero based indexed elements from the
 *	start of the memory region.
 */
/*------------------------------------------------------------------------*/

class RankSpecification
{
	private:
		enum MemoryLayout { RS_UNDEFINED = 0, RS_ISCONTIG = 1, RS_ISFIXED = 2};

	private:
		nxWORD					m_flags;
		size_t					m_rank;									//!< The rank or number of dimensions
		size_t					m_nelements;							//!< The number of individual elements in the array
		size_t*					m_dims;									//!< Pointer to the array of dimension sizes
		size_t*					m_strides;								//!< Pointer to the array of strides. Distance between each element of this dimension in bytes.
		size_t					m_staticdims   [NX_RANKSPEC_MAXDIM];	//!< Upto 4 dimensions and 4 strides are statically allocated
		size_t					m_staticstrides[NX_RANKSPEC_MAXDIM];	//!< as this keeps the number of memory allocations down

	private:
		void					init						();
		bool					AllocateSpace				();
		bool					CopyOther					( const RankSpecification& other );
		void					ReleaseResources			();
		bool					CopyMismatchedDimensions	( size_t mandatoryrank, const RankSpecification& source);
		bool					CopyMismatchedStrides		( size_t mandatoryrank, const RankSpecification& source);

	public:
								RankSpecification			();
							   ~RankSpecification			();
								RankSpecification			( const RankSpecification& other );
		bool					AllocateUninitializedSpace	( size_t rank);					//!< Allocates space to hold the descriptors for the specifued array
		RankSpecification&		operator=					( const RankSpecification& other );
		bool					operator==					( const RankSpecification& other ) const;
		bool					Configure					( size_t rank, const size_t* dims, size_t elementsize, const size_t* strides );	//!< Configue the array layout descriptor
		void					Erase						();																			//!< Erase this array layout descriptor
		nxBYTE*					LogicalToPointer			( size_t logicalposition, nxBYTE* baseptr) const;						//!< Convert the logical linear poisition to a pointer to the element.
		bool					IsSameLayout				( size_t nrank, const size_t* dims, const size_t* strides) const;												//!< See if this array descriptor is the same size as the given rank and dimensions
		size_t					GetContiguousStorageSize	( size_t elementsize) const;											//!< Returns the storage size of this array in bytes
		bool					ReshapeToMandatoryRank		(  bool doshallowcopy, int mandatoryrank, size_t elementsize, const RankSpecification& source );


	public:				// ---- Properties
		inline bool				IsContiguous				() const { return BitIsSet (m_flags,RS_ISCONTIG);}	//!< Returns TRUE if the entire array is contiguous in memory
		inline bool				IsFixedStride				() const { return BitIsSet (m_flags,RS_ISFIXED);}	//!< Returns TRUE if each element of every dimension is a fixed stride from the next element
		inline size_t				Rank						() const { return m_rank;}							//!< Returns the Rank of this array
		inline const size_t*	Dims						() const { return m_dims;}							//!< Returns a pointer to the number of elements in each dimension. There are #Rank dimensions.
		inline const size_t*	Strides						() const { return m_strides;}						//!< Returns a pointer to the stride, in bytes of each dimension. There are #Rank dimensions.
		inline size_t*			UnsafeDims					()		 { return m_dims;}							//!< Returns a non constant pointer to the internal dimensions array. This is considered unsafe and should only be used where absolutely necessary
		inline size_t*			UnsafeStrides				()       { return m_strides;}						//!< Returns a non-constant pointer to the internal strides array. This is considered unsafe and should only be used where absolutely necessary
		inline size_t			Stride						() const { return *m_strides;}						//!< Returns the stride of the first dimension. Useful for fixed stride arrays
		inline size_t			N_Elements					() const { return m_nelements;}				//!< Returns the number of elements in this array.
		bool					Reshape						( size_t nrank, const size_t* dims );
};


/*---------------------------------------------------------------------------
 *					InxMemoryManager								2002-5-27*/
/**	\ingroup nxLinearArray_Internals
 *	This is a class that manages linear chunks of memory.  The coding provides
 *	a mechanism that:
 *
 *	- a) allows different memory allocation schemes to
 *	    be used eg. new, malloc, idl arrays, matlab arrays, python arrays,
 *	    COM safearrays and static arrays.  Simply create derived classes
 *		and override the InternalFreemem and InternalAllocate functions.  Note
 *		the functions don't need to do much parameter checking as it has all be done
 *	- b)  provides a mechanism so chunks of allocated memory can be quickly shifted
 *		or shared between array objects.  This helps eliminate the allocation
 *		and deallocation of temporary arrays when executing code with lots of
 *		operator overloading. The Memory managers lifetime is managed by AddRef
 *		and Release (like a COM object, but it is NOT a COM object).
 *
 *	Note that InxmemoryManager objects must be created on the heap before they
 *	can be succesfully shared between different array objects.
 */
 /*------------------------------------------------------------------------*/

#if defined (NXDEBUG)
extern size_t g_nxMemoryManager_numinstance;
extern bool g_nxMemoryManager_isinitialized;
void nxMemoryManagerInitCode();
#endif

template <class T>
class InxMemoryManager
{
	private:
		size_t				m_lockcount;			//!< Number of "locks" on this object. I.E. Level of sharing
		size_t				m_memorylockcount;		//!< Number of locks on memory
		size_t				m_linearsize;			//!< This is the number of elements in the array
		T*					m_pointer;				//!< This is the pointer to the start of memory
		T*					m_endpointer;			//!< This is the pointer to the start of memory
		size_t				m_reservesize;			//!< The amount of memory actually allocated

	private:
		void				Erase();

	protected:
		virtual void		InternalFreemem			();
		virtual T*			InternalAllocate		(  size_t linearsize);
		virtual bool		InternalAllowRealloc	()	{ return true;}

	public:
							InxMemoryManager	();
		virtual			   ~InxMemoryManager	();
		size_t				AddRef				()			{ return ++m_lockcount;}
		size_t				Release				()			{ if (--m_lockcount > 0) return m_lockcount; delete this; return 0;}

	public:
		bool				AllocateAndLock		(  const RankSpecification* rankspecs, bool isallowedtogrow, T** ptr );		// Allocate to this number of elements

	public:
		bool				CheckMemoryRange	( T* baseptr, T* lastptr ) { return (baseptr >= m_pointer) && (baseptr < m_endpointer) && (lastptr >= baseptr) && (lastptr < m_endpointer);}
		ULONG				LockCount			() const	{ return m_lockcount;}
		size_t				LinearSize			() const	{ return m_linearsize;}
		void				LockMemory			()			{ if (m_pointer != NULL) ++m_memorylockcount;}
		void				UnlockMemory		();
		bool				MemoryIsShared		() const	{ return (m_memorylockcount > 1);}

};



/*---------------------------------------------------------------------------
 *'					InxMemoryManagerIDL                             2002-9-11
 *	Here is the basic template for attaching to an IDL array.
 *-------------------------------------------------------------------------*/

template <class T>
class InxMemoryManagerIDL : public InxMemoryManager<T>
{
	protected:
		virtual void		InternalFreemem ();
		virtual T*			InternalAllocate(  const RankSpecification* rankspecs);
		virtual bool		InternalAllowRealloc	()	{ return false;}


	public:
							InxMemoryManagerIDL	()			{ }
		virtual			   ~InxMemoryManagerIDL	()			{ }
};


/*-----------------------------------------------------------------------------
 *					 InxMemoryManagerSafeArray::InxMemoryManager<T>		2009-7-24*/
/** **/
/*---------------------------------------------------------------------------*/

#if defined(NX_WINDOWS)

template <class T>
class InxMemoryManagerSafeArray	: public InxMemoryManager<T>
{
	private:
		SAFEARRAY*			m_safearray;
		VARIANT				m_var;

	protected:
		virtual void		InternalFreemem			();
		virtual T*			InternalAllocate		(  const RankSpecification* rankspecs);
		virtual bool		InternalAllowRealloc	()	{ return false;}


	public:
							InxMemoryManagerSafeArray	();
		virtual			   ~InxMemoryManagerSafeArray	();
};
#endif

/*---------------------------------------------------------------------------
 *'					class nxArrayIterCore                            2003-6-6*/
/**	\ingroup nxLinearArray_Internals
 *	This class is a base class for implementing iterators for the arrays.
 *	It is provides the implementation for iterators of contiguous arrays and
 *	provides a base class for #nxArrayIterfixedStride and #nxArrayIterVariableStride.
*/
/*-------------------------------------------------------------------------*/


template <class T>
class nxArrayIterCore
{
	protected:
		T*						m_pointer;

	public:
								nxArrayIterCore()										{ m_pointer = NULL; }				//m_beginpointer = NULL; m_endpointer = NULL; }
								nxArrayIterCore( const nxArrayIterCore& other )			{ m_pointer = other.m_pointer;}
		nxArrayIterCore<T>&		operator=      ( const nxArrayIterCore& other )			{ m_pointer = other.m_pointer; return *this;}
		virtual					~nxArrayIterCore()										{}

		inline T*				Pointer					() const { return m_pointer;}
		inline T&				operator *				() 										{ return *m_pointer;}
		inline bool				operator <				( const nxArrayIterCore& other) const	{ return m_pointer < other.m_pointer;}
		inline bool				operator >				( const nxArrayIterCore& other) const	{ return m_pointer > other.m_pointer;}
		inline bool				operator <=				( const nxArrayIterCore& other) const	{ return m_pointer <= other.m_pointer;}
		inline bool				operator >=				( const nxArrayIterCore& other) const	{ return m_pointer >= other.m_pointer;}
		inline bool				operator ==				( const nxArrayIterCore& other) const	{ return m_pointer == other.m_pointer;}
		inline bool				operator !=				( const nxArrayIterCore& other) const	{ return m_pointer != other.m_pointer;}

	public:
		virtual void			Configure( nxBYTE* basepointer, const RankSpecification* rank )			{ m_pointer = (T*)basepointer; if (!rank->IsContiguous() )  throw("nxArrayIterCore::Configure, Rank is not contiguous"); }
		virtual void			operator ++				()									{++m_pointer;}					// Pre  increment
		virtual void			operator --				()									{--m_pointer;}					// Pre  decrement
		virtual void			operator +=				(size_t n)							{m_pointer += n; }				// Addition
		virtual void			operator -=				(size_t n)							{m_pointer -= n; }				// Subtraction
		virtual	T&				operator []				(size_t n)							{return m_pointer[n];}			// Indexing
		virtual	size_t			operator -				(const nxArrayIterCore<T>& other)	{return (size_t)(m_pointer-other.m_pointer);}			// Indexing
};


/*<--------------------------------------------------------------------------
 *'					class nxArrayIterFixedStride                     2002-9-7*/
/**	\ingroup nxLinearArray_Internals
 *	Class used to iterate over fixed stride arrays in a quick and efficient
 *	manner.
*/
/*------------------------------------------------------------------------*/

template <class T>
class nxArrayIterFixedStride : public nxArrayIterCore<T>
{
	private:
		size_t					m_stride;

	public:
									nxArrayIterFixedStride	()	{ m_stride = 0;}
		virtual					   ~nxArrayIterFixedStride	( ) {}
								    nxArrayIterFixedStride	(const nxArrayIterFixedStride& other)	{ nxArrayIterCore<T>::operator=(other); m_stride = other.m_stride;}
		nxArrayIterFixedStride&		operator=				(const nxArrayIterFixedStride& other )	{ nxArrayIterCore<T>::operator=(other); m_stride = other.m_stride; return *this;}

	public:
		virtual void			Configure( nxBYTE* basepointer, const RankSpecification* rank ){ this->m_pointer = (T*)basepointer; m_stride  = *(rank->Strides()); NXASSERT( rank->IsFixedStride() );}
		virtual void			operator ++ ()									{ this->m_pointer = (T*)((nxBYTE* )this->m_pointer + m_stride);}	// Pre increment
		virtual void			operator --	()									{ this->m_pointer = (T*)((nxBYTE* )this->m_pointer - m_stride);}	// Pre decrement
		virtual void			operator +=	(size_t n)							{ this->m_pointer = (T*)((nxBYTE* )this->m_pointer + n*m_stride); }
		virtual void			operator -=	(size_t n)							{ this->m_pointer = (T*)((nxBYTE* )this->m_pointer - n*m_stride); }
		virtual	T&				operator []	(size_t n) 							{ T* p =(T*)((nxBYTE* )this->m_pointer + n*m_stride); return *p;}
		virtual	size_t			operator -	(const nxArrayIterCore<T>& other)	{return ((nxBYTE *)this->m_pointer - (nxBYTE *)other.Pointer())/m_stride;}			// Indexing
};

/*<--------------------------------------------------------------------------
 *'					class nxArrayIterVariableStride                  2002-9-7*/
/**	\ingroup nxLinearArray_Internals
 *	Class used to iterate over variable stride arrays. This is quite a bit
 *	slower than fixed stride arrays as there is more book-keeping to track.
 *	The code keeps track of the logicalposition in the array (ie. number of elements
 *	from the beginning) unfortunately every time we update the logical position
 *	we must also update the internal pointer which takes time.
 */
/*------------------------------------------------------------------------*/

template <class T>
class nxArrayIterVariableStride : public nxArrayIterCore<T>
{
	private:
		size_t							m_logicalposition;			// The logical position in the array.
		nxBYTE*							m_basepointer;				// The start of the array
		const RankSpecification*		m_rankspec;					// The memory layout of the array

	private:
		void							UpdatePointer()				{ this->m_pointer = (T*)m_rankspec->LogicalToPointer( m_logicalposition, m_basepointer);}

	public:
										nxArrayIterVariableStride()	{}
									   ~nxArrayIterVariableStride() {}
										nxArrayIterVariableStride	( const nxArrayIterVariableStride<T>& other );
		nxArrayIterVariableStride<T>&	operator =					( const nxArrayIterVariableStride<T>& other );

	public:
		virtual void		Configure		(nxBYTE* basepointer, const RankSpecification* rank);
		virtual void		operator ++		()									{m_logicalposition++;    UpdatePointer();}
		virtual void		operator --		()									{m_logicalposition--;    UpdatePointer();}
		virtual void		operator +=		(size_t n)							{m_logicalposition += n; UpdatePointer();}
		virtual void  		operator -=		(size_t n)							{m_logicalposition -= n; UpdatePointer();}
		virtual	T&			operator []    	(size_t n)							{size_t pos = m_logicalposition; m_logicalposition += n; UpdatePointer(); T* p = this->m_pointer; m_logicalposition = pos; return *p;}
		virtual	size_t		operator -		(const nxArrayIterCore<T>& other)	{const nxArrayIterVariableStride<T>* p = ( const nxArrayIterVariableStride<T>*)&other; return (size_t)(m_logicalposition - p->m_logicalposition);}
};

/*---------------------------------------------------------------------------
 *'					class nxArrayIter                                2002-9-7*/
/**	\ingroup nxLinearArray
 *	Class used to iterate over nxArrayLinear arrays. The iteration handles
 *	contiguous, fixed and variable stride arrays in a quick and efficient manner.
 *	The fixed stride pre-increments are fast (close to standard C pre-increment
 *	speeds) while the variable stride pre-increments have a bit more overhead
 *	and will be somewhat slower.
 *
 *	The class is derived from STL random access iterators and provide
 *	the user with much of the search and iteration algorithms found in the
 *	STL.  Note that we have not made any attempt to support element insertion
 *	or deletion (ie anything that changes the array size)
 *
 */
/*------------------------------------------------------------------------*/

template <class T>
class nxArrayIter : public std::iterator< std::random_access_iterator_tag, T, size_t, T*, T&>
{
	private:
		nxArrayIterCore<T>					m_contiguous;				//!< The very fast contiguous implementation
		nxArrayIterFixedStride<T>			m_fixedstride;				//!< The quite fast fixed stride implementation
		nxArrayIterVariableStride<T>		m_variablestride;			//!< The slower variable stride implementation.
		nxArrayIterCore<T>*					m_iterimp;					//!< Pointer to the implementation used by this instance

//	public:	// ---- Iterator Traits for STL implementation
//		typedef _Category						iterator_category;
//		typedef T								value_type;
//		typedef _Diff_type						difference_type;
//		typedef _Pointer_type					pointer;
//		typedef _Reference_type					reference;

	public:
								nxArrayIter				()	{ m_iterimp = &m_fixedstride;}							//!< Default C++ constructor.
								nxArrayIter				( nxBYTE* basepointer, const RankSpecification* rank);		//!< The constructor called by all of the nxArrayLinear code.
								nxArrayIter				( const nxArrayIter& other):std::iterator< std::random_access_iterator_tag, T, size_t, T*, T&>(other) {*this = other;}	//!< The default copy constructor
		nxArrayIter&			operator =				( const nxArrayIter& other);																				//!< The asignment operator
		T&						operator *				() const						{ return m_iterimp->operator*();}										//!< The dereferencing operator
		nxArrayIter&			operator ++				()								{ m_iterimp->operator++();    return *this;}	//!< The pre-increment operator (fast)
		nxArrayIter&			operator --				()								{ m_iterimp->operator--();    return *this;}	//!< The pre-decrement operator (fast)
		nxArrayIter				operator ++				(int /*d*/)						{ nxArrayIter<T> a(*this); m_iterimp->operator++();    return a;}		//!< The post increment operator (slow). Consider using pre-increment for efficiency
		nxArrayIter				operator --				(int /*d*/)						{ nxArrayIter<T> a(*this); m_iterimp->operator--();    return a;}		//!< The post increment operator (slow). Consider using pre-decrement for efficiency.
		nxArrayIter&			operator +=				(size_t n)						{                          m_iterimp->operator+=(n);   return *this;}	//!< Randomly move forward N elements in the array (fast)
		nxArrayIter&			operator -=				(size_t n)						{                          m_iterimp->operator-=(n);   return *this;}	//!< Randomly move backward N elements in the array (fast)
		nxArrayIter				operator +				(size_t n)						{ nxArrayIter<T> a(*this); a.m_iterimp->operator+=(n); return a;}		//!< Randomly move forward N elements. Slow. Consider using the += operator.
		nxArrayIter				operator -				(size_t n)	const				{ nxArrayIter<T> a(*this); a.m_iterimp->operator-=(n); return a;}		//!< Randomly move backward N elements. Slow. Consider using the -= operator
		size_t					operator -				( const nxArrayIter<T>& other) const { return m_iterimp->operator-( *(other.m_iterimp) );}					//!< Get the number of elements between two iterators.
		T&						operator []				(size_t n) 						{ return m_iterimp->operator[](n);}										//!< Get the element N locations from the current position

		bool					operator <				( const nxArrayIter& other) const	{ return m_iterimp->operator< ( *other.m_iterimp);}						//!< Return TRUE if this iterator is less than the other iterator
		bool					operator >				( const nxArrayIter& other) const	{ return m_iterimp->operator> ( *other.m_iterimp);}						//!< Return TRUE if this iterator is greater than the other iterator
		bool					operator <=				( const nxArrayIter& other) const	{ return m_iterimp->operator<=( *other.m_iterimp);}						//!< Return TRUE if this iterator is less than or equal to the other iterator
		bool					operator >=				( const nxArrayIter& other) const	{ return m_iterimp->operator>=( *other.m_iterimp);}						//!< Return TRUE if this iterator is greater than or equal to the other iterator
		bool					operator ==				( const nxArrayIter& other) const	{ return m_iterimp->operator==( *other.m_iterimp);}						//!< Return TRUE if this iterator equals to the other iterator
		bool					operator !=				( const nxArrayIter& other) const	{ return m_iterimp->operator!=( *other.m_iterimp);}						//!< Return TRUE if this iterator is not equal greater than or equal to the other iterator
};

#define NXARRAY_STARSELECT ((size_t)(-1))
/*-----------------------------------------------------------------------------
 *'					nxArrayLinear									2002-5-27*/
/** \ingroup nxLinearArray
 *	This is the base class for all of the array classes. All arrays can be
 *	considered as linear arrays.
 *
 *	The array paradigm is a collection of dimensions (without dimensional limit)
 *	where the spacing of elements in each dimension is a fixed (byte) stride.
 *  This should allow the base class to extract arrays from arrays of structures
 *	and to efficiently slice and dice existing arrays.
 *
 *	\par Element Types
 *	The array class is templated but I have made the tacit assumption that the basic
 *	elemental types are numeric in nature. In particular I have assumed the mathematical
 *	operators, +,-,*,/, >, < etc. are defined for the elements (type T).
 *
 *	\par Memory Management
 *	The class has moved towards using a memory manager to allocate and deallocate
 *	memory.  This is particularly useful for IDL, python, Matlab and COM
 *	development  as all of these technologies need to use their own memory
 *	allocation routines.
 *
 *  Support is still provided for user managed chunks of memory where the
 *	user just passes in a pointer.  This is particularly useful
 *	for mapping the array onto chunks of memory within  local variables.
 *
 *  The memory manager objects created by one nxArrayLinear instance can be shared
 *	with other instances.  This is the mechanism by which regions of managed memory
 *	can be quickly passed around.  The memory manager has a lifetime controlled by
 *	calls to AddRef and Release. However the memory controlled by the memory manager
 *	also has a lock count controlled by the number of nxArrayLinear instances using
 *	the controlled memory.
 *
 *	\par Dimensions and Strides
 *	The array is specified as dimension "m_rank". Internally the m_dims member stores
 *	the number of elements in each dimension. Member m_stride stores the stride in bytes of
 *	each element of each dimension. The stride is the distance in bytes between
 *	consecutive members of the same dimension. Contiguous arrays set the stride of the
 *	first dimension to sizeof(T) and calculate the stride of subsequent dimensions (N*sizeof(T))
 *
 */
 /*------------------------------------------------------------------------*/

template <class T>
class nxArrayLinear
{
	private:
		typedef T*					(nxArrayLinear::* INDEXTOPOINTER)( const size_t* dimsindex) const;

	private:													// Protected so nxArrayTemporary has direct access
		InxMemoryManager<T>*		m_manager;					// The memory manager controlling this object, NULL if not using a memory manager (ie attached)
		RankSpecification			m_rankspecs;				// The rank, dims and strides of this array.
		INDEXTOPOINTER				m_indextopointer;			// Points to a member function to conver an index to a nxBYTE* pointer
		T*							m_baseptr;					// The pointer to the base location (usually first location) of this array
		T*							m_endptr;					// The pointer to the last  location
		bool						m_checkbounds;				// If true then check indexing bounds when calling At(..)
		bool						m_allowmemoryreuse;			// If true then allow this array to shrink to smaller sizes using current memory allocation (default is false).
		bool						m_arrayisattached;			// True if this array is attached to someone else's memory. The array cannot change size, layout or base pointe. It can only Detach.
		bool						m_istemporary;				// only nxArrayTemporary arrays set this value.

	private:
		bool						IsAttachedToPointerWithoutManager () const { return (m_manager == NULL) && (m_baseptr != NULL);}
		void						UnlockMemoryIfLocked		(bool force);
		void						ReleaseMemoryManager		();
		bool						CheckBounds					( const size_t* dimsindex ) const;
		T*							IndexToPointer_EmptyArray	( const size_t* dimsindex ) const;
		T*							IndexToPointer_General		( const size_t* dimsindex ) const;
		T*							IndexToPointer_1D_Contiguous( const size_t* dimsindex ) const;
		T*							IndexToPointer_1D_Fixed		( const size_t* dimsindex ) const;
		T*							IndexToPointer_2D_Contiguous( const size_t* dimsindex ) const;
		T*							IndexToPointer_2D_Fixed		( const size_t* dimsindex ) const;
		void						ConfigureIndexingOperator	();
		bool						ConfigureMemoryLayout		( T* ptr );
		T*							LastAddressableLocation		() const;
		bool						InternalAttach				( size_t nrank, const size_t* dims, T* ptr, const size_t* strides, InxMemoryManager<T>* manager  );
		void						ConstructorInit				();
		bool						CheckSlicingDims			( const size_t indexlo[], const size_t indexhi[], size_t rank, size_t* numsubsetdims ) const;

	protected:
		bool						SetMemoryManager		( InxMemoryManager<T>* manager, T* baseptr );
		bool						SetTemporaryArrayParams	( InxMemoryManager<T>* manager, const RankSpecification* rank, T* baseptr );

	public:
									nxArrayLinear();
									nxArrayLinear(const nxArrayLinear<T>& other);
		virtual					   ~nxArrayLinear()									{ Detach(); }
virtual int							MandatoryRankValue		() { return 0;}																		// Returns the mandatory required rank.

		bool						Slice					( const size_t indexlo[], const size_t indexhi[], size_t nrank, nxArrayLinear<T>* column) const;			//!< extract a slice from this arrays
		void						SetTemporary			()									{ m_istemporary = true;}							//!< Flag this array as temporary. Reserved for internal usage only. Used when evaluationg operator overloads
		bool						SetSize					( size_t nrank, const size_t* dims, const size_t* strides = NULL);							//!< Set the size and layout of this array
		bool						Attach					( size_t nrank, const size_t* dims, T* ptr, InxMemoryManager<T>* manager = NULL, const size_t* strides = NULL );	//!< Attach this array to a memory manager defined area of memory
		void						Detach					( );																			//!< Detach this array from memory
		inline T&					At						( const size_t* dims )	      { return * ( (this->*m_indextopointer)(dims) ) ;}	//!< Return the element at the indexed location
		inline const T&				At						( const size_t* dims )	const { return * ( (this->*m_indextopointer)(dims) ) ;} //!< Return the element  at the indexed location (const version)
		inline T&					at						( const size_t* dims )	      { return * ( (this->*m_indextopointer)(dims) ) ;}	//!< Return the element at the indexed location
		inline const T&				at						( const size_t* dims )	const { return * ( (this->*m_indextopointer)(dims) ) ;} //!< Return the element  at the indexed location (const version)
		bool						At						( nxArrayLinear<size_t>& logicalindex, nxArrayLinear<T>* answer);				//!< Perform a logical indexing lookup.
		T							At						( size_t logicalindex ) const { return *(begin() + logicalindex);}

		void						SetTo					( const T& constant );														//!< Set all of the elements in the array to the specified value
		bool						ShallowCopy				( const nxArrayLinear<T>& other);											//!< Configure this array so it shares the same meory and layout as the otehr array
		bool						DeepCopy				( const nxArrayLinear<T>& other, bool copycontents = true);				//!< Copy the other array so we have a completely new instance in contiguous memory.
		void						SetReuseMemory			( bool allowreuse)			{ m_allowmemoryreuse = allowreuse;}			//!< Configure this array so it will reuse allocated memory when resizing rather than deallocating and allocating new memory. Net effect is that the array tends to grow.
		void						SetBoundsChecking		( bool checkbounds)			{ m_checkbounds     = checkbounds;}			//!< If TRUE then check bounds when using the #At methods
		bool						TrimSize				( size_t nx);																//!< Trims the array to the new size. Useful when shrinking arrays in the WHERE functions.
		bool						InputColumnMajorText	( const char * filename, size_t numcolumns = 0, size_t numdatalines = 0 );													//!< Read in 2-D compatible array from column major text from a file
		bool						InputColumnMajorText	( std::istream& instr,   size_t numcolumns = 0, size_t numdatalines = 0 );													//!< Read in 2-D compatible array from column major text from an input stream
		bool						IsSameDimensionalSize	( const nxArrayLinear<T>& other ) const { return m_rankspecs == other.m_rankspecs;}


		bool						IsTemporary				( ) const						{ return m_istemporary;}					//!< Return TRUE if this array is temporary. Used when evaluating operator overloads
		const RankSpecification*	ArrayRankSpecs			( ) const						{ return &m_rankspecs;}						//!< Return the layout description of this array
		InxMemoryManager<T>*		ArrayManager			( ) const						{ return m_manager;}						//!< Return the memory manger object used by this array, if any.
		T*							UnsafeArrayBasePtr		( ) const						{ return (T *)m_baseptr;}					//!< return a pointer to the base memory location of this array (usually the first location if offsets are zero)
		const T*					ArrayBasePtr			( ) const						{ return (T *)m_baseptr;}					//!< return a pointer to the base memory location of this array (usually the first location if offsets are zero)
		bool						IsEmpty					( ) const						{ return m_baseptr == NULL;}				//!< Return true if the array is empty, no dimensions and no memory
		nxArrayIter<T>				begin					( ) const	{ return nxArrayIter<T>( (nxBYTE *)m_baseptr, &m_rankspecs);}	//!< Return an STL random access, compliant iterator for the beginning of the array
		nxArrayIter<T>				end						( ) const	{ return nxArrayIter<T>( (nxBYTE *)m_endptr,  &m_rankspecs);}	//!< Return an STL random access, compliant iterator for the end of the array. Note it may not be possible to dereference the end iterator
		size_t						size					( ) const	{ return m_rankspecs.N_Elements();}								//!< Return the number of elements in this array
		void						erase					( );																		//!< empty the array and free any allocated resources
		T							back                    ( ) const   { return *LastAddressableLocation();}
		T&							back                    ( )         { return *LastAddressableLocation();}
		T							front					( ) const   { return *m_baseptr;}											//!< Return the first value
		T&							front					( )         { return *m_baseptr;}

		nxArrayLinear<T>&			operator =				( const nxArrayLinear<T>& other )	{ DeepCopy(other); return *this; }
		nxArrayLinear<T>&			operator =				( const T& constant)				{ SetTo(constant); return *this; }

};

template <class T> std::ostream& operator << ( std::ostream& outputstream, nxArrayLinear<T>& userdata );
template <class T> std::istream& operator >> ( std::istream& inputstream, nxArrayLinear<T>& userdata );




/*<--------------------------------------------------------------------------
 *'					nxArrayTemporary                                2002-5-27*/
/**	\ingroup nxLinearArray_Internals
 *	A array that is specifically designed to allow efficient movement of
 *  large temporary arrays through the code when executing binary and
 *	assignment operators.  It always allocates memory dynamically.
 *	NxArrayTemporary is allowed exclusive control of the m_istemporary
 *	variable in nxArrayLinear.
 *
 *	Note:  nxArrayTemporary will never attach to memory areas.
 *		   nxArrayTemporary will only share to with other nxArrayTemporary objects.
 *		   nxArrayTemporary will only allocate dynamic memory managers.
 *		   nxArrayTemporary will only allocate memory with a stride of 1.
 */
 /*------------------------------------------------------------------------*/



template <class T >
class nxArrayTemporary : public nxArrayLinear<T>
{
	public:
							nxArrayTemporary(){ this->SetTemporary(); }
							nxArrayTemporary( const nxArrayTemporary<T>& other);
						   ~nxArrayTemporary(){};
		bool				PrepareBinaryOperator( const nxArrayLinear<T>& a, const nxArrayLinear<T>& b);
};



/*<--------------------------------------------------------------------------
 *'					class nx1dArray									2002-9-10 */
/** \ingroup nxLinearArray
 *	Class for managing 1D arrays.  Provides support for arrays which are
 *	a fixed stride in memory.
 *
 */
/*------------------------------------------------------------------------*/

template< class T >
class nx1dArray : public nxArrayLinear<T>
{

	protected:
	virtual int				MandatoryRankValue()							{ return 1;}			//!< Returns the mandatory required rank.

	public:
							nx1dArray	()									{ }						//!< Default constructor creates an empty array
							nx1dArray	(size_t nx)							{ SetSize(nx);}			//!< Constructor creates an array of the specifed number of elements
							nx1dArray	(size_t nx, T*ptr )					{ Attach( nx, ptr);}	//!< Constructor attaches this array to the specified memory location with the specified number of elements
							nx1dArray	( const nx1dArray& other ) : nxArrayLinear<T>(other) { }	//!< Copy constructor.  Attaches this array to the meory used by the other array.
		virtual			   ~nx1dArray	(){}

		nx1dArray<T>&		Indgen		( size_t nx);												//!< Create the array of size \b nx elements filled with ascending integeres 0 to \b nx-1
		bool				Slice		( size_t lo, size_t hi, nxArrayLinear<T>* column) const;		//!< Extract a subset from the 1d array use an array passed in by user
		nx1dArray<T>		Slice		( size_t lo, size_t hi) const;								//!< Extract a subset from the 1d array. Slightly slower than other Slice methods but maybe more convenient

		bool				SetSize		( size_t nx  )						{ return nxArrayLinear< T>::SetSize( 1, &nx );}			//!< Set the size of this array to \b nx elements
		bool				Attach		( size_t nx, T* ptr )				{ return nxArrayLinear< T>::Attach ( 1, &nx, ptr );}	//!< Attach this array to \b nx contiguous elements  at memory location \b ptr
		inline T&			At			( size_t nx )						{ return nxArrayLinear< T>::At(&nx);}					//!< Get the modifiable value of the \b nx'th element
        inline const T&		At			( size_t nx ) const					{ return nxArrayLinear< T>::At(&nx);}					//!< Get the non-modifiable value of the \b nx'th element
		inline T&			at			( size_t nx )						{ return nxArrayLinear< T>::At(&nx);}					//!< Get the modifiable value of the \b nx'th element
        inline const T&		at			( size_t nx ) const					{ return nxArrayLinear< T>::At(&nx);}					//!< Get the non-modifiable value of the \b nx'th element
		inline T&			operator [] ( size_t nx )						{ return nxArrayLinear< T>::At(&nx);}					//!< Get the modifiable value of the \b nx'th element. Same as #At.
        inline const T&		operator [] ( size_t nx ) const					{ return nxArrayLinear< T>::At(&nx);}					//!< Get the non-modifiable value of the \b nx'th element. Same as #At.
		nx1dArray&			operator =  (const nx1dArray<T>&     other )	{ this->DeepCopy(other); return *this; }						//!< Assign the other array to this one using a deep copy.
		nx1dArray&			operator =  (const nxArrayLinear<T>& other )	{ this->DeepCopy(other); return *this; }						//!< Assign the other array to this one using a deep copy.
		nx1dArray&			operator =  (const T& constant )				{ this->SetTo(constant); return *this; }						//!< Assign the scalar constant to each element.
		bool				WriteToTextFile	( const char* filename, bool columnmajor, const char* formatstr );
		std::vector<T>		STLVector	( ) const;

};

template <class T> std::ostream& operator << ( std::ostream& outputstream, nx1dArray<T>& userdata ) { return operator<<(outputstream, dynamic_cast<nxArrayLinear<T>&>(userdata) );}
template <class T> std::istream& operator >> ( std::istream& inputstream,  nx1dArray<T>& userdata ) { return operator>>(inputstream,  dynamic_cast<nxArrayLinear<T>&>(userdata) );}


/*<--------------------------------------------------------------------------
 *'					class nx2dArray									2002-9-10*/
/** \ingroup nxLinearArray
 *	Class for managing 2D arrays. The arrays are typically 2-D X/Y arrays where
 *	the X is the horizontal and Y is the vertical.  In this case the user would
 *	index the array with something like \c array->At(x,y). The class uses 0
 *	based indexing and the X axis changes most rapidly in memory.
 */
/*------------------------------------------------------------------------*/

template< class T >
class nx2dArray : public nxArrayLinear<T>
{

	protected:
		virtual int		MandatoryRankValue()								{ return 2;}			//!< Returns the mandatory required rank.

	public:
							nx2dArray		(){}															//!< Default constructor creates an empty array
							nx2dArray		(const nx2dArray& other )  : nxArrayLinear<T>(other)  {}		//!< Makes a shallow copy of the other array. Allows for fairly efficient passing of arrays as value arguments although the user may be surprised that the array is changed.
							nx2dArray		(size_t nx, size_t ny)				{ SetSize(nx,ny);}			//!< Constructor creates an array of the specifed number of elements
							nx2dArray		(size_t nx, size_t ny, T*ptr )		{ Attach( nx,ny, ptr);}		//!< Constructor attaches this array to the specified memory location with the specified number of elements
		virtual			   ~nx2dArray		()									{ }
		bool				Slice			( size_t lo0, size_t hi0, size_t lo1, size_t hi1, nxArrayLinear<T>* column) const;
		nx2dArray<T>		Slice			( size_t lo0, size_t hi0, size_t lo1, size_t hi1 ) const;
		nx1dArray<T>& 		XSlice			( size_t y,  nx1dArray<T>* row ) const		{ Slice( NXARRAY_STARSELECT, NXARRAY_STARSELECT,   y, y,   row); return *row;}
		nx1dArray<T>& 		YSlice			( size_t x,  nx1dArray<T>* row ) const		{ Slice(  x,  x,  NXARRAY_STARSELECT,NXARRAY_STARSELECT,   row); return *row;}
		std::vector<T>& 	XSlice			( size_t y,  std::vector<T>*  arow ) const	{ nx1dArray<T> row; Slice( NXARRAY_STARSELECT, NXARRAY_STARSELECT,   y, y,   &row); arow->assign( row.begin(), row.end()); return *arow;}
		std::vector<T>& 	YSlice			( size_t x,  std::vector<T>*  arow ) const	{ nx1dArray<T> row; Slice(  x,  x,  NXARRAY_STARSELECT, NXARRAY_STARSELECT,   &row); arow->assign( row.begin(), row.end()); return *arow;}
		bool				SetSize			( size_t nx, size_t ny  )			{ size_t dims[2]={nx,ny}; return nxArrayLinear< T>::SetSize( 2, dims );}		//!< Set the size of this array to \b nx by \b ny elements
		bool				Attach			( size_t nx, size_t ny, T* ptr )	{ size_t dims[2]={nx,ny}; return nxArrayLinear< T>::Attach ( 2, dims, ptr);}	//!< Attach this array to \b nx by \b ny contiguous elements  at memory location \b ptr
		inline T&			At				( size_t nx, size_t ny )			{ size_t dims[2]={nx,ny}; return nxArrayLinear< T>::At     ( dims);}			//!< Get the modifiable value of the \b (nx,ny)'th element
		inline const T&		At				( size_t nx, size_t ny ) const		{ size_t dims[2]={nx,ny}; return nxArrayLinear< T>::At     ( dims);}			//!< Get the non-modifiable value of the \b (nx,ny)'th element
		inline T&			at				( size_t nx, size_t ny )			{ size_t dims[2]={nx,ny}; return nxArrayLinear< T>::At     ( dims);}			//!< Get the modifiable value of the \b (nx,ny)'th element
		inline const T&		at				( size_t nx, size_t ny ) const		{ size_t dims[2]={nx,ny}; return nxArrayLinear< T>::At     ( dims);}			//!< Get the non-modifiable value of the \b (nx,ny)'th element
		size_t				XSize			() const							{ return this->ArrayRankSpecs()->Dims()[0];}										//!< Return the size of the first (X) dimension
		size_t				YSize			() const							{ return this->ArrayRankSpecs()->Dims()[1];}										//!< Return the size of the second (Y) dimension
//		nx1dArray<T>		operator[]		(size_t row)						{ nx1dArray<T> answer; XSlice( (int)row, &answer); answer.SetTemporary(); return answer;}				//!< Returns the specified row. Same functionality as #GetRow but slower as an array copy is involved.
		nx2dArray&			operator =      (const nx2dArray<T>& other )		{ this->DeepCopy(other); return *this; }											//!< Assign the other array to this one using a deep copy.
		nx2dArray&			operator =      (const nxArrayLinear<T>& other )	{ this->DeepCopy(other); return *this; }											//!< Assign the other array to this one using a deep copy.
		nx2dArray&			operator =		(const T& constant )				{ this->SetTo(constant); return *this; }											//!< Assigns the constant value to every element in this array
		bool				WriteToTextFile	( const char* filename, bool columnmajor, const char* formatstr );
};

/*---------------------------------------------------------------------------
 *'					class nx3dArray									2002-9-10*/
/** \ingroup nxLinearArray
 *	Class for managing 3D arrays. The arrays are typically X/Y/Z arrays. In
 *	this case the user would index the array with something like \c array->At(x,y,z).
 *	The class uses 0 based indexing. The X axis changes most rapidly in memory and
 *	the Z axis changes least frequently.
 *
 */
/*-------------------------------------------------------------------------*/

template< class T >
class nx3dArray : public nxArrayLinear<T>
{

	protected:
	virtual int			MandatoryRankValue()							{ return 3;}			//!< Returns the mandatory required rank.

	public:
								nx3dArray	(){}																		//!< Default constructor creates an empty array
								nx3dArray	(size_t nx, size_t ny, size_t nz)					{ SetSize(nx,ny,nz);}			//!< Constructor creates an array of the specifed number of elements
								nx3dArray	(size_t nx, size_t ny, size_t nz, T*ptr )			{ Attach (nx,ny,nz, ptr);}		//!< Constructor attaches this array to the specified memory location with the specified number of elements
								nx3dArray	( const nx3dArray& other )  : nxArrayLinear<T>(other)  {}					//!< Makes a shallow copy of the other array. Makefs for fairly quick passing of arrays as arguments although the user may be suprised that the mmeory values can be changed.
		virtual				   ~nx3dArray		(){}
		bool					Slice( size_t lo0, size_t hi0, size_t lo1, size_t hi1, size_t lo2, size_t hi2, nxArrayLinear<T>* column) const;
		nx3dArray<T>			Slice( size_t lo0, size_t hi0, size_t lo1, size_t hi1, size_t lo2, size_t hi2 ) const;

		bool					SetSize		( size_t nx, size_t ny, size_t nz  )		{ size_t dims[3]={nx,ny,nz}; return nxArrayLinear< T>::SetSize( 3, dims );}			//!< Set the size of the array
		bool					Attach		( size_t nx, size_t ny, size_t nz, T* ptr )	{ size_t dims[3]={nx,ny,nz}; return nxArrayLinear< T>::Attach ( 3, dims, ptr);}	//!< Attach this array to \b nx by \b ny by \b nz contiguous elements at memory location \b ptr
		T&						At			( size_t nx, size_t ny, size_t nz )			{ size_t dims[3]={nx,ny,nz}; return nxArrayLinear< T>::At     ( dims);}			//!< Get the modifiable value of the \b (nx,ny,nz)'th element
		const T&				At			( size_t nx, size_t ny, size_t nz ) const	{ size_t dims[3]={nx,ny,nz}; return nxArrayLinear< T>::At     ( dims);}			//!< Get the non-modifiable value of the \b (nx,ny)'th element
		T&						at			( size_t nx, size_t ny, size_t nz )			{ size_t dims[3]={nx,ny,nz}; return nxArrayLinear< T>::At     ( dims);}			//!< Get the modifiable value of the \b (nx,ny,nz)'th element
		const T&				at			( size_t nx, size_t ny, size_t nz ) const	{ size_t dims[3]={nx,ny,nz}; return nxArrayLinear< T>::At     ( dims);}			//!< Get the non-modifiable value of the \b (nx,ny)'th element
		nx1dArray<T>& 			XSlice		( int y,int z, nx1dArray<T>* row   ) const		{ return Slice( NXARRAY_STARSELECT, NXARRAY_STARSELECT,   y, y,   z, z, row);}
		nx1dArray<T>& 			YSlice		( int x,int z, nx1dArray<T>* row   ) const		{ return Slice(  x,  x,  NXARRAY_STARSELECT, NXARRAY_STARSELECT,   z, z, row);}
		nx1dArray<T>& 			ZSlice		( int x,int y, nx1dArray<T>* row   ) const		{ return Slice(  x,  x,   y, y,  NXARRAY_STARSELECT, NXARRAY_STARSELECT, row);}
		size_t					XSize		() const								{ return this->ArrayRankSpecs()->Dims()[0];}		//!< Return the size of the first (X) dimension
		size_t					YSize		() const								{ return this->ArrayRankSpecs()->Dims()[1];}		//!< Return the size of the second (Y) dimension
		size_t					ZSize		() const								{ return this->ArrayRankSpecs()->Dims()[2];}		//!< Return the size of the third (Z) dimension
		nx3dArray&				operator =  (const nx3dArray<T>& other )			{ this->DeepCopy(other); return *this; }			//!< Assign the other array to this one using a deep copy.
		nx3dArray&				operator =  (const nxArrayLinear<T>& other )		{ this->DeepCopy(other); return *this; }			//!< Assign the other array to this one using a deep copy.
};

#include "nxarrayiter.hpp"
#include "nxmemorymanager.hpp"
#include "nxarraylinear.hpp"
#include "nxarraytemporary.hpp"



/*---------------------------------------------------------------------------
 *	Here are a set of macros used to quickly code up the differnt operators
 *	for the arry classes
 *-------------------------------------------------------------------------*/

// --- this MACRO is for binary mathematical operators like a+b where a and b are both arrays

#define nxArrayLinearBinaryOperator( assignment_operator)							\
{																					\
	nxArrayTemporary<A>		temp;													\
	bool					ok;														\
																					\
	ok = temp.PrepareBinaryOperator(a,b);											\
	if (!ok)  nxLog::Record(NXLOG_WARNING, "nxArrayLinearBinaryOperator, Error preparing arrays are not the same size, operation not performed"); \
	else																			\
	{																				\
		nxArrayIter<A>	titer	   = temp.begin();									\
		nxArrayIter<A>	aiter      = a.begin();										\
		nxArrayIter<A>	aend       = a.end();										\
		nxArrayIter<A>	biter      = b.begin();										\
																					\
		while ( aiter != aend ){*titer = *aiter assignment_operator *biter; ++aiter; ++biter; ++titer;} \
	}																				\
	return temp;																	\
}

// --- this MACRO is for binary mathematical operators like a+b where a is an array and b is a scalar

#define nxArrayLinearBinaryScalar( assignment_operator)								\
{																					\
	nxArrayTemporary<A>		temp;													\
	bool					ok;														\
																					\
	temp.DeepCopy(a, false);														\
	ok = (temp.size() == a.size());													\
	if (!ok)  nxLog::Record(NXLOG_WARNING, "nxArrayLinearBinaryScalar, temp array is not correct size, operation not performed"); \
	else																			\
	{																				\
		nxArrayIter<A>	titer	   = temp.begin();									\
		nxArrayIter<A>	aiter      = a.begin();										\
		nxArrayIter<A>	aend       = a.end();										\
																					\
		while ( aiter != aend ){ *titer = *aiter assignment_operator b; ++aiter; ++titer;} \
	}																				\
	return temp;																	\
}

// --- this MACRO is for binary mathematical operators like b+a where a is an array and b is a scalar

#define nxArrayLinearBinaryReverseScalar( assignment_operator)						\
{																					\
	nxArrayTemporary<A>		temp;													\
	bool					ok;														\
																					\
	temp.DeepCopy(a, false);														\
	ok = (temp.size() == a.size());													\
	if (!ok)  nxLog::Record(NXLOG_WARNING, "nxArrayLinearBinaryReverseScalar, temp array is not correct size, operation not performed"); \
	else																			\
	{																				\
		nxArrayIter<A>	titer	   = temp.begin();									\
		nxArrayIter<A>	aiter      = a.begin();										\
		nxArrayIter<A>	aend       = a.end();										\
																					\
		while ( aiter != aend ){ *titer = b assignment_operator *aiter ; ++aiter; ++titer;} \
	}																				\
	return temp;																	\
}


// --- this MACRO is for comparison operators like max (a, b) where a is an array and b is a scalar

#define nxArrayLinearBinaryComparator( comparefunc )                                \
{																					\
	nxArrayTemporary<A>		temp;													\
	bool					ok;														\
																					\
	temp.DeepCopy(a, false);														\
	ok = (temp.size() == a.size());													\
	if (!ok)  nxLog::Record(NXLOG_WARNING, "nxArrayLinearBinaryScalar, temp array is not correct size, operation not performed"); \
	else																			\
	{																				\
		nxArrayIter<A>	titer	   = temp.begin();									\
		nxArrayIter<A>	aiter      = a.begin();										\
		nxArrayIter<A>	aend       = a.end();										\
																					\
		while ( aiter != aend ){ *titer = comparefunc(*aiter,b);++titer; ++aiter;}	\
	}																				\
	return temp;																	\
}


// --- this MACRO is for unary mathematical assignment operators like a += b where a and b are both arrays

#define nxArrayLinearUnaryAssignment( assignment_operator)								\
{																						\
	bool					ok;															\
																						\
	ok = (a.size() == b.size());														\
	if (!ok)  nxLog::Record(NXLOG_WARNING, "nxArrayLinearUnaryAssignment, arrays are not the same size, operation not performed"); \
	else																				\
	{																					\
		nxArrayIter<A>	aiter      = a.begin();											\
		nxArrayIter<A>	aend       = a.end();											\
		nxArrayIter<A>	biter      = b.begin();											\
																						\
		while ( aiter != aend ){ *aiter assignment_operator *biter; ++aiter; ++biter;}	\
	}																					\
	return a;																			\
}

// --- this MACRO is for unary mathematical assignment operators like a += b where a is array and b is a scalar

#define nxArrayLinearScalarUnaryAssignment( assignment_operator)					\
{																					\
																					\
	nxArrayIter<A>	aiter      = a.begin();											\
	nxArrayIter<A>	aend       = a.end();											\
																					\
	while ( aiter != aend ) {*aiter assignment_operator b; ++aiter;	}				\
	return a;																		\
}

// --- this MACRO is for mathematical functions line cos(a) where a is an array userfunction is a simple function


#define nxArrayLinearUnaryFunction( userfunction )								\
{																				\
	nxArrayTemporary<A>		temp;												\
	bool					ok;													\
																				\
	temp.DeepCopy(a, false);													\
	ok = (temp.size() == a.size());												\
	if (!ok)  nxLog::Record(NXLOG_WARNING, "nxArrayLinearUnaryFunction, temp array is not correct size, operation not performed");\
	else																		\
	{																			\
		nxArrayIter<A>	titer	   = temp.begin();								\
		nxArrayIter<A>	aiter      = a.begin();									\
		nxArrayIter<A>	aend       = a.end();									\
		while ( aiter != aend ) {*titer = userfunction(*aiter); ++titer; ++aiter;}		\
	}																			\
	return temp;																\
}

//---------------------------------------------------------------------------
// Here we define a whole range of different operators for our arrays
//---------------------------------------------------------------------------

template <class A> nxArrayLinear<A>	operator +  ( const nxArrayLinear<A>& a, const nxArrayLinear<A>& b ){ nxArrayLinearBinaryOperator( + ) }
template <class A> nxArrayLinear<A>	operator -  ( const nxArrayLinear<A>& a, const nxArrayLinear<A>& b ){ nxArrayLinearBinaryOperator( - ) }
template <class A> nxArrayLinear<A>	operator *  ( const nxArrayLinear<A>& a, const nxArrayLinear<A>& b ){ nxArrayLinearBinaryOperator( * ) }
template <class A> nxArrayLinear<A>	operator /  ( const nxArrayLinear<A>& a, const nxArrayLinear<A>& b ){ nxArrayLinearBinaryOperator( / ) }
template <class A> nxArrayLinear<A>	operator |  ( const nxArrayLinear<A>& a, const nxArrayLinear<A>& b ){ nxArrayLinearBinaryOperator( | ) }
template <class A> nxArrayLinear<A>	operator &  ( const nxArrayLinear<A>& a, const nxArrayLinear<A>& b ){ nxArrayLinearBinaryOperator( & ) }
template <class A> nxArrayLinear<A>	operator ^  ( const nxArrayLinear<A>& a, const nxArrayLinear<A>& b ){ nxArrayLinearBinaryOperator( ^ ) }
template <class A> nxArrayLinear<A>	operator %  ( const nxArrayLinear<A>& a, const nxArrayLinear<A>& b ){ nxArrayLinearBinaryOperator( % ) }

template <class A> nxArrayLinear<A>	operator +  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryScalar( + ) }
template <class A> nxArrayLinear<A>	operator -  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryScalar( - ) }
template <class A> nxArrayLinear<A>	operator *  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryScalar( * ) }
template <class A> nxArrayLinear<A>	operator /  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryScalar( / ) }
template <class A> nxArrayLinear<A>	operator |  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryScalar( | ) }
template <class A> nxArrayLinear<A>	operator &  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryScalar( & ) }
template <class A> nxArrayLinear<A>	operator ^  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryScalar( ^ ) }
template <class A> nxArrayLinear<A>	operator %  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryScalar( % ) }
template <class A> nxArrayLinear<A>	operator >  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryComparator( nxmax<A> ) }
template <class A> nxArrayLinear<A>	operator <  ( const nxArrayLinear<A>& a, A b )	{ nxArrayLinearBinaryComparator( nxmin<A> ) }


template <class A> nxArrayLinear<A>	operator +  ( A b, const nxArrayLinear<A>& a )	{ nxArrayLinearBinaryReverseScalar( + ) }
template <class A> nxArrayLinear<A>	operator -  ( A b, const nxArrayLinear<A>& a )	{ nxArrayLinearBinaryReverseScalar( - ) }
template <class A> nxArrayLinear<A>	operator *  ( A b, const nxArrayLinear<A>& a )	{ nxArrayLinearBinaryReverseScalar( * ) }
template <class A> nxArrayLinear<A>	operator /  ( A b, const nxArrayLinear<A>& a )	{ nxArrayLinearBinaryReverseScalar( / ) }
template <class A> nxArrayLinear<A>	operator |  ( A b, const nxArrayLinear<A>& a )	{ nxArrayLinearBinaryReverseScalar( | ) }
template <class A> nxArrayLinear<A>	operator &  ( A b, const nxArrayLinear<A>& a )	{ nxArrayLinearBinaryReverseScalar( & ) }
template <class A> nxArrayLinear<A>	operator ^  ( A b, const nxArrayLinear<A>& a )	{ nxArrayLinearBinaryReverseScalar( ^ ) }
template <class A> nxArrayLinear<A>	operator %  ( A b, const nxArrayLinear<A>& a )	{ nxArrayLinearBinaryReverseScalar( % ) }


template <class A> nxArrayLinear<A>&	operator += ( nxArrayLinear<A>& a, const nxArrayLinear<A>& b )		{ nxArrayLinearUnaryAssignment  ( += ) }
template <class A> nxArrayLinear<A>&	operator -= ( nxArrayLinear<A>& a, const nxArrayLinear<A>& b )		{ nxArrayLinearUnaryAssignment  ( -= ) }
template <class A> nxArrayLinear<A>&	operator *= ( nxArrayLinear<A>& a, const nxArrayLinear<A>& b )		{ nxArrayLinearUnaryAssignment  ( *= ) }
template <class A> nxArrayLinear<A>&	operator /= ( nxArrayLinear<A>& a, const nxArrayLinear<A>& b )		{ nxArrayLinearUnaryAssignment  ( /= ) }
template <class A> nxArrayLinear<A>&	operator |= ( nxArrayLinear<A>& a, const nxArrayLinear<A>& b )		{ nxArrayLinearUnaryAssignment  ( |= ) }
template <class A> nxArrayLinear<A>&	operator &= ( nxArrayLinear<A>& a, const nxArrayLinear<A>& b )		{ nxArrayLinearUnaryAssignment  ( &= ) }
template <class A> nxArrayLinear<A>&	operator ^= ( nxArrayLinear<A>& a, const nxArrayLinear<A>& b )		{ nxArrayLinearUnaryAssignment  ( ^= ) }

template <class A> nxArrayLinear<A>&	operator += ( nxArrayLinear<A>& a, const A b )	{ nxArrayLinearScalarUnaryAssignment  ( += ) }
template <class A> nxArrayLinear<A>&	operator -= ( nxArrayLinear<A>& a, const A b )	{ nxArrayLinearScalarUnaryAssignment  ( -= ) }
template <class A> nxArrayLinear<A>&	operator *= ( nxArrayLinear<A>& a, const A b )	{ nxArrayLinearScalarUnaryAssignment  ( *= ) }
template <class A> nxArrayLinear<A>&	operator /= ( nxArrayLinear<A>& a, const A b )	{ nxArrayLinearScalarUnaryAssignment  ( /= ) }
template <class A> nxArrayLinear<A>&	operator |= ( nxArrayLinear<A>& a, const A b )	{ nxArrayLinearScalarUnaryAssignment  ( |= ) }
template <class A> nxArrayLinear<A>&	operator &= ( nxArrayLinear<A>& a, const A b )	{ nxArrayLinearScalarUnaryAssignment  ( &= ) }
template <class A> nxArrayLinear<A>&	operator ^= ( nxArrayLinear<A>& a, const A b )	{ nxArrayLinearScalarUnaryAssignment  ( ^= ) }


#include "nxarray_algorithms.hpp"
#endif
