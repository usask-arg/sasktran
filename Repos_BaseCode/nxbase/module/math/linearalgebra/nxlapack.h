//#include "lapackdefs.h"

class LapackFactorLU;
class LapackFactorCholesky;
class LapackFactorSymEigen;
class LapackSymmetricMatrix;
class LapackMatrix;

/*-----------------------------------------------------------------------------
 *					class LapackBLAS								2004-11-23*/
/** \ingroup Math_Matrix
 *	A base class that wraps the BLAS (Basic Linear Algebra Subroutines) functions.
 *	The BLAS is available from the INtel Math kernel Library and is also
 *	available as Fortran source code.  This intent is to grow class so it has much
 *	better coverage of the entireBLAS.
**/
/*---------------------------------------------------------------------------*/

class LapackBLAS
{
	public:
	
	static char				UpcaseChar( char c );
	static nxBOOL			Gemm(												//!<  Matrix Matrix multiplications
										double					dalpha,			// Typically 1.0
										const LapackMatrix&		A,				// Matrix A
										char					transa,			// ' ' use A  'T' use Transpose(A)
										const LapackMatrix&		B,				// Matrix B
										char					transb,			// ' ' use B  'T' use Transpose(B)
										LapackMatrix*			C,				// Answer matrix
										double					dbeta = 0.0		// typically 0.0
									 );

	static nxBOOL			Dtrsm(												//!< Solve AX=B where A is triangular 
										char							side,
										const LapackSymmetricMatrix&	A,
										char							transa,
										const LapackMatrix&				B,
										LapackMatrix*					X,
										double							dalpha
									  );
};

/*-----------------------------------------------------------------------------
 *					class LapackMatrix								2004-11-23*/
/** \ingroup Math_Matrix
 *	A class specifically designed for matrix mathematics. The class builds upon the
 *	Linear Algebra Package (LAPACK) which is available in the Intel Math kernel Library
 *	and is also available as Fortran Source code.  This class is meant to represent 
 *	matrices. Although this class derives from nx2dArray<double> it does have some very
 *	specific differences so we can make the code \i look like the formulae typed in text books:
 *		- indexing is 1 based rather than the normal 'C' 0 based
 *		- The matrices are stored internally in  column major so it is consistent
 *		   with normal FORTRAN definitions
 *		-  the code makes unashamed, extensive use of the LAPACK fortran libraries.
 *
 * \par Operators
 *	The LapackMatrix class overloads most of the sensible mathematical operators, +,-,/,*
 *	as well as the function operator (i,j). This allows one to write A(i,j) as a shorthand
 *	notation for A.At(i,j) with very little penalty.
 **/
/*---------------------------------------------------------------------------*/

class LapackMatrix: public nx2dArray<double>
{

	public:
		LapackMatrix			LUInverse		(); 											//!< LAPACK dgetrf followed by dgetri
		LapackSymmetricMatrix	CholeskyInverse	(); 											//!< LAPACK dgetrf followed by dgetri
		LapackMatrix			SolveSymmetricEigen( nx1dArray<double>* eigenvalues );			//!< Get Eigenvectors and eigenvalues of a symmetric matrix

	public:
							LapackMatrix	(){};																			//!< create empty matrix
							LapackMatrix	( size_t numrows, size_t numcols) : nx2dArray<double>( numrows, numcols) {};			//!< Create the array with the specified size
							LapackMatrix	( size_t numrows, size_t numcols, double* storage) : nx2dArray<double>( numrows, numcols, storage) {}	//!< Create the array with the specified size and storage area 
							LapackMatrix	( const LapackMatrix& copy):nx2dArray<double>(copy){}							//!< Copy the other matrix to this instance
						   ~LapackMatrix	(){};																			//!< default destructor

		nxBOOL				Create					( size_t numrows, size_t numcols, double initialvalue = 0.0 );			//!< Create the specified matrix with the specified intial value
		bool				CreateFromVectorVector	( const std::vector< std::vector<double> >&	x);
		bool				CreateFromRowVector		( const std::vector< double>&				x);
		bool				CreateFromColumnVector	( const std::vector< double>&				x);
		nxBOOL				CreateUnity				( size_t numrows );																//!< Create a square unit matrix, 1's on the diagonal
		void				Dump			( const char* message );											//!< Dump the matrix to stdout
		nxBOOL				Initialize		( size_t numrows, size_t numcols, ... );
		void				Clear			()									{ SetTo(0);}									//!< Set all elements of the matrix to 0
		size_t				NumRows			() const							{ return (int)XSize();}								//!< Return the vertical size of the matrix
		size_t				NumColumns		() const							{ return (int)YSize();}								//!< Return the horizontal size of the matrix.
		double&				At				( size_t row, size_t col)					{ return nx2dArray<double>::At(--row, --col);}	//!< Get the element at the (row,column), using 1's based indexing (i.e. like indexing used in math formulae)
		double				At				( size_t row, size_t col)	const			{ return nx2dArray<double>::At(--row, --col);}	//!< Get the element at the (row,column), using 1's based indexing (i.e. like indexing used in math formulae)
		LapackMatrix&		Slice			( size_t lo1, size_t hi1, size_t lo2, size_t hi2, LapackMatrix* slice) ;
		LapackMatrix&		Row				( size_t row, LapackMatrix* slice)		 { return Slice(row,row,  NXARRAY_STARSELECT, NXARRAY_STARSELECT,   slice);}
		LapackMatrix&		Col				( size_t col, LapackMatrix* slice)		 { return Slice( NXARRAY_STARSELECT, NXARRAY_STARSELECT,   col,col, slice);}
		LapackMatrix		Transpose		();
		void				CopyData        ( const nxArrayLinear<double>& copy) { DeepCopy(copy);}								//!< Make a fresh copy of the other matrix into this one
		LapackMatrix&		operator =		( const nxArrayLinear<double>& copy) { this->DeepCopy(copy); return *this;}				//!< Make a fresh copy of the otehr matrix into this one.
		double&				operator ()		( size_t row, size_t col)			{ return At(row,col);}							//!< A(1,3) is the same as A.At(1,3)
		LapackMatrix&		operator  =     ( double scalar )					{  nx2dArray<double>::operator= (scalar);  return *this;}
		LapackMatrix&		operator +=		( double scalar )					{ (nxArrayLinear<double> &)(*this) += scalar; return *this;}
		LapackMatrix&		operator -=		( double scalar )					{ (nxArrayLinear<double> &)(*this) -= scalar; return *this;}
		LapackMatrix&		operator /=		( double scalar )					{ (nxArrayLinear<double> &)(*this) /= scalar; return *this;}
		LapackMatrix&		operator *=		( double scalar )					{ (nxArrayLinear<double> &)(*this) *= scalar; return *this;}
		LapackMatrix&		operator -=		( const LapackMatrix& B)			{ (nxArrayLinear<double> &)(*this) -= B;      return *this;}
		LapackMatrix&		operator +=		( const LapackMatrix& B)			{ (nxArrayLinear<double> &)(*this) += B;	  return *this;}
							operator double*() const							{ return UnsafeArrayBasePtr();}

	public:									// ---- interfaces to LAPACK routines
		bool				Dgels			( const nxArrayLinear<double>& user_y, nxArrayLinear<double>* x );

};

LapackMatrix operator +( const LapackMatrix& A, double scalar );//					{ LapackMatrix other = *this; other.m_array += scalar; return other; }
LapackMatrix operator -( const LapackMatrix& A, double scalar );//					{ LapackMatrix other = *this; other.m_array -= scalar; return other; }
LapackMatrix operator /( const LapackMatrix& A, double scalar );//					{ LapackMatrix other = *this; other.m_array /= scalar; return other; }
LapackMatrix operator *( const LapackMatrix& A, double scalar );//					{ LapackMatrix other = *this; other.m_array *= scalar; return other; }
//LapackMatrix operator *( double scalar,         const  LapackMatrix& other );		//needs to be written
LapackMatrix operator *( const LapackMatrix& A, const  LapackMatrix& B);//			{ LapackMatrix C; LapackBLAS::GEMM(  1.0, *this, ' ',B,' ',&C,0.0); return C;}


/*-----------------------------------------------------------------------------
 *					class LapackFactorLU		2004-11-23*/
/** \ingroup Math_Matrix
 *	A specialty class specifically designed to determine and hold the matrix from
 *	LU factorization of a regular matrix.  This is also valuable when solving
 *	matrix equations of the form \b AX=B
 *
 *  LU factorization of a general M-by-N matrix A
 *  using partial pivoting with row interchanges.
 *  The factorization has the form
 *  \code A = P * L * U \endcode
 *  where P is a permutation matrix, L is lower triang  ular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
**/
/*---------------------------------------------------------------------------*/

class LapackFactorLU 
{
	private:
		LapackMatrix		m_A;
		nx1dArray<int>		m_pivots;
		nxBOOL				m_isfactored;

	private:
		nxBOOL				AllocatePivots	( int numpivots ) 			{ m_pivots.SetSize( (size_t)numpivots); return ( (int)m_pivots.size() == numpivots);}	// Allocate	pivots
		int*				PivotArray		()							{ return m_pivots.UnsafeArrayBasePtr();}

	public:
							LapackFactorLU() {m_isfactored = nxFALSE;}


		nxBOOL				GetInverse ( LapackMatrix *  inverse); 			//!< Get the inverse of the original matrix that is factorized in this instance 
		nxBOOL				LUFactorize( const LapackMatrix& A );			//!< LU Factorize matrix A and store results internally.
		nxBOOL				Solve( const LapackMatrix& B, LapackMatrix* X);	//!< Solve the matrix equation \b AX=B where this instance holds te LU factorization of \b A
		nxBOOL				Solve( const LapackMatrix& A, const LapackMatrix& B, LapackMatrix* X);	//!< Solve the matrix equation \b AX=B where this instance holds te LU factorization of \b A
};


/*-----------------------------------------------------------------------------
 *					class LapackSymmetricMatrix						2004-11-23*/
/** \ingroup Math_Matrix
 *	A class used to represent symmetric matrices used in the LAPACK library.
 *	These matrices have a different memory layout to regular matrices. This
 *	class is used by other parts of the library
**/
/*---------------------------------------------------------------------------*/

class LapackSymmetricMatrix : public LapackMatrix
{
	protected:
		char				m_uplo;

	private:
		void					Inituplo( nxBOOL upper) { m_uplo = upper ? 'U' : 'L';} 
	public:
								LapackSymmetricMatrix	( nxBOOL upper=nxTRUE )	{ Inituplo(upper);}
								LapackSymmetricMatrix   ( int numrows, int numcolumns, double* storage, nxBOOL upper=nxTRUE) : LapackMatrix( numrows, numcolumns, storage) {Inituplo(upper);}

		void					FillOtherHalf			(  );
		char					UpLo					() const				{ return m_uplo;}
		nxBOOL					GetTriangularInverse	( LapackSymmetricMatrix* inverse);
		void					ZeroUnusedTriangle();
		LapackSymmetricMatrix&	operator=( const LapackSymmetricMatrix& other ) { *(LapackMatrix *)this = (LapackMatrix)other; m_uplo = other.m_uplo; return *this;}
};


/*---------------------------------------------------------------------------
 *						class LapackFactorCholesky
 *
 *  Cholesky factorization of a real symmetric positive definite matrix A.
 *
 *  The factorization has the form
 *     A = U**T * U,  if UPLO = 'U', or
 *     A = L  * L**T,  if UPLO = 'L',
 *  where U is an upper triangular matrix and L is lower triangular.
 *
 *-------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
 *					class LapackFactorCholesky						2004-11-23*/
/**	\ingroup Math_Matrix
 *	A specialty class specifically designed to determine and hold the matrix from
 *	Cholesky factorization of a regular matrix.  This is also valuable when solving
 *	matrix equations of the form \b AX=B.
 *
 *  Cholesky factorization of a real symmetric positive definite matrix A.
 *	The factorization has the form
 *  \code
 *		A = U<SUP>T</SUP>*U,		if UPLO = 'U' or
 *		A = L*L<SUP>T</SUP>,		if UPLO = 'L',
 *	\endcode
 *  where U is an upper triangular matrix and L is lower triangular.
**/
/*---------------------------------------------------------------------------*/

class LapackFactorCholesky : public LapackSymmetricMatrix
{
	private:
		nxBOOL				m_isfactored;

	public:
							LapackFactorCholesky		( int numrows, int numcolumns, double* storage, nxBOOL upper=nxTRUE) : LapackSymmetricMatrix( numrows, numcolumns, storage, upper) {m_isfactored = nxFALSE;}
							LapackFactorCholesky		( nxBOOL upper=nxTRUE ) : LapackSymmetricMatrix(upper){ m_isfactored = nxFALSE;}

		nxBOOL				CholeskyFactorize			( const LapackMatrix* A = NULL );
		nxBOOL				GetInverse					( LapackSymmetricMatrix* inverse); 									// LAPACK dgetri
		nxBOOL				LinearSolve					( LapackMatrix& B, LapackMatrix* X);
		void				SetIsFactored( nxBOOL val)  {m_isfactored = val;}
};

/*-----------------------------------------------------------------------------
 *					class LapackGeneralizedSymEigen					2004-11-23*/
/**	\ingroup Math_Matrix
 *  A class used to solve the generalized eigenvalue equation
 *  \code
 *  Ax = lamda*B*x
 *	\endcode
 *	where A is symmetric. In  this example \b B is this class.  It is symmetric
 *	and definite positive as it has been factored by Cholesky.
 **/
/*---------------------------------------------------------------------------*/

class LapackGeneralizedSymEigen : public LapackFactorCholesky
{
	public:
							LapackGeneralizedSymEigen	( int numrows, int numcolumns, double* storage, nxBOOL upper=nxTRUE ) : LapackFactorCholesky( numrows, numcolumns, storage, upper) {}
							LapackGeneralizedSymEigen	( nxBOOL upper=nxTRUE ) : LapackFactorCholesky(upper){}
		nxBOOL				SolveEigenvector			( const LapackMatrix& userA, LapackMatrix* eigenvectors, nx1dArray<double>* eigenvalues );
};

/*---------------------------------------------------------------------------
 *						class LapackFactorSymEigen

 *-------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
 *					LapackSymmetricMatrix		2004-11-23*/
/**	\ingroup Math_Matrix
 *  Reduces a real symmetric matrix A to real symmetric
 *  tridiagonal form T by an orthogonal similarity transformation:\n
 *  Q<SUP>T</SUP>*A*Q = T.
 *	\n
 *	This is then used to solve the Standard eigenvalue problem\n
 *	A*X = lambda*X\n
 *	where A is symmetric
**/
/*---------------------------------------------------------------------------*/

class LapackFactorSymEigen : public LapackSymmetricMatrix
{
	private:
		nx1dArray<double>		m_D;
		nx1dArray<double>		m_E;
		nx1dArray<double>		m_TAU;

	private:
		void				ExtractQ		( LapackMatrix*	Q);
		void				ExtractH		( int i, LapackMatrix* H, nx1dArray<double>& v);

	public:
							LapackFactorSymEigen( nxBOOL upper=nxTRUE ): LapackSymmetricMatrix(upper){}
							LapackFactorSymEigen( int numrows, int numcolumns, double* storage, nxBOOL upper=nxTRUE) : LapackSymmetricMatrix( numrows, numcolumns, storage, upper) {}

		nxBOOL				GetEigenVectors		(LapackMatrix* vectors, nx1dArray<double>* eigenvalues); 							
		nxBOOL				Factorize			(LapackMatrix* matrix = NULL);								//!< LAPACK DSYTRD
};

