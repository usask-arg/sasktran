#include "nxbase_linearalgebra.h"
#include "sklapack.h"

LapackMatrix operator +( const LapackMatrix& A, double scalar )
{
	LapackMatrix ans( A );
	ans.SetTemporary();
	return (ans += scalar);
}

LapackMatrix operator -( const LapackMatrix& A, double scalar )
{ 
	LapackMatrix ans( A );
	ans.SetTemporary();
	return (ans -= scalar);
}

LapackMatrix operator /( const LapackMatrix& A, double scalar )
{
	LapackMatrix ans( A );
	ans.SetTemporary();
	return (ans /= scalar);
}

LapackMatrix operator*( double scalar, LapackMatrix& other )
{
	LapackMatrix ans( other );
	ans.SetTemporary();
	return (ans *= scalar );
}

LapackMatrix operator *( const LapackMatrix& A, double scalar )
{ 
	LapackMatrix ans( A );
	ans.SetTemporary();
	return (ans *= scalar);
}

LapackMatrix operator *( const LapackMatrix& A, const  LapackMatrix& B)
{
	LapackMatrix C;
	LapackBLAS::Gemm(  1.0, A, ' ',B,' ',&C,0.0);
	C.SetTemporary();
	return C;
}



/*-----------------------------------------------------------------------------
 *					LapackMatrix::Slice		2007-7-20*/
/** **/
/*---------------------------------------------------------------------------*/

LapackMatrix& LapackMatrix::Slice( size_t lo1, size_t hi1, size_t lo2, size_t hi2, LapackMatrix* slice)
{
	size_t	loidx[2] = {lo1-1,lo2-1};
	size_t	hiidx[2] = {hi1-1,hi2-1};

	nxArrayLinear<double>::Slice( loidx, hiidx, 2, slice );
	return *slice;
}

/*---------------------------------------------------------------------------
 *						LapackMatrix::SolveSymmetricEigen
 *	Gets the eigen values and eigen vectors of this matrix assuming it is 
 *	symmetric real.  Uses LAPACK DSYTRD	and LAPACK DSTEDC
 *
 *	Solve Ax = lamda*x
 *
 *	where A is Symmetric and non-singular.
 *-------------------------------------------------------------------------*/

LapackMatrix LapackMatrix::SolveSymmetricEigen( nx1dArray<double>* eigenvalues )
{
	LapackMatrix		 eigenvectors;
	LapackFactorSymEigen factor;
	nxBOOL				 ok;

	ok =       factor.Factorize(this);
	ok = ok && factor.GetEigenVectors(&eigenvectors, eigenvalues);
	if (!ok)
	{
		eigenvectors.Clear();
		eigenvalues->erase();
	}
	eigenvectors.SetTemporary();
	return eigenvectors;
}


/*---------------------------------------------------------------------------
 *						LapackMatrix::LUInverse
 *-------------------------------------------------------------------------*/

LapackMatrix LapackMatrix::LUInverse()
{
	nxBOOL			ok;
	LapackMatrix	inverse; 								// LAPACK dgetrf followed by dgetri
	LapackFactorLU	factor;

	ok =       factor.LUFactorize( *this );							// Factorize this matrix
	ok = ok && factor.GetInverse( &inverse);
	if (!ok) inverse.Clear();
	return inverse;
}

/*---------------------------------------------------------------------------
 *                    LapackMatrix::Transpose                     2019-10-09 */
/** **/
/*---------------------------------------------------------------------------*/

LapackMatrix LapackMatrix::Transpose()
{
	const size_t*	oldstrides = ArrayRankSpecs()->Strides();
	const size_t*	olddims    = ArrayRankSpecs()->Dims();
	size_t			strides[2];
	size_t			dims[2];
	LapackMatrix	transpose;


	strides[0] = oldstrides[1];
	strides[1] = oldstrides[0];
	dims[0]    = olddims[1];
	dims[1]    = olddims[0];
	transpose.nxArrayLinear<double>::Attach( 2, dims, UnsafeArrayBasePtr(), ArrayManager(), strides );
	return transpose;
}

/*-----------------------------------------------------------------------------
 *					LapackMatrix::CreateFromVectorVector		2013-11-21*/
/** Copy the vector<vector> into this Matrix. The code checks to make sure
 *	that each row of the  vector<vector> has the same number of columns.
 **/
/*---------------------------------------------------------------------------*/

bool LapackMatrix::CreateFromVectorVector( const std::vector<vector<double> > &x )
{
	size_t	i;
	size_t	j;
	size_t	ncolumns;
	size_t	nrows;
	bool	ok;

	nrows = x.size();
	ok =  (nrows == 0);
	if (ok)
	{
		Clear();
	}
	else
	{
		ncolumns = x[0].size();
		ok = SetSize( nrows, ncolumns );
		if (ok)
		{
			for (i=0; i < x.size(); i++ )
			{
				ok = ok && ( x[i].size() == ncolumns );
				if (ok)
				{
					for (j=0; j < ncolumns; j++ )
					{
						At(i+1,j+1) = x[i][j];			//lapack matrices are 1 indexed
					}
				}
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "LapackMatrix::CreateFromVectorVector, The number of columns in row of the vector vector was inconsistent. Thats not good");
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"LapackMatrix::CreateFromVectorVector, Error allocating 2D matrix vector with %d rows and %d columns", (int)x.size(), (int)x[0].size());
		Clear();
	}
	return ok;

}

/*-----------------------------------------------------------------------------
 *					LapackMatrix::CreateFromRowVector		2013-11-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool LapackMatrix::CreateFromRowVector( const std::vector<double> &x )
{
	size_t j;
	bool ok;
	
	ok = SetSize( 1, x.size() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"LapackMatrix::CreateFromRowVector, Error allocating row vector with %d rows ", (int)x.size() );
		Clear();
	}
	else
	{
		for (j=0; j<x.size(); j++ )
			At(1,j+1) = x[j];			//lapack matrices are 1 indexed
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					LapackMatrix::CreateFromColumnVector		2013-11-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool LapackMatrix::CreateFromColumnVector( const std::vector<double> &x )
{
	size_t i;
	bool ok;
	
	ok = SetSize( x.size(), 1 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"LapackMatrix::CreateFromColumnVector, Error allocating row vector with %d rows ", (int)x.size() );
		Clear();
	}
	else
	{
		for (i=0; i<x.size(); i++ )
				At(i+1,1) = x[i];			//lapack matrices are 1 indexed
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *						LapackMatrix::CholeskyInverse
 *-------------------------------------------------------------------------*/

LapackSymmetricMatrix LapackMatrix::CholeskyInverse()
{
	nxBOOL					ok;
	LapackSymmetricMatrix	inverse; 								// LAPACK dgetrf followed by dgetri
	LapackFactorCholesky	factor(nxTRUE);

	ok =       factor.CholeskyFactorize( this);							// Factorize this matrix
	ok = ok && factor.GetInverse(&inverse);
	if (!ok) inverse.Clear();
	return inverse;
}

/*---------------------------------------------------------------------------
 *						LapackGeneralizedSymEigen::GeneralizedEigenVector

 *  Solve the generalized eigenvalue equation

	Ax = lamda*B*x

	A is symmetric
	B is this object.  It is symmetric and definite positive. It has been factored by Cholesky)

	See LAPACK users guide: Generalized Symmetric Definite Eigenproblems for more
	specific info.

      SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int            INFO, ITYPE, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYGST reduces a real symmetric-definite generalized eigenproblem
*  to standard form.
*
*  If ITYPE = 1, the problem is A*x = lambda*B*x,
*  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
*
*  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
*  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
*
*  B must have been previously factorized as U**T*U or L*L**T by DPOTRF.
*
*  Arguments
*  =========
*
*  ITYPE   (input) int
*          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
*          = 2 or 3: compute U*A*U**T or L**T*A*L.
*
*  UPLO    (input) CHARACTER
*          = 'U':  Upper triangle of A is stored and B is factored as
*                  U**T*U;
*          = 'L':  Lower triangle of A is stored and B is factored as
*                  L*L**T.
*
*  N       (input) int
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the transformed matrix, stored in the
*          same format as A.
*
*  LDA     (input) int
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
*          The triangular factor from the Cholesky factorization of B,
*          as returned by DPOTRF.
*
*  LDB     (input) int
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) int
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
 *-------------------------------------------------------------------------*/

nxBOOL LapackGeneralizedSymEigen::SolveEigenvector( const LapackMatrix& userA, LapackMatrix* eigenvectors, nx1dArray<double>* eigenvalues )
{
	nxBOOL							ok;
	blas::f77integer				ITYPE = 1;
	blas::f77integer				N;
	blas::f77integer				LDA;
	blas::f77integer				LDB;
	blas::f77integer				INFO;
	blas::f77char					UPLO = UpLo();
	LapackFactorSymEigen			C(UPLO == 'U');

	C.CopyData(userA);
	N   = (int)C.NumColumns();
	LDA = N;
	LDB	= N;
	ok  = ((size_t)N == NumColumns());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "LapackGeneralizedSymEigen::SolveEigenvector, array dimensions mismatch");
	}
	else
	{
		LAPACK(dsygst)( &ITYPE, &UPLO, &N, C.UnsafeArrayBasePtr(), &LDA, UnsafeArrayBasePtr(), &LDB, &INFO );
		ok = (INFO == 0);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"LapackGeneralizedSymEigen::SolveEigenvector,  Parameter %d had an illegal value", (int)-INFO);
		}
		else
		{
			LapackSymmetricMatrix InvL;
			LapackMatrix Y;

			C.Factorize(NULL);									// Now factorize The C matrix
			ok = C.GetEigenVectors( &Y, eigenvalues );			// Get the eigenvector and values from the standard symmetric eigen  problem.
			ok = ok && GetTriangularInverse( &InvL );			// Get the inverse of this (SDP) matrix so we can get the generalized  eigenvectors
/*
			LapackSymmetricMatrix QB;							// Now check that this * Invl = unity
			QB = *this;											// Copy this matrix (which is CHOLESKY factored into L*L**T or U**T*U)
			QB.ZeroUnusedTriangle();							// Zero the half which shoul dbe zero
			LapackMatrix QS;									// another matrix
			QS = InvL*QB;										// do the multiplication
			QS.Dump("QS");										// This is UNITY (it is!)
*/
			if (ok)												// So far so good
			{													// convert the STD eigenvectors to generalized eigen vectors see LAPACK USers Guide for details
				if (UPLO == 'U') ok = LapackBLAS::Gemm( 1.0, InvL,' ', Y, ' ', eigenvectors,  0.0 );
				else             ok = LapackBLAS::Gemm( 1.0, InvL,'T', Y, ' ', eigenvectors , 0.0 );
			}													// and that is that
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"LapackGeneralizedSymEigen::SolveEigenvector,  Error solving standard eigenvalue problem");
			}
		}
	}
	if (!ok)
	{
		eigenvectors->Clear();
		eigenvalues->erase();
	}

	return ok;
}

/*---------------------------------------------------------------------------
 *						LapackMatrix::Factorize
 *	Factorize a symmteric matix using Cholesky factorization

      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOTRF computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U**T * U,  if UPLO = 'U', or
*     A = L  * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the block version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) int
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T.
*
*  LDA     (input) int
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) int
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
 *-------------------------------------------------------------------------*/

nxBOOL LapackFactorCholesky::CholeskyFactorize( const LapackMatrix* A)
{
	if (A != NULL) CopyData(*A);

	blas::f77integer	N    = (int)NumColumns();
	blas::f77integer	LDA  = N;
	blas::f77char		UPLO = UpLo();
	blas::f77integer	info;

	nxBOOL  ok;

	ok  = ((size_t)N == NumRows() && (N > 0));
	if (!ok) 
	{
		nxLog::Record(NXLOG_WARNING, "LapackFactorCholesky::Factorize, the matrix is not square  or is undefined (%dx%d)",(int)NumRows(), (int)NumColumns());
	}
	else
	{
		LAPACK(dpotrf)( &UPLO, &N, UnsafeArrayBasePtr(), &LDA, &info);
		ok = (info == 0);
		if (!ok)
		{
			if (info > 0)
			{
				nxLog::Verbose(NXLOG_WARNING, "LapackMatrix::CholeskyFactorize, The leading minor of order %d is not positive definite", (int)info );
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "LapackMatrix::CholeskyFactorize, The %d'the parameter had an illegal value", (int)-info);
			}
		}

	}
	if (!ok)
	{
		Clear();
	}
	m_isfactored = ok;
	return ok;
}


/*---------------------------------------------------------------------------
 *						LapackFactorCholesky::FillUnusedTriangle
 *	Fills the unused half of the cholesky factorization
 *-------------------------------------------------------------------------*/

void LapackSymmetricMatrix::ZeroUnusedTriangle()
{
	char UPLO = UpLo();

	if (UPLO == 'U')													// Fill the other half of the array
	{																	// with zeroes.
		for (size_t row=2; row <= NumRows(); row++)
		{
			for (size_t col = 1; col <row; col++)At(row,col)=0;
		}
	}
	else
	{
		for (size_t col=2; col <= NumColumns(); col++)
		{
			for (size_t row = 1; row < col; row++) At(row,col)=0;
		}
	}
}


/*---------------------------------------------------------------------------
*						LapackMatrix::LUFactorize
*
*      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
*      int            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
*      int            IPIV( * )
*      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) int
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) int
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) int
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) int array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) int
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
 *-------------------------------------------------------------------------*/


nxBOOL LapackFactorLU::LUFactorize( const LapackMatrix& A)
{

	nxBOOL						ok;
	blas::f77integer			M;
	blas::f77integer			N;
	blas::f77integer			LDA;
	blas::f77integer			INFO;
	blas::f77integer			k;

	m_A.CopyData(A);
	M    = (int)m_A.NumColumns();
	N 	 = (int)m_A.NumRows();
	LDA  = M;
	k    = nxmin(M,N);
	AllocatePivots(k);
	ok = (M >0) && (N > 0);
	if (ok)
	{
		LAPACK(dgetrf)( &M, &N, m_A.UnsafeArrayBasePtr(), &LDA, PivotArray(), &INFO );
		ok = (INFO == 0);
		if (!ok)
		{
			if (INFO < 0)
			{
				nxLog::Record(NXLOG_WARNING, "LapackMatrix::LUFactorize(), DGETRF parameter %d had an illegal value", (int)-INFO);
			}
			else
			{
				nxLog::Verbose(NXLOG_INFO,"LapackMatrix::LUFactorize(), The factorized matrix is singular");
			}
			m_A.Clear();
		}
	}
	m_isfactored = ok;
	return ok;
}

/*---------------------------------------------------------------------------
 *						LapackFactorLU::LULinearSolve
       SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      int            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      int            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) int
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) int
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) int
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) int array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) int
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) int
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*-------------------------------------------------------------------------*/

nxBOOL LapackFactorLU::Solve( const LapackMatrix& B, LapackMatrix* X)
{
	nxBOOL	ok;
	blas::f77char			TRANS = 'N';
	blas::f77integer		N;
	blas::f77integer		NRHS;
	blas::f77integer		LDA;
	blas::f77integer		LDB;
	blas::f77integer		INFO;

	ok = m_isfactored;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "LapackFactorLU::LULinearSolve, Matrix A not yet defined, use method LUFactorize or other Solve method");
	}
	else
	{
		N     = (int)m_A.NumColumns();
		NRHS  = (int)B.NumColumns();
		LDA   = (int)m_A.NumRows();
		LDB   = (int)B.NumRows();
		INFO  = -1;
		X->CopyData( B );
		LAPACK(dgetrs)( &TRANS, &N, &NRHS, m_A.UnsafeArrayBasePtr(), &LDA, m_pivots.UnsafeArrayBasePtr(), X->UnsafeArrayBasePtr(), &LDB, &INFO );
		ok = (INFO == 0);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "LapackFactorLU::Solve, Parameter %d had an illegal value ", (int)-INFO);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					LapackFactorLU::Solve		2007-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

nxBOOL LapackFactorLU::Solve( const LapackMatrix& A, const LapackMatrix& B, LapackMatrix* X)
{
	bool	ok;
	
	ok = LUFactorize( A );
	ok = ok && Solve( B, X );
	return ok;
}




/*---------------------------------------------------------------------------
 *						LapackMatrix::Create
 *-------------------------------------------------------------------------*/

nxBOOL LapackMatrix::Create( size_t numrows, size_t numcols, double initialvalue)
{
	nxBOOL ok;

	ok =SetSize( numrows, numcols );
	SetTo(initialvalue);
	return ok;
}

/*---------------------------------------------------------------------------
 *						LapackMatrix::Initialize
 *-------------------------------------------------------------------------*/

nxBOOL LapackMatrix::Initialize( size_t numrows, size_t numcols, ...)
{
	nxBOOL	ok;
	nxArrayIter<double>	ptr;
	nxArrayIter<double>	ptrend;
	va_list ArgPtr;

	ok = SetSize( numrows, numcols );
	if (ok)
	{
		va_start(ArgPtr, numcols);
		ptr = begin();
		ptrend = end();
		while (ptr < ptrend)
		{
			*ptr = va_arg(ArgPtr, double);
			++ptr;
		}
		va_end(ArgPtr);
	}
	return ok;
}
/*---------------------------------------------------------------------------
 *						LapackMatrix::CreateUnity
 *-------------------------------------------------------------------------*/

nxBOOL LapackMatrix::CreateUnity( size_t numrows )
{
	nxBOOL ok;

	ok = Create( numrows, numrows, 0.0);
	if (ok)
	{
		for (size_t i=1; i <= numrows; i++) (*this)(i,i) = 1.0;
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *						LapackMatrix::Dump
 *-------------------------------------------------------------------------*/

void LapackMatrix::Dump(const char* message)
{
	printf("%s\n", (const char*)message );
	for (size_t row = 1; row <= NumRows(); row++)
	{
		for (size_t col=1; col <= NumColumns(); col++)
		{
			printf("%14.6e ", (double)At(row,col) );
		}
		printf("\n");
	}
}

/*---------------------------------------------------------------------------
 *						
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      int            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      int            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRI computes the inverse of a matrix using the LU factorization
*  computed by DGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) int
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) int
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) int array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) int
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) int
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
 *-------------------------------------------------------------------------*/

nxBOOL LapackFactorLU::GetInverse(LapackMatrix* inverse)
{
	nxBOOL	    ok = false;
	int			N;
	int			LDA;
	double		dummywork;
	double		stackwork[640];
	double*		work;
	int			LWORK;
	int			INFO;
	nxBOOL		allocmem;

	if (!m_isfactored)
	{
		nxLog::Record(NXLOG_WARNING, "LapackFactorLU::GetInverse, cannot get inverse of submitted matrix until you have successfully submitted a matrix to LUFactorize()");
	}
	else
	{
		*inverse = m_A;																	// Copy the pivots over to a new arra
		N    = (int)m_A.NumColumns();
		LDA  = (int)m_A.NumRows();
		ok = (m_pivots.size() >=(size_t)N);												// and make sure we have sensible set of elements
		if (ok)																			// We are ok so far
		{																				// so
			LWORK = -1;																	// first do a workspace query
			LAPACK(dgetri)( &N, *inverse, &LDA, PivotArray(), &dummywork, &LWORK, &INFO );		// do it
			if (INFO < 0) LWORK = 64*N;													// check the status, MKL does not support this query based mode so check for it. and use there recommended size.
			else          LWORK = (int)dummywork;										// Get the optimal value of the size of the WORK arrau
			if (LWORK < N) LWORK = N;													// do a quick sanity check
			allocmem = LWORK > (int)((sizeof(stackwork)/sizeof(stackwork[0])));			// Is our stack space ok to work with
			if (allocmem) work = new double[LWORK];										// Nope so allocate the work memory
			else          work = &stackwork[0];											// otherwise use our stack memory
			ok = (work != NULL);														// make sure everything is ok
			if (ok)																		// it is
			{																			// so
				LAPACK(dgetri)( &N, *inverse, &LDA, PivotArray(), work, &LWORK, &INFO );		// Get the inverse
				if (allocmem) delete [] work;											// and delete the work space if necessary
				ok = (INFO >= 0);														// make sure the status is ok
			}
		}
	}
	if (!ok) inverse->Clear();														// if we had a problem then setthe inverse to zero.
	return ok;																		// return the status.
}
		

/*---------------------------------------------------------------------------
 *						LapackSymmetricMatrix::GetTriangularInverse
 *	Get the inverse of a triangular matrix.  This is done by solving the
 *	The matrix equation AX = B
 *
 *	where:
 *	A = inverse
 *	X = this matrix (triangular)
 *	B = Unit Matrix;
 *-------------------------------------------------------------------------*/

nxBOOL LapackSymmetricMatrix::GetTriangularInverse( LapackSymmetricMatrix* inverse)
{
	nxBOOL ok;
	LapackMatrix	unity;
	unity.CreateUnity(NumColumns() );

	ok = LapackBLAS::Dtrsm( 'R', *this, ' ',unity, inverse, 1.0);	// Now solve InvL*this = I (ie get the inverse 

	/*
	LapackSymmetricMatrix QC;
	QC = *this;
	QC.ZeroUnusedTriangle();
	QC.CopyData( QC*(*inverse) );
	QC.Dump("QC");
	*/
	return ok;
}
/*---------------------------------------------------------------------------
 *						LapackSymmetricMatrix::FillOtherHalf
 *-------------------------------------------------------------------------*/

void LapackSymmetricMatrix::FillOtherHalf()
{
	size_t row;
	size_t	col;

	if (m_uplo == 'U')										// We have been using the upper half
	{														// so we have to fill in the lower half
		for (col=2; col <= NumColumns(); col++)
		{
			for (row = 1; row < col; row++ )
			{
				At(col,row) = At(row,col);
			}
		}
	}
	else													// Otherwise we are using the lower half
	{														// so fill in the upper half
		for (row=2; row <= NumRows(); row++)
		{
			for (col = 1; col < row ; col++ )
			{
				At(col,row) = At(row,col);
			}
		}
	}
}

/*---------------------------------------------------------------------------
 *						LapackFactorCholesky::GetInverse

  SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOTRI computes the inverse of a real symmetric positive definite
*  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
*  computed by DPOTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) int
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T, as computed by
*          DPOTRF.
*          On exit, the upper or lower triangle of the (symmetric)
*          inverse of A, overwriting the input factor U or L.
*
*  LDA     (input) int
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) int
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the (i,i) element of the factor U or L is
*                zero, and the inverse could not be computed.
*
 *-------------------------------------------------------------------------*/

nxBOOL LapackFactorCholesky::GetInverse( LapackSymmetricMatrix* inverse )
{
	int N = (int)NumRows();
	int LDA = (int)NumColumns();
	int info=-1;
	nxBOOL	ok;

	if (!m_isfactored) CholeskyFactorize();
	*inverse = *(LapackSymmetricMatrix *)this;
	LAPACK(dpotri)( &m_uplo, &N, inverse->UnsafeArrayBasePtr(), &LDA, &info );
	ok = (info ==0);
	if (!ok)
	{
		if (info > 0)
		{
			nxLog::Verbose(NXLOG_WARNING, "LapackFactorCholesky::Inverse, The (%d,%d) element is zero and inverse cannot be computed", (int)info, (int)info);
		}
		else
		{
			nxLog::Record(NXLOG_WARNING, "LapackFactorCholesky::Inverse, Parameter %d had an illegal value", (int)-info);
		}
		inverse->Clear(); 
	}
	inverse->FillOtherHalf();				// Complete the other half of the inverse
	return ok;
}

/*---------------------------------------------------------------------------
*						LapackFactorCholesky::LinearSolve

      SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOTRS solves a system of linear equations A*X = B with a symmetric
*  positive definite matrix A using the Cholesky factorization
*  A = U**T*U or A = L*L**T computed by DPOTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) int
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) int
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The triangular factor U or L from the Cholesky factorization
*          A = U**T*U or A = L*L**T, as computed by DPOTRF.
*
*  LDA     (input) int
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) int
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) int
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
 *-------------------------------------------------------------------------*/

nxBOOL LapackFactorCholesky::LinearSolve( LapackMatrix& B, LapackMatrix* X)
{
	int	N;
	int	NRHS;
	int	LDA;
	int	LDB;
	int	INFO;
	nxBOOL	ok;

	N    = (int)NumRows();
	NRHS = (int)B.NumColumns();
	LDB  = (int)B.NumRows();
	LDA  = (int)NumColumns();
	INFO =-1;

	if (!m_isfactored) CholeskyFactorize();
	X->CopyData(B);
    LAPACK(dpotrs)( &m_uplo, &N, &NRHS, UnsafeArrayBasePtr(), &LDA, X->UnsafeArrayBasePtr(), &LDB, &INFO );
	ok = (INFO==0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "LapackFactorCholesky::LinearSolve, Parameter %d had an illegal value", (int)-INFO);
		X->Clear();
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *						LapackFactorSymEigen:Factorize

      SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYTRD reduces a real symmetric matrix A to real symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*  Q**T * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) int
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the orthogonal
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the orthogonal matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) int
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) int
*          The dimension of the array WORK.  LWORK >= 1.
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) int
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
 *-------------------------------------------------------------------------*/


nxBOOL LapackFactorSymEigen::Factorize(LapackMatrix* matrix)						// LAPACK DSYTRD
{
	if (matrix != NULL) CopyData(*matrix);			// Copy the matrix into the destination matrix

	nxBOOL			ok;
	int			M    = (int)NumColumns();
	int			N 	 = (int)NumRows();
	int			LDA  = N;
	int			INFO;
	double			dummywork;
	double*			work;
	double			stackwork[640];
	nxBOOL			allocmem;
	int			LWORK;
	char			UPLO = UpLo();

	ok = (M == N) && (N > 0);
	if (ok)
	{
		m_D.SetSize(N);				// Must be N
		m_E.SetSize(N);				// Only need N-1 here but DSTEGR asks for N
		m_TAU.SetSize(N);			// Only need N-1 here but DSTEGR asks for N
		LWORK = -1;
        
		LAPACK(dsytrd)( &UPLO, &N, UnsafeArrayBasePtr(), &LDA, m_D.UnsafeArrayBasePtr(), m_E.UnsafeArrayBasePtr(), m_TAU.UnsafeArrayBasePtr(), &dummywork, &LWORK, &INFO );		// Try a query to find optimal setting for WORK array
		if (INFO != 0) LWORK = N*64;						// The Intel MKL does not support the query based routine  so just allocate a whole bunch of blocked memory
		else           LWORK = (int)dummywork;			// Otherwise allocate the memory recommended by LAPACK
		if (LWORK < N) LWORK = N;							// Finally make sure everything is sensible

		allocmem = LWORK > (int)((sizeof(stackwork)/sizeof(stackwork[0])));				// Is our stack space ok to work with
		if (allocmem) work = new double[LWORK];										// Nope so allocate the work memory
		else          work = &stackwork[0];											// otherwise use our stack memory
		LAPACK(dsytrd)( &UPLO, &N, UnsafeArrayBasePtr(), &LDA, m_D.UnsafeArrayBasePtr(), m_E.UnsafeArrayBasePtr(), m_TAU.UnsafeArrayBasePtr(), work, &LWORK, &INFO );
		if (allocmem) delete [] work;
		ok = (INFO == 0);
		if (INFO < 0)
		{
			nxLog::Record(NXLOG_WARNING, "LapackMatrix::SymmetricEigenFactorize(), DSYTRD parameter %d had an illegal value", (int)-INFO);
		}
	}
	if (!ok)
	{
		Clear();
		m_D.SetTo(0);
		m_E.SetTo(0);
		m_TAU.SetTo(0);
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *						LapackFactorSymEigen::H
 *
 *	Get the H matrix from the factorized array. This is described in some
 *	detail in DSYTRD documentation.
 *-------------------------------------------------------------------------*/

void LapackFactorSymEigen::ExtractH( int i, LapackMatrix* H, nx1dArray<double>& v)
{
	double			tau;
	int				k;
	int				N = (int)NumColumns();
	char			UPLO = UpLo();

	v.SetSize(N+1);									// Make sure v has sufficient space (use 1 based indexing for ease)
	v.SetTo(0);										// Make sure its completely zero
	H->Create(N,N,0.0);								// Create the matrix and make sure it is zero
	tau  = -m_TAU[i-1];								// Get tau and watch out for zero based indexing

	if (UPLO == 'U') 
	{
		v[i] = 1;
		for (k=1; k < i; k++) v[k] = At(k,i+1);
	}
	else
	{
		v[i+1] = 1;
		for (k=i+2; k <= k; k++) v[k] = At(k,i);
	}

	for (int row=1; row <= N; row++)					// Now do the matrix multiplication
	{
		for (int col=1; col <= N; col++)
		{
			H->At(row,col) = tau*v[row]*v[col];
		}
		H->At(row,row) += 1.0;
	}
}

/*---------------------------------------------------------------------------
 *						LapackFactorSymEigen::GetQ
 *-------------------------------------------------------------------------*/

void LapackFactorSymEigen::ExtractQ(LapackMatrix*	QFINAL)
{
	LapackMatrix	H;
	LapackMatrix    TEMPM;
	LapackMatrix*   TEMP  = &TEMPM;
	LapackMatrix*   Q     = QFINAL;
	LapackMatrix*   dummy;
	int				N = (int)NumColumns();
	nx1dArray<double>	work;
	char			UPLO = UpLo();

	if (UPLO == 'U') 
	{	
		ExtractH( N-1, Q, work );										// Extract the H matrix for element i
		for (int i= N-2; i >= 1; i--)									// Now do the matrix multiplication
		{																// So
			ExtractH( i, &H, work );									// Get the next H
			LapackBLAS::Gemm( 1.0, *Q,' ',H,' ',TEMP, 0.0);	// Get Q = Q*H
			dummy = Q;
			Q = TEMP;
			TEMP = dummy;
		}
	}
	else
	{
		ExtractH( 1, Q, work );
		for (int i=2; i <= N-1; i++)
		{
			ExtractH( i, &H, work );
			LapackBLAS::Gemm( 1.0, *Q,' ',H,' ',TEMP, 0.0);
			dummy = Q;
			Q = TEMP;
			TEMP = dummy;
		}
	}
	if (Q != QFINAL) *QFINAL = *TEMP;
}


/*---------------------------------------------------------------------------
 *						LapackFactorSymEigen::GetEigenVectors

      SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,
     $                   LIWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      int            INFO, LDZ, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      int            IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSTEDC computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the divide and conquer method.
*  The eigenvectors of a full or band real symmetric matrix can also be
*  found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this
*  matrix to tridiagonal form.
*
*  This code makes very mild assumptions about floating point
*  arithmetic. It will work on machines with a guard digit in
*  add/subtract, or on those binary machines without guard digits
*  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
*  It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.  See DLAED3 for details.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*          = 'V':  Compute eigenvectors of original dense symmetric
*                  matrix also.  On entry, Z contains the orthogonal
*                  matrix used to reduce the original matrix to
*                  tridiagonal form.
*
*  N       (input) int
*          The dimension of the symmetric tridiagonal matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the subdiagonal elements of the tridiagonal matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*          On entry, if COMPZ = 'V', then Z contains the orthogonal
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original symmetric matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If  COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) int
*          The leading dimension of the array Z.  LDZ >= 1.
*          If eigenvectors are desired, then LDZ >= max(1,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array,
*                                         dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) int
*          The dimension of the array WORK.
*          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
*          If COMPZ = 'V' and N > 1 then LWORK must be at least
*                         ( 1 + 3*N + 2*N*lg N + 3*N**2 ),
*                         where lg( N ) = smallest int k such
*                         that 2**k >= N.
*          If COMPZ = 'I' and N > 1 then LWORK must be at least
*                         ( 1 + 4*N + N**2 ).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) int array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) int
*          The dimension of the array IWORK.
*          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
*          If COMPZ = 'V' and N > 1 then LIWORK must be at least
*                         ( 6 + 6*N + 5*N*lg N ).
*          If COMPZ = 'I' and N > 1 then LIWORK must be at least
*                         ( 3 + 5*N ).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) int
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  The algorithm failed to compute an eigenvalue while
*                working on the submatrix lying in rows and columns
*                INFO/(N+1) through mod(INFO,N+1).
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Jeff Rutter, Computer Science Division, University of California
*     at Berkeley, USA
*  Modified by Francoise Tisseur, University of Tennessee.
*
 *-------------------------------------------------------------------------*/


nxBOOL LapackFactorSymEigen::GetEigenVectors(LapackMatrix* Z, nx1dArray<double>* eigenvalues)
{

	char				COMPZ = 'V';				// Get eigenvectors of original dense matrix also.
	int					N     = (int)NumRows();
	int					LDZ	  = N;
	int					INFO;
	nx1dArray<double>	E;
	nx1dArray<int>		IWORK;
	nx1dArray<double>	WORK;
	int					LWORK;
	int					LIWORK;
	double				dummywork;
	int					dummyiwork;
	nxBOOL				ok;
	
	*eigenvalues = m_D;							// Copy m_D as it will be overwritten by LAPACK
	E			 = m_E;							// copy m_E as it will be overwritten by LAPACK
	ExtractQ(Z);								// Get the orthogonal matrix used in the reduction to tridiagonal form.
	LIWORK = -1;
	LWORK  = -1;
	INFO   = -1;								// Try and query for an optimal size
	LAPACK(dstedc)( &COMPZ, &N, eigenvalues->UnsafeArrayBasePtr(), E.UnsafeArrayBasePtr(), Z->UnsafeArrayBasePtr(), &LDZ, &dummywork, &LWORK, &dummyiwork, &LIWORK, &INFO );	// Try an optimal size query
	if (INFO != 0)								// IF the query did not work
	{											// eg. Intel MKL 5.0 does not support query calls
		LWORK  = 1 + 3*N + 2*N*32 + 3*N*N;		// So set work sizes greater than indicated in documentation
		LIWORK = 6 + 6*N + 5*N*32;				// for both of the work spaces
	}	
	else										// otherwise the query worked
	{											// so
 		LWORK  = (int)dummywork;				// get the recommended sizes
		LIWORK = (int)dummyiwork;				// and store them
	}

	WORK.SetSize ( LWORK );
	IWORK.SetSize( LIWORK);
	LAPACK(dstedc)( &COMPZ, &N, eigenvalues->UnsafeArrayBasePtr(), E.UnsafeArrayBasePtr(), Z->UnsafeArrayBasePtr(), &LDZ, WORK.UnsafeArrayBasePtr(), &LWORK, IWORK.UnsafeArrayBasePtr(), &LIWORK, &INFO );
	ok = (INFO ==0);
	if (!ok)
	{
		*eigenvalues = m_D;
		if (INFO > 0)
		{
			nxLog::Record( NXLOG_WARNING, "LapackFactorSymEigen::GetEigenVectors  The algorithm failed to compute an eigenvalue while working on the submatrix lying in rows and columns %d/%d through mod(%d,%d)", (int)INFO, (int)(N+1), (int)INFO, (int)(N+1) );
		}
		else
		{
			nxLog::Record( NXLOG_WARNING, "LapackFactorSymEigen::GetEigenVectors, Parameter %d was invalid ", (int)-INFO);
		}
	}
	if (!ok) Z->Clear();
	return ok;
}
