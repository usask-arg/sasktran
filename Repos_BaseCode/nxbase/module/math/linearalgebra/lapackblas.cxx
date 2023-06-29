#include "nxbase_linearalgebra.h"
#include "sklapack.h"
#include <ctype.h>
/*---------------------------------------------------------------------------
 *						LapackBLAS::UpcaseChar
 *-------------------------------------------------------------------------*/

char LapackBLAS::UpcaseChar( char c )
{
	return (char)toupper(c);
}
/*---------------------------------------------------------------------------
 *						LapackBLAS::BLAS_DTRSM
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      int            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - int.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - int.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - int.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - int.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
 *-------------------------------------------------------------------------*/

nxBOOL LapackBLAS::Dtrsm( char								side,
							     const LapackSymmetricMatrix&	A,
								 char							transa,
								 const LapackMatrix&			B,
								 LapackMatrix*					X,
								 double							dalpha)
{
	
	int	N;
	int M;
	int LDA;
	int LDB;
	nxBOOL	ok = nxTRUE;
	char    DIAG  = 'N';
	char	UPLO  = A.UpLo();

	

	N      = (int)B.NumColumns();
	M      = (int)B.NumRows();
	LDB    = M;
	side   = UpcaseChar(side);
	transa = UpcaseChar(transa);
	if (side   != 'L') side   = 'R';				// Make sure we only have limited set
	if (transa != 'T') transa = 'N';

	LDA = (int)A.NumRows();								// Get the number of rows
	ok  = ((size_t)LDA == A.NumColumns());					// and Make sure A is square
	if (ok)											// if it is
	{												// then
		if (side == 'L') ok = (LDA == M);			// Make sure A is compatible with B
		else             ok = (LDA == N);
	}

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "LapackBLAS::DTRSM, dimensions of op(A) and op(B) are incompatible");
		X->Clear();
	}
	else
	{
		X->CopyData(B);						// Copy the B matrix over so it does not get overwritten
		BLAS(dtrsm)( &side, &UPLO, &transa, &DIAG, &M, &N, &dalpha, A.UnsafeArrayBasePtr(), &LDA, X->UnsafeArrayBasePtr(), &LDB );
	}
	return ok;
}

/*---------------------------------------------------------------------------
*						LapackMatrix::BLAS_gemm
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - int.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - int.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - int.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - int.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - int.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - int.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
 *-------------------------------------------------------------------------*/

nxBOOL LapackBLAS::Gemm(		double				dalpha,
									const LapackMatrix&	A,
								    char				transa,
									const LapackMatrix&	B,
									char				transb,
									LapackMatrix*		C,
									double				dbeta)
{
	
	int	m;
	int n;
	int k;
	int lda;
	int ldb;
	int ldc;
	nxBOOL	ok;

	switch (transa )
	{
		case 'T' :
		case 't' : 
		case 'C' :
		case 'c' :	m = (int)A.NumColumns();
					k = (int)A.NumRows();
					break;

		default  :	transa = 'N';
					m = (int)A.NumRows();
					k = (int)A.NumColumns();
					break;
	}

	switch (transb )
	{
		case 'T' :
		case 't' : 
		case 'C' :
		case 'c' :	n  = (int)B.NumRows();
					ok = ((int)B.NumColumns() == k);
					break;

		default  :	transb = 'N';
					n = (int)B.NumColumns();
					ok = (B.NumRows() == (size_t)k);
					break;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "LapackBLAS::GEMM, dimensions of op(A) and op(B) are incompatible");
	}
	else
	{
		lda = (int)A.NumRows();
		ldb = (int)B.NumRows();
		ldc = m;
		if ( ( C->NumRows() != (size_t)m)  || (C->NumColumns() != (size_t)n) ) ok = C->Create( m,n, 0.0);
		if (ok)
		{
			BLAS(dgemm)( &transa, &transb, &m, &n, &k, &dalpha, A, &lda, B, &ldb, &dbeta, *C, &ldc);
		}
	}
	return ok;
}


