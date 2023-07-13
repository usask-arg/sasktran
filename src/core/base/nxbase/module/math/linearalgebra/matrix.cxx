#include "nxbase_linearalgebra.h"

#if defined (NX_WINDOWS)
#pragma message( "matrix.cxx is DEPRECATED. It needs porting so it uses nxLapackMatrix and MKL" )
#endif

#if 0

/*@doc CurveFitting */

/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/


//----------------------------------------------------------------------------
//			NR_PYTHAG
//----------------------------------------------------------------------------

static double NR_PYTHAG( double  a, double b)
{
   double at,bt,ct;

    at=fabs(a);
    bt=fabs(b);
    return ( (at > bt) ? (ct=bt/at,at*sqrt(1.0+ct*ct))
		       : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0) );
}

//----------------------------------------------------------------------------
//			NR_MAX
//----------------------------------------------------------------------------

static double NR_MAX (double maxarg1, double maxarg2 )
{
 return (maxarg1 > maxarg2) ? maxarg1 : maxarg2;
}

//----------------------------------------------------------------------------
//			NR_SIGN
//----------------------------------------------------------------------------

static double NR_SIGN( double a, double b)
{
   return ( b >= 0.0 ? fabs(a) : -fabs(a));
}

//----------------------------------------------------------------------------
//			nrerror
//----------------------------------------------------------------------------

static void nrerror( char *error_text)
{
	nxLog::Record(NXLOG_WARNING,"MatrixArray::nrerror, Numerical Recipes error [%s]",error_text);
}


//----------------------------------------------------------------------------
//			NR_vector
//----------------------------------------------------------------------------

static double *NR_vector( int nh)
{
   double *v;

   v=(double *)malloc((unsigned) (nh+1)*sizeof(double));
   if (!v)
   {
      nrerror("NR_Vector, allocation failure in NR_vector()");
      return NULL;
   }
   return v;
}

//----------------------------------------------------------------------------
//			NR_free_vector
//----------------------------------------------------------------------------

static void NR_free_vector( double *v )
{
	free((char*) (v));
}


//-----------------------------------------------------------------------------
//			svdcmp
//
//	Given a matrix a[0..m-1, 0..n-1] compute singular value decomposition
//	a=u.w.vT.  The matrix u replaces matrix a on output.  The diagonal
//	matrix of singular values w is output as a vector  w[0..n-1].  The
//	matrix v (not the transpose vT) is output as v[0..n-1,0..n-1].
//
//     	Upon Input
//	----------
//	a	double a[][],  MxN matrix to decompose
//	m	number or rows     (ie vertical size of array)
//	n	number of colums   (ie horizontal size of array)
//	w	double w[N]        matrix for storing diagonal weights.
//	v	double v[N][N]     matrix for storing v (not transpose v)
//
//
//	On Output
//	---------
//	a	double u[M][N]
//
//	If there is any problem with the singular value decomposition then
//	it returns nxFALSE, else returns TRUE.
//
//	Numerical recipes, Chapter 2.9, Singular Value Decomposition.
//----------------------------------------------------------------------------

static nxBOOL NR_svdcmp( double **a, int m, int n, double *w, double **v )
{
   int flag,i,its,j,jj,k;
   int l = 0;
   int nm = 0;
   double c,f,h,s,x,y,z;
   double anorm=0.0,g=0.0,scale=0.0;
   double *rv1;

	if (m < n)
	{
		nrerror("SVDCMP: More columns than rows!, You must augment A with extra zero rows");
		return nxFALSE;
	}
	rv1=NR_vector(n);
	if (rv1 == NULL) return nxFALSE;

	for (i=0; i < n; i++)
	{
		l=i+1;
		rv1[i]=scale*g;
		g = s = scale = 0.0;
		if (i < m)
		{
			for (k=i; k < m; k++) scale += fabs(a[k][i]);
			if (scale)
			{
				for (k=i;k < m;k++)
				{
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f = a[i][i];
				g = -NR_SIGN(sqrt(s),f);
				h = f*g-s;
				a[i][i]=f-g;
				if (i != (n-1))
				{
					for ( j=l; j < n; j++)
					{
						for (s=0.0,k=i;k < m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i; k<m; k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if ( (i < m) && (i != (n-1)) )
		{
			for (k=l;k<n;k++) scale += fabs(a[i][k]);
			if (scale)
			{
				for (k=l;k<n;k++)
				{
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f = a[i][l];
				g = -NR_SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
				if (i != (m-1) )
				{
					for (j=l;j<m;j++)
					{
						for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=NR_MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--)
	{
		if (i < n-1)
		{
			if (g)
			{
				for (j=l;j<n;j++)
				v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++)
				{
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n-1;i>=0;i--)
	{
		l=i+1;
		g=w[i];
		if (i < n-1) for (j=l;j<n;j++) a[i][j]=0.0;
		if (g)
		{
			g=1.0/g;
			if (i != (n-1))
			{
				for (j=l;j<n;j++)
				{
					for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		}
		else
		{
			for (j=i;j<m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n-1;k>=0;k--)
	{
		for (its=1;its<=30;its++)
		{
			flag=1;
			for (l=k;l>=0;l--)
			{
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm)
				{
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag)
			{
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++)
				{
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=NR_PYTHAG(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s=(-f*h);
					for (j=0;j<m;j++)
					{
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k)
			{
				if (z < 0.0)
				{
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30)
			{
				nrerror("NR_svdcmp, No convergence in 30 SVDCMP iterations");
				NR_free_vector( rv1);
				return nxFALSE;
			}
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=NR_PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+NR_SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++)
			{
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=NR_PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=0;jj<n;jj++)
				{
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=NR_PYTHAG(f,h);
				w[j]=z;
				if (z)
				{
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=0;jj<m;jj++)
				{
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	NR_free_vector(rv1);
	return nxTRUE;
}

//----------------------------------------------------------------------------
//			NR_svbksb
//	Solves A.x = b for a vector X, where A is specified by the arrays,
//	u[0..m-1, 0..n-1], w[0..n-1], v[0..n-1, 0..n-1] as returned by svdcmp.
//	m and n are dimensions of A and will be equal for square matrices.
//	b[0..m-1] is the input right-hand side. x[0..n-1] is the output
//	solution vector. No input quantities are destroyed so the routine may be
//	called sequentially with different b's.
//
//	Note the caller should modify the singular diagonal elements before
//	calling this routine, if you wish to use the full power of singular
//	value decomposition
//
//	Numerical Recipes, Chapter 2.9, Singular Value Decomposition.
//-----------------------------------------------------------------------------


static nxBOOL NR_svbksb( double *u[], double w[], double *v[], int m, int n, double b[], double x[])
{
   int jj,j,i;
   double s,*tmp;

   tmp= NR_vector(n);
   if (!tmp)
   {
      nrerror("NR_svbksb, Unable to allocate temporary array");
      return nxFALSE;
   }
   for (j=0; j < n;j++)
   {
      s=0.0;
      if (w[j])
      {
         for (i=0;i < m;i++) s += u[i][j]*b[i];
	 s /= w[j];
      }
      tmp[j]=s;
   }
   for (j=0;j < n; j++)
   {
      s=0.0;
      for (jj=0; jj < n;jj++) s += v[j][jj]*tmp[jj];
      x[j]=s;
   }
   NR_free_vector(tmp);
   return nxTRUE;
}

//----------------------------------------------------------------------------
//			MaxtrixArray::default Constructor
//@mfunc
//	Default constructor for the MatrixArray class
//----------------------------------------------------------------------------

MatrixArray::MatrixArray()
{
}

//----------------------------------------------------------------------------
//			MatrixArray::Constructor 1
// @mfunc
//	Create a new array, we allocate this memory so set flag to delete
//	when the destructor is called.
//----------------------------------------------------------------------------

MatrixArray::MatrixArray(  
						   const int w,			// @parm width of array (X Dimension) 
						   const int h)			// @parm height of array (Y Dimension)
{
   SetSize( w,h );
}

//----------------------------------------------------------------------------
//			MaxtrixArray::Copy Constructor
// @mfunc
//	The copy constructor.  Make a complete copy of the other object.
//----------------------------------------------------------------------------

MatrixArray::MatrixArray(
						 MatrixArray& copy )
{
   *this = copy;
}


//----------------------------------------------------------------------------
//			MatrixArray::operator[] (row vector)
//
//	@mfunc
//		Returns the i'th row of the row major matrix
//
//----------------------------------------------------------------------------

nxArray<double>&	MatrixArray::operator [] (
											  int i )		// @parm Index of row to return, 0 based index
{
	return nx2DArray<double>::operator[](i);
}
		
//----------------------------------------------------------------------------
//			MatrixArray::operator=	(x,y)
// @mfunc
//    Assignment operator.  Replaces the matrix with the 2D array passed in.
//----------------------------------------------------------------------------

MatrixArray& MatrixArray::operator= ( 
									 const nx2DArray<double>& RHValue ) // @parm The array to copy
{
	nx2DArray<double>::operator= ( RHValue );
	return *this;
}

/*---------------------------------------------------------------------------
 *						MatrixArray::isColumnVector
 @mfunc
  Returns nxTRUE if the object is a 1-D vector
 *-------------------------------------------------------------------------*/

nxBOOL MatrixArray::isColumnVector()
{
	return (XSize() == 1) || (YSize() == 1 );
}

//----------------------------------------------------------------------------
//			MatrixArray::Transpose
// @mfunc
// calculates the transposes of this array and returns the result
// in transpose.
//
// @rdesc
//	Returns nxTRUE if successful.
//----------------------------------------------------------------------------

nxBOOL MatrixArray::Transpose(
							   MatrixArray *transposed )	// @parm returns the transpose of this array
{
    int i,j;
	nxBOOL	ok;

	ok = transposed->SetSize( YSize(), XSize() );
	if (ok)
	{
		for( i=0; i < YSize(); i++)
		{
			for (j=0; j < XSize(); j++)
			{
				transposed->At(i,j) = At(j,i);
			}
		}
	}
	return ok;
}

//----------------------------------------------------------------------------
//			MatrixArray::SVD
//
//@mfunc
//	performs Singular Value decomposition on the current array
//	and returns the SVD in three arrays u, w, v
//
//	w is a one-dimensional array representing the diagonal elements
//	of a square array.  I have not expanded it out so it is compatible
//	with the SVDSolveLinear which expects "w" to be one-dimensional.
//
//	Uses a modified form of the Numerical recipes svdcmp routine to
//	to the actual decomposition.
//
//@rdesc
//	Returns nxTRUE if successful.
//----------------------------------------------------------------------------

nxBOOL MatrixArray::SVD(
						 MatrixArray* u,	// @parm returns u array from SVD
						 MatrixArray* w,	// @parm returns w array from SVD, returned as a 1-D array
						 MatrixArray* v)	// @parm returns v array from SVD
{
   int			m,n;
   nxBOOL		ok;
   double**		uptr;
   double**		vptr;
   double*	    wptr;
   int			i;

	if (isNull())												// test to make sure
	{															// member array is defined
		nxLog::Record(NXLOG_WARNING, "MatrixArray::SVD, member matrix is undefined");	// error.
		return nxFALSE;
	}
 	m  = YSize();											// get the number of rows
	n  = XSize();											// and number of columns

	if ( m < n )											// columns
	{														// because svdcmp bitterly complains
		nxLog::Record(NXLOG_WARNING, "MatrixArray::SVD, member matrix cannot have more columns than rows");
		return nxFALSE;
	}

	u->operator=(*this);										// assign this array to u for SVD input
	w->SetSize( n, 1 );											// make the 1-D vector for w
	v->SetSize( n, n );											// make the nxn square array for v
	uptr = new double* [m];										// allocate space to buffer the double** pointers for u
	vptr = new double* [n];										// and v
	wptr = (double *)((*w)[0]);
	for (i=0;i < m; i++) uptr[i] = (double*)((*u)[i]);			// copy the pointers over
	for (i=0;i < n; i++) vptr[i] = (double*)((*v)[i]);			// using the type casting

	ok = NR_svdcmp( uptr, m, n, wptr, vptr );						// do the SVD
	if (!ok)													// and check the status.
	{
		nxLog::Record( NXLOG_WARNING,"MatrixArray::SVD. Numerical recipes algorithm failed");
		w->Empty();
		v->Empty();
		u->Empty();
	}
	delete [] vptr;
	delete [] uptr;
	return ok;
}

//----------------------------------------------------------------------------
//			MatrixArray::SVDSolveLinear
//
//@mfunc
//	Solves the linear equation A.x = b where "this" matrix is b.
//	Three matrices, u,w,v are input from the member function SVD.
//
//	Uses a modified form of the Numerical recipes svdcmp routine to
//	do the actual decomposition.
//
//@rdesc
//	Returns nxTRUE if success.
//----------------------------------------------------------------------------

nxBOOL MatrixArray::SVDSolveLinear(
								    MatrixArray* x,		// @parm x matrix in AX=B
									MatrixArray &u,		// @parm u array from SVD of A
									MatrixArray &w,		// @parm w array from SVD of A
									MatrixArray &v)		// @parm v array from SVD of A
{
   int			m,n;
   double**		uptr;
   double**		vptr;
   nxBOOL		ok;
   int			i;

   m = u.YSize();
   n = u.XSize();

   ok =      ( m != 0) && (n != 0 ) && ( m >= n)					// ensure we have proper dimensions.
         &&  ( isColumnVector()     && ( N_Elements() == m))		// ensure b is  b[0..m-1] b is "this"
		 &&  ( w.isColumnVector()   && ( w.N_Elements() == n))		// ensure w is  w[0..n-1]
	     &&  ((v.XSize() == n)      && (v.YSize()== n));			// ensure v is  v[0..n-1][0..n-1]

   if (!ok)
   {
		nxLog::Record( NXLOG_WARNING, "MatrixArray::SVDSolveLinear, Matrix sizes are inconsistent for analysis");
		x->Empty();
		return nxFALSE;
   }

   x->SetSize( n );													// Set up a buffer to receive the solution
   uptr = new double* [m];											// Create a buffer to pass over the twod array
   vptr = new double* [n];											// Create another buffer for the other two-D array
   for (i=0; i < m; i++) uptr[i] = u[i];							// Copy the pointers into the array
   for (i=0; i < n; i++) vptr[i] = v[i];							// and do same for this
   ok = NR_svbksb( uptr, w, vptr, m, n, *this, *x );					// call the SVD solution
   delete [] vptr;													// delete temporary allocation
   delete [] uptr;													// for both arrays

   if (!ok)
   {
		nxLog::Record( NXLOG_WARNING, "MatrixArray:SVDSolveLinear, Error in Numerical Recipes svbksb");
		x->Empty();
   }
   return ok;
}

//----------------------------------------------------------------------------
//			MatrixArray::SquarefromDiagonal
//
//@mfunc
//	Takes a 1-D array whose elements represent the diagonal of a square
//	matrix, and makes that matrix.  Useful for SVD where array "w" is
//	a diagonal nxn matrix returned as a 1-D array.
//
//@rdesc
//	Returns TRUE if success.
//----------------------------------------------------------------------------

/*<----------------------------------------------------------------------------
 *'			MatrixArray::SquarefromDiagonal
 *	Takes a 1-D array whose elements represent the diagonal of a square
 *  matrix, and makes that matrix.  Useful for SVD where array "w" is
 *	a diagonal nxn matrix returned as a 1-D array.
 *
 *''RETURNS:
 *	nxBOOL, Returns TRUE if success.
 *>----------------------------------------------------------------------------*/



nxBOOL MatrixArray::SquareFromDiagonal(
									   MatrixArray* square) // @parm Returns the square matrix.
{
	nxBOOL	ok;
	if (!isColumnVector())
	{
	   nxLog::Record(NXLOG_WARNING, "MatrixArray::SquareFromDiagonal, Matrix is not a column vector");
	   square->Empty();
	   ok = nxFALSE;
	}
	else
	{
		int n = N_Elements();
		ok = square->SetSize(n,n);
		square->SetTo(0);
		for (int i=0; i < n; i++) square->At(i,i) = m_array[i];
	}
	return ok;
}

//----------------------------------------------------------------------------
//			MatrixArray::Multiply
//@mfunc
//	Performs matrix multiplication.  answer = [this]*[other]
//
//@rdesc
//	Returns TRUE if success.
//----------------------------------------------------------------------------


nxBOOL MatrixArray::Multiply(
							  MatrixArray& other,		// @parm the [other] matrix.
							  MatrixArray* answer )		// @parm Returns the answer of [this]*[other]
{
	nxBOOL				ok;
	int					row;
	int					col;
	int					k;
	double				sum;

	ok = (NumColumns() == other.NumRows());							// Make sure the arrays are of compatible size
	if (ok)															// If the are then proceed
	{																// so  (l == Number of rows, n = Number of columns )
		answer->SetSize( other.NumColumns(), NumRows() );					// allocate space for this array  A(l x m)*B(m x n)  ==> C( l x n )
		int otherstride = other.NumColumns();						// Get the stride between points
		for (row = 0; row < NumRows(); row++)							// Now iterate over the rows 
		{
			nxArray<double>& Arow = (*this)[row];					// Get the pointer to this row
			for (col = 0; col < other.NumColumns(); col++)                        
			{
				sum = 0.0;
				double*	Aptr =  (double *)Arow;
				double* Bptr = ((double *)other) + col;
				for ( k = 0;  k < NumColumns(); k++ )         
				{
					sum  += (*Aptr++) * (*Bptr);					// Get the summation
					Bptr += otherstride;							// Point to the next element
				}
				answer->At( col, row) = sum;
			}
		}
	}                                            
   return ok;
}

#endif




