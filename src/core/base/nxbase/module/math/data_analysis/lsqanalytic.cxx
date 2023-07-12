
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_math.h"
#include "nxbase_linearalgebra.h"

#if defined(NX_WINDOWS)
#pragma message( "lsqanalytic.cxx is DEPRECATED. It needs porting so it uses nxLapackMatrix" )
#endif

#if 0

//---------------------------------------------------------------------------
//						LSQ_ANALYTIC_CHI_SQUARE
// @func
// Deprecated function.  Calculates chi square.  It will eventually
// disappear from the library
//---------------------------------------------------------------------------


 double  LSQ_ANALYTIC_CHI_SQUARE( const nxdblarr& y, const nxdblarr& yfit, int lim1, int lim2 ) 
{ 
	int          I;
	double       X;

	X = 0; 
	for (I = lim1; I <= lim2; I++)
	{
		X += nxmath::sqr(y[I] - yfit[I]); 
	}
	return X; 
} 


//---------------------------------------------------------------------------
//					LSQ_ANALYTIC_MATRIX_INVERT
// @func
// Solves the matrix equation necessary to obtain solutions for
// the fit TOPPER. The routine uses the class MatrixArray to solve
// the matrix equation.
//
// A deprecated function. It will eventually disappear 
//
//---------------------------------------------------------------------------

 nxBOOL  LSQ_ANALYTIC_MATRIX_INVERT(	const nxdblarr&			 BETA,
											const nx2dArray<double>& ALPHA,
											double					 lamda,
											nxdblarr*				 NEWA,
											const nxdblarr&			 OLDA)
{ 

	int          n = BETA.N_Elements();
	MatrixArray  A(n,n);
	MatrixArray  V(n,1);
	MatrixArray  X(n,1);
	MatrixArray	 u;
	MatrixArray	 v;
	MatrixArray  w;
	double       FACT;
	nxBOOL		 ok;

	for (int J = 0; J < n; J++) 
	{ 
		V[0][J] = BETA[J]; 
		for (int K = 0; K < n; K++) 
		{ 
			if (K == J) FACT = (1 + lamda); 
			else        FACT = 1; 
			A[K][J] = ALPHA[K][J] * FACT; 
		} 
	} 

	ok  = A.SVD( &u, &w, &v );
	ok  = ok && V.SVDSolveLinear( &X, u,w,v );
	*NEWA  = X[0];
	*NEWA += OLDA;
	return ok;
}

//---------------------------------------------------------------------------
//				ALPHA_AND_BETA
// @func
//	A deprecated function.  Calculates the arrays associated with linear
//	least square fitting using derivative expansion about a solution.
//---------------------------------------------------------------------------

 void   ALPHA_AND_BETA( nx2DArray<double>&	 derivs,
							 nxdblarr*	         BETA,
							 nx2DArray<double>*	 ALPHA,
							 const nxdblarr&	 Y,
							 const nxdblarr&     yfit,
							 const nxdblarr&	 parameters,
							 int				 lim1,
							 int				 lim2)
{
	int         I;
	double      DEL;
	int			J, K;
	double*		deriv;
	double		dj;
	int			np = parameters.N_Elements();

	BETA->SetSize(np);
	ALPHA->SetSize(np,np);
	BETA->SetTo(0.0);
	ALPHA->SetTo(0.0);

	for (I = lim1; I <= lim2; I++)											// scan over valid data
	{																			// for each data point 
		DEL = Y[I] - yfit[I];													// get local difference 
		deriv = derivs[I];
		for (J = 0; J < np; J++)												// with each set of derivs
		{																		// increment the
			dj           = deriv[J];
			BETA->At(J) += (DEL * dj);											// beta array
			for (K = 0; K <= J; K++)											// and then do the
			{																	// alpha array
				ALPHA->At(K,J) += (dj * deriv[K]); 
			}																	// do across alpha row 
		}																		// do down matrix columns
	}																			// do all data points
	for (J = 0; J < np; J++)													// with each set of derivs
	{																			// symmetric so set all
		for (K = 0; K <= J; K++) ALPHA->At(J,K) = ALPHA->At(K,J);				// values at end
	}																			// of calculation
} 
#endif
