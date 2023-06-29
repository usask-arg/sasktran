#include "nxbase_math.h"
#include "nxbase_linearalgebra.h"



/*-----------------------------------------------------------------------------
 *					nxPolynomial::nxPolynomial		2007-7-23*/
/** **/
/*---------------------------------------------------------------------------*/

nxPolynomial::nxPolynomial()
{
}


/*-----------------------------------------------------------------------------
 *					nxPolynomial::~nxPolynomial		2007-7-23*/
/** **/
/*---------------------------------------------------------------------------*/

nxPolynomial::~nxPolynomial()
{
}


/*-----------------------------------------------------------------------------
 *					nxPolynomial::SetCoeffs		2007-7-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxPolynomial::SetCoeffs( nxArrayLinear<double>& coeffs )
{
	m_coeffs = coeffs;
	return (m_coeffs.size() == coeffs.size() );
}


/*-----------------------------------------------------------------------------
 *					nxPolynomial::SetLinear		2007-7-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxPolynomial::SetLinear( double m, double c )
{
	double				buffer[2] = {m,c};
	nx1dArray<double>	buf(2, &buffer[0]);

	m_coeffs = buf;
	return (m_coeffs.size() == 2);
}


/*-----------------------------------------------------------------------------
 *					nxPolynomial::SetQuadratic		2007-7-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxPolynomial::SetQuadratic( double a, double b, double c )
{
	double				buffer[5] = {a,b,c};
	nx1dArray<double>	buf(3, &buffer[0]);

	m_coeffs = buf;
	return (m_coeffs.size() == 3);
}


/*-----------------------------------------------------------------------------
 *					nxPolynomial::SetCoeffsFromPolyFit		2007-7-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxPolynomial::SetCoeffsFromPolyFit	( const nxArrayLinear<double>& x, const nxArrayLinear<double>& y, int norder  )
{
	int					M      = (int)x.size();
	int					N      = norder+1;
	bool				ok;
    nx1dArray<double>	work;
	LapackMatrix		A;
	LapackMatrix		column;
	int					i;

	ok = (M == (int)y.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxPolynomial::SetCoeffsFromPolyFit, The input x and y array are not the same size, I have not made the fit");
	}
	else
	{
		ok = (M >= N);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "nxPolynomial::SetCoeffsFromPolyFit, The input arrays arew too small (%d) to fit a %d order polynomial", (int)x.size(), (int)norder);
		}
		else
		{
			ok =       work.SetSize(M);						// Now set its size down to the number of rows 
			ok = ok && A.SetSize(M,N);						// Set the size of the 
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "nxPolynomial::SetCoeffsFromPolyFit, Error allocating buffers for analysis");
			}
			else
			{
				work.SetTo(1.0);									// Set up the matrix A
				for (i = 0; i < N; i++)								// it has columns of x**N, X**N-1,.... X**2, X, 1 
				{													// So start with the 1 column
					A.Slice( 1, M, N-i, N-i, &column );				// assign the column here
					column = work;
					work *= x;										// NOw point to the next column
				}
				ok = A.Dgels( y, &m_coeffs );
				if (!ok)
				{
					nxLog::Record( NXLOG_WARNING, "nxPolynomial::SetCoeffsFromPolyFit, Error performing linear least squares fit");
				}
			}
		}
	}
	if (!ok) m_coeffs.erase();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxPolynomial::Evaluate		2008-8-1*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxPolynomial::Evaluate(  double x, double* y)
{
	bool	ok;
	size_t	norder;
//	size_t	npts;

	norder = m_coeffs.size();
	ok     = (norder > 0 );
	if (!ok)
	{
		nxLog::Verbose(NXLOG_WARNING, "nxPolynomial::Evaluate, nothing to evalaute as there are no  polynomial coeffiecients are ");
	}
	else
	{
		*y = m_coeffs[0];
		for (size_t i = 1; i < norder; i++)
		{
			*y *= x;
			*y += m_coeffs.At(i);
		}
	}
	if (!ok) *y = 0.0;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxPolynomial::Evaluate		2007-7-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxPolynomial::Evaluate( nxArrayLinear<double>& x, nxArrayLinear<double>* y )
{
	bool	ok;
	size_t	norder;
	size_t	npts;

	norder = m_coeffs.size();
	ok     = (norder > 0 );
	if (!ok)
	{
		nxLog::Verbose(NXLOG_WARNING, "nxPolynomial::Evaluate, nothing to evalaute as there are no  polynomial coeffiecients are ");
	}
	else
	{
		npts = x.size();
		ok   = y->SetSize( 1, &npts );							//!< Set the size and layout of this array
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "nxPolynomial::Evaluate, Error allocating memeory to hold %d points", (int)npts );
		}
		else
		{
			y->SetTo( m_coeffs[0] );
			for (size_t i = 1; i < norder; i++)
			{
				*y *= x;
				*y += m_coeffs.At(i);
			}
		}
	}
	if (!ok) y->erase();
	return ok;
}
