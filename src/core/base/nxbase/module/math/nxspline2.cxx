#include "nxbase_math.h"


/*-----------------------------------------------------------------------------
 *					nxSpline2::nxSpline2		2002-10-7						*/
/** Default Constructor **/
/*---------------------------------------------------------------------------*/

nxSpline2::nxSpline2()
{
	m_padvalue = std::numeric_limits<double>::quiet_NaN();
	m_x.SetReuseMemory(true);
	m_y.SetReuseMemory(true);
	m_y2.SetReuseMemory(true);
	m_u.SetReuseMemory(true);
}


/*---------------------------------------------------------------------------
 *                      nxSpline2::nxSpline2                      2020-01-13 */
/** **/
/*---------------------------------------------------------------------------*/

nxSpline2::nxSpline2( const nxSpline2& other)
{
	m_padvalue = other.m_padvalue;
	m_x.SetReuseMemory(true);
	m_y.SetReuseMemory(true);
	m_y2.SetReuseMemory(true);
	m_u.SetReuseMemory(true);
	DeepCopy(other);
}

/*-----------------------------------------------------------------------------
 *					nxSpline2::~nxSpline2								2002-10-7*/
/** Destructor **/
/*---------------------------------------------------------------------------*/

nxSpline2::~nxSpline2()
{
}

/*---------------------------------------------------------------------------
 *'					nxSpline2::ReleaseResources                      2002-10-7
 *-------------------------------------------------------------------------*/

void nxSpline2::ReleaseResources()
{
	m_x.SetSize(0);
	m_y.SetSize(0);
	m_y2.SetSize(0);
	m_u.SetSize(0);
}

/*---------------------------------------------------------------------------
 *'					nxSpline2::Allocate                              2002-10-7
 *-------------------------------------------------------------------------*/

nxBOOL nxSpline2::Allocate(size_t ncoefs )
{
	nxBOOL	ok = nxTRUE;

	ok = ok && m_x.SetSize( ncoefs );
	ok = ok && m_y.SetSize( ncoefs );
	ok = ok && m_y2.SetSize( ncoefs );
	ok = ok && m_u.SetSize( ncoefs );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxSpline2::Allocate, Error allocating %d elements for the spline curve", (int)ncoefs);
		ReleaseResources();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxSpline2::DeepCopy		2008-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::DeepCopy( const nxSpline2& other )
{
	bool	ok;

	ok =       m_x.DeepCopy( other.m_x);
	ok = ok && m_y.DeepCopy( other.m_y);
	ok = ok && m_y2.DeepCopy( other.m_y2);
	ok = ok && m_u.DeepCopy( other.m_u);
	m_padvalue			= other.m_padvalue;
	m_startsplineindex  = other.m_startsplineindex;
	m_endsplineindex    = other.m_endsplineindex;

	return ok;
}



/*---------------------------------------------------------------------------
 *                  nxSpline2::CheckIsAscending                   2020-01-13 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::CheckIsAscending(const double* x,  size_t n )
{
	bool	ok = true;
	double  lastx = x[0];

	for (size_t i = 0; i < n; i++)
	{
		ok = ok && (x[i] >= lastx);
		lastx = x[i];
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxSpline2::CheckIsAscending, the input x-array is not in ascending order. This creates problems and is not supported.");
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *                  nxSpline2::CheckNotBadValues                  2020-01-13 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::CheckNotBadValues(const double* y,  size_t n )
{
	bool	ok = true;

	for (size_t i = 0; i < n; i++)
	{
		ok = ok && ( (y[i] != m_padvalue) && isfinite(y[i]) );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxSpline2::CheckNotBadValues, the input array contains either bad values (%e) or NaNs. This creates problems and is not supported.", (double)m_padvalue);
	}
	return ok;
}



/*---------------------------------------------------------------------------
 *                       nxSpline2::CopyXY                        2020-01-14 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::CopyXY( const double*x, const double* y,  size_t n )
{
	for (size_t i = 0; i < n; i++)
	{
		m_x.at(i) = x[i];
		m_y.at(i) = y[i];
	}
	m_y2.SetTo(0.0);
	m_u.SetTo(0.0);
	return true;
}


/*---------------------------------------------------------------------------
 *                 nxSpline2::GetStartSplineIndex                 2020-01-13 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::GetStartSplineIndex(const double* y,  size_t n)
{
	bool	same = true;
	double  firstvalue = y[0];
	size_t	i;
	
	m_startsplineindex = m_y.size();

	for (i = 1; same && (i < n); i++)																		// Scan forwards through the array
	{																										// and 
		same = (y[i] == firstvalue) || ( std::isnan	(firstvalue) && std::isnan(y[i]));						// find if this element is the same as the first element
		if (!same)																							// if its not then we have found the start of the array
		{																									// so if we have the first set of constant values are the pad value then start at first non-pad value
			m_startsplineindex = (i > 1) ? (i): (i-1);													// else start at the element efore
		}
	}
	return true;
}

/*---------------------------------------------------------------------------
 *                  nxSpline2::GetEndSplineIndex                  2020-01-13 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::GetEndSplineIndex   (const double* y,  size_t n)
{
	bool	same = true;
	double  lastvalue = y[n-1];
	size_t	i;
	
	m_endsplineindex = m_y.size();

	for (i = n-1; same && (i >= m_startsplineindex); i--)														// Scan backwards through the array
	{																										// and 
		same = (y[i] == lastvalue) || ( std::isnan	(lastvalue) && std::isnan(y[i]));						// find if this element is the same as the last element
		if (!same)																							// if its not then we have found the start of the array
		{																									// so if we have the first set of constant values are the pad value then start at first non-pad value
			m_endsplineindex = ( i < n-2) ? (i+1): (i+2);													// else start at the element efore
		}
	}
	return true;

}


/*---------------------------------------------------------------------------
 *                     nxSpline2::CheckBounds                     2020-01-13 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::CheckBounds(const double* x,  const double* y, size_t n)
{
	bool	ok;

	ok = (n > 0);
	ok = ok && CheckIsAscending   ( x, n);
	ok = ok && CheckNotBadValues  ( y, n);
	ok = ok && GetStartSplineIndex( y, n);
	ok = ok && GetEndSplineIndex  ( y, n);
	ok = ok && CopyXY( x,y,n);
	return ok;
}


/*---------------------------------------------------------------------------
 *               nxSpline2::ConfigureSplineSegment                2020-01-13 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::ConfigureSplineSegment( double* x,  double* y, double* y2, double* u, size_t n, double yp0, double ypn)
{
	int					i,k;
	double				p,qn,sig,un;			//,*u,*vector();
	bool				ok;
	double				xm1;
	double				xp1;
	double				xi;
	double				yp1;
	double				ym1;
	double				yi;
	double				um1;
	double				xp1mxm1;
	double				ximxm1;
	double				ui;

	ok = ( n > 2);
	{
		if (yp0 >= 1.0e30)
		{
			y2[0] = 0.0;
			u [0] = 0.0;
		}
		else
		{
			y2[0] = -0.5;
			u [0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp0);
		}

		um1 = u[0];
		xm1 = x[0];
		ym1 = y[0];
		xi  = x[1];
		yi  = y[1];
		for (i=1; i< n-1; i++)
		{
			xp1     = x[i+1];
			yp1     = y[i+1];
			xp1mxm1 = xp1-xm1;
			ximxm1  = xi-xm1;
//			sig     = (x[i] - x[i-1])/ (x[i+1]-x[i-1]);									// THIS CAN GO
			sig     = ximxm1/xp1mxm1;
			p       = sig*y2[i-1]+2.0;
			y2[i]   = (sig-1.0)/p;
			ui      = (yp1   - yi )/(xp1-xi) - (yi-ym1)/(ximxm1);
			ui      = (6.0*ui/(xp1mxm1)-sig*um1)/p;
			u[i]    = ui;
			um1     = ui;
			xm1     = xi;
			ym1     = yi;
			xi      = xp1;
			yi      = yp1;
		}

		if (ypn >= 1.0e30)
		{
			qn = 0.0;
			un = 0.0;
		}
		else
		{
			qn = 0.5;
			un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
		}
		y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
		for (k=(int)n-2; k >=0; k--)
		{
			y2[k]= y2[k]*y2[k+1]+u[k];
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxSpline2::Configure								2002-10-7*/
/**	Calculates the coefficients for future spline interpolation. Also caches the
 *	x and y data points internally. Based upon the code given in Numerical Recipes 
 *
 *
 *	\param x
 *		The x coordinates of the n tabulated data points. Unchanged
 *
 *	\param y
 *		The y coordinates of the n tabulated data points. Unchanged
 *			
 *	\param n
 *		The number of tabulated data points
 *
 *	\param yp0
 *		The gradient of spline at first point. If > 1.0E30 then use a natural spline at start
 *
 *	\param ypn
 *		The gradient of spline at last point. If > 1.0E30 then use a natural spline at start
 *			
 *
 *	\return 
 *		nxTRUE if successful
 *
 **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::Configure(const double* x,  const double* y, size_t n, double yp0, double ypn)
{
	bool				ok;
	size_t				numspline;

	ok = ( n > 1);
	if (ok)
	{
		ok = Allocate( n );
		ok = ok && CheckBounds(x, y, n);
		if (ok)
		{
			numspline =  m_endsplineindex - m_startsplineindex;
			ok = ConfigureSplineSegment( m_x.UnsafeArrayBasePtr()  + m_startsplineindex,
										 m_y.UnsafeArrayBasePtr()  + m_startsplineindex,
										 m_y2.UnsafeArrayBasePtr() + m_startsplineindex,
										 m_u.UnsafeArrayBasePtr()  + m_startsplineindex,
										 numspline,			yp0, ypn);
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxSpline2::Configure								2002-10-7*/
/**	Calculates the coefficients for future spline interpolation. Also caches the
 *	x and y data points internally. Based upon the code given in Numerical Recipes 
 *
 *
 *	\param x
 *		The x coordinates of the n tabulated data points. Unchanged
 *
 *	\param y
 *		The y coordinates of the n tabulated data points. Unchanged
 *			
 *	\param n
 *		The number of tabulated data points
 *
 *	\param yp0
 *		The gradient of spline at first point. If > 1.0E30 then use a natural spline at start
 *
 *	\param ypn
 *		The gradient of spline at last point. If > 1.0E30 then use a natural spline at start
 *			
 *
 *	\return 
 *		nxTRUE if successful
 *
 **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::Configure( const nx1dArray<double>& x,  const nx1dArray<double>& y, double yp0, double ypn)
{
	const double*	xptr =x.ArrayBasePtr();
	const double*	yptr =y.ArrayBasePtr();
	bool			ok;

	ok = x.ArrayRankSpecs()->IsContiguous() && (y.ArrayRankSpecs()->IsContiguous());
	ok = ok && (x.size() == y.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxSpline2::Configure for nx1dArray only works for contiguous nxArrays of the same size");
	}
	else
	{
		ok = Configure( xptr,  yptr, x.size(), yp0, ypn );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxSpline2::Configure		2013-11-1*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline2::Configure( const std::vector<double>& userx,  const std::vector<double>& usery, double yp0, double ypn)
{
	bool				ok;

	ok =  (userx.size() == usery.size());
	ok = ok && Configure( &userx.front(),&usery.front(), userx.size(), yp0, ypn );
	return ok;
}
/*-----------------------------------------------------------------------------
 *					nxSpline2::Interpolate							2002-10-7*/
/** Interpolates the piecewise cubic spline setup by a previous call to Configure
 *	The valid interval is from the first point to just less than the last point.
 *	Based upon the code given in Numerical Recipes 
 *
 *	\param x
 *		The value at which interpolation is required
 *			
 *	\param badvalue
 *		The value to return if there are problems
 *
 *	\return
 *		the interpolated value
 **/
/*---------------------------------------------------------------------------*/

double nxSpline2::Interpolate( double x ) const
{
	const double*		begin = m_x.ArrayBasePtr();
	const double*		end   = m_x.ArrayBasePtr() + m_x.size();
	const double*		khi;
	const double*		klo;
	intptr_t			ilo;
	intptr_t			ihi;
	bool				ok;
	double				y = 0;
	double		h,a,b;
	
	ok = (m_x.size() > 1) && (x >= *begin) && ( x <= m_x.back());
	if (ok)
	{
		khi = std::upper_bound( begin, end, x );			//find  location where  (begin < x < end )
		ok = ((khi >= begin) && (khi < end));				// see if we are out of range
		if (ok)
		{
			if (khi == begin) khi++;
			klo = khi-1;
			ilo = (klo-begin);
			ihi = ilo+1;
			if ( (ilo >= (intptr_t)m_startsplineindex) && (ihi < (intptr_t)m_endsplineindex))
			{
				h   =  *khi - *klo;
 				a   = (*khi - x)/h;
				b   = (x - *klo)/h;

				y   =   a*m_y[ilo]
					  + b*m_y[ihi]
					  +( (a*a*a-a)*m_y2[ilo] + (b*b*b-b)*m_y2[ihi] )*(h*h)/6.0;
			}
			else
			{
				y = (ilo < (intptr_t)m_startsplineindex) ? m_y.front() : m_y.back();
			}
		}
	}
	if (!ok)
	{
		y = m_padvalue;
	}
	return y;
}
