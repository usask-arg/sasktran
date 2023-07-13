#include "nxbase_math.h"


/*-----------------------------------------------------------------------------
 *					nxPiecewiseLinear::nxPiecewiseLinear		2002-10-7						*/
/** Default Constructor **/
/*---------------------------------------------------------------------------*/

nxPiecewiseLinear::nxPiecewiseLinear()
{
	m_x        = NULL;
	m_y        = NULL;
	m_numcoefs = 0;
	m_outofboundsaction = nxLinearInterpolate::ENUM_TRUNCATE;
	m_missingval        = std::numeric_limits<double>::quiet_NaN();
}

/*---------------------------------------------------------------------------
 *              nxPiecewiseLinear::nxPiecewiseLinear              2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/
nxPiecewiseLinear::nxPiecewiseLinear(const nxPiecewiseLinear& other)
{
	m_x        = NULL;
	m_y        = NULL;
	m_numcoefs = 0;
	DeepCopy(other);
}


/*-----------------------------------------------------------------------------
 *					nxPiecewiseLinear::~nxPiecewiseLinear								2002-10-7*/
/** Destructor **/
/*---------------------------------------------------------------------------*/

nxPiecewiseLinear::~nxPiecewiseLinear()
{
	ReleaseResources();
}

/*---------------------------------------------------------------------------
 *'					nxPiecewiseLinear::ReleaseResources                      2002-10-7
 *-------------------------------------------------------------------------*/

void nxPiecewiseLinear::ReleaseResources()
{
	if (m_x != NULL) delete [] m_x;
	m_x        = NULL;
	m_y        = NULL;
	m_numcoefs = 0;
}

/*---------------------------------------------------------------------------
 *'					nxPiecewiseLinear::Allocate                              2002-10-7
 *-------------------------------------------------------------------------*/

nxBOOL nxPiecewiseLinear::Allocate(size_t ncoefs )
{
	nxBOOL	ok = nxTRUE;

	if (ncoefs != m_numcoefs)
	{
		ReleaseResources();
		if (ncoefs > 0)
		{
			m_x        = new double [2*ncoefs];
			m_y        = m_x  + ncoefs;
			m_numcoefs = ncoefs;
			ok         = (m_x != NULL);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "nxPiecewiseLinear::Allocate, Error allocating %d elements for the spline curve", (int)ncoefs);
				ReleaseResources();
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxPiecewiseLinear::DeepCopy		2008-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxPiecewiseLinear::DeepCopy( const nxPiecewiseLinear& other )
{
	bool		ok;
	size_t		i;

	ok = Allocate( other.m_numcoefs );
	if (ok)
	{
		for (i = 0; i < other.m_numcoefs; i++)
		{
			m_x[i]  = other.m_x[i];
			m_y[i]  = other.m_y[i];
		}
	}
	m_outofboundsaction = other.m_outofboundsaction;
	m_missingval        = other.m_missingval;

	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxPiecewiseLinear::Configure								2002-10-7*/
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
 *			
 *
 *	\return 
 *		nxTRUE if successful
 *
 **/
/*---------------------------------------------------------------------------*/

bool nxPiecewiseLinear::Configure(const double* x,  const double* y, size_t n)
{
	bool				ok;
	double				lastval = x[0];
	bool				isascending = true;
	ok = ( n > 1);
	if (ok)
	{
		ok = Allocate( (int)n );
		if (ok)
		{
			for (size_t i = 0; i < n; i++)
			{
				isascending = isascending && (lastval <= x[i]);
				lastval = x[i];
				m_x[i]  = x[i];
				m_y[i]  = y[i];
			}
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "nxPiecewiseLinear::Configure. The x values provided are not in ascending order. Thats a problem.");
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxPiecewiseLinear::Configure								2002-10-7*/
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
 *	\return 
 *		nxTRUE if successful
 *
 **/
/*---------------------------------------------------------------------------*/

bool nxPiecewiseLinear::Configure( const nx1dArray<double>& x,  const nx1dArray<double>& y)
{
	const double*	xptr =x.ArrayBasePtr();
	const double*	yptr =y.ArrayBasePtr();
	bool			ok;

	ok = x.ArrayRankSpecs()->IsContiguous() && (y.ArrayRankSpecs()->IsContiguous());
	ok = ok && (x.size() == y.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxPiecewiseLinear::Configure for nx1dArray only works for contiguous nxArrays of the same size");
	}
	else
	{
		ok = Configure( xptr,  yptr, x.size());
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxPiecewiseLinear::Configure		2013-11-1*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxPiecewiseLinear::Configure( const std::vector<double>& userx,  const std::vector<double>& usery)
{
	bool				ok;

	ok =  (userx.size() == usery.size());
	ok = ok && Configure( &userx.front(),&usery.front(), userx.size());
	return ok;
}
/*-----------------------------------------------------------------------------
 *					nxPiecewiseLinear::Interpolate							2002-10-7*/
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

double nxPiecewiseLinear::Interpolate( double x, double badvalue ) const
{
	return nxLinearInterpolate::EvaluateYatX( x, m_x, m_y, m_numcoefs, m_outofboundsaction, m_missingval  );

}
