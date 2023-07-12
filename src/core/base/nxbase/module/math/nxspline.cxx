#include "nxbase_math.h"


/*-----------------------------------------------------------------------------
 *					nxSpline::nxSpline		2002-10-7						*/
/** Default Constructor **/
/*---------------------------------------------------------------------------*/

nxSpline::nxSpline()
{
	m_x        = NULL;
	m_y        = NULL;
	m_y2       = NULL;
	m_numcoefs = 0;
}


/*-----------------------------------------------------------------------------
 *					nxSpline::~nxSpline								2002-10-7*/
/** Destructor **/
/*---------------------------------------------------------------------------*/

nxSpline::~nxSpline()
{
	ReleaseResources();
}

/*---------------------------------------------------------------------------
 *'					nxSpline::ReleaseResources                      2002-10-7
 *-------------------------------------------------------------------------*/

void nxSpline::ReleaseResources()
{
	if (m_x != NULL) delete [] m_x;
	m_x        = NULL;
	m_y        = NULL;
	m_y2       = NULL;
	m_u        = NULL;
	m_numcoefs = 0;
}

/*---------------------------------------------------------------------------
 *'					nxSpline::Allocate                              2002-10-7
 *-------------------------------------------------------------------------*/

nxBOOL nxSpline::Allocate(int ncoefs )
{
	nxBOOL	ok = nxTRUE;

	if (ncoefs != m_numcoefs)
	{
		ReleaseResources();
		if (ncoefs > 0)
		{
			m_x        = new double [4*ncoefs];
			m_y        = m_x  + ncoefs;
			m_y2       = m_y  + ncoefs;
			m_u        = m_y2 + ncoefs;
			m_numcoefs = ncoefs;
			ok         = (m_x != NULL);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "nxSpline::Allocate, Error allocating %d elements for the spline curve", (int)ncoefs);
				ReleaseResources();
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxSpline::DeepCopy		2008-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline::DeepCopy( const nxSpline& other )
{
	bool	ok;
	int		i;

	ok = Allocate( other.m_numcoefs );
	if (ok)
	{
		for (i = 0; i < other.m_numcoefs; i++)
		{
			m_x[i]  = other.m_x[i];
			m_y[i]  = other.m_y[i];
			m_y2[i] = other.m_y2[i];
			m_u[i]  = other.m_u [i];
		}
	}
	return ok;
}




/*-----------------------------------------------------------------------------
 *					nxSpline::Configure								2002-10-7*/
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

bool nxSpline::Configure(const double* x,  const double* y, size_t n, double yp0, double ypn)
{
	int					i,k;
	double				p,qn,sig,un;			//,*u,*vector();
//	double*				u;
	nxBOOL				ok;
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

	ok = ( n > 1);
	if (ok)
	{
		ok = Allocate( (int)n );
		if (ok)
		{
			for (i=0; i < (int)n; i++)
			{
				m_x[i]  = x[i];
				m_y[i]  = y[i];
			}

			if (yp0 >= 1.0e30)
			{
				m_y2[0] = 0.0;
				m_u [0] = 0.0;
			}
			else
			{
				m_y2[0] = -0.5;
				m_u [0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp0);
			}

			um1 = m_u[0];
			xm1 = x[0];
			ym1 = y[0];
			xi  = x[1];
			yi  = y[1];
			for (i=1; i< (int)n-1; i++)
			{
				xp1     = x[i+1];
				yp1     = y[i+1];
				xp1mxm1 = xp1-xm1;
				ximxm1  = xi-xm1;
	//			sig     = (x[i] - x[i-1])/ (x[i+1]-x[i-1]);									// THIS CAN GO
				sig     = ximxm1/xp1mxm1;
				p       = sig*m_y2[i-1]+2.0;
				m_y2[i] = (sig-1.0)/p;
				ui      = (yp1   - yi )/(xp1-xi) - (yi-ym1)/(ximxm1);
				ui      = (6.0*ui/(xp1mxm1)-sig*um1)/p;
				m_u[i]  = ui;

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
			m_y2[n-1]=(un-qn*m_u[n-2])/(qn*m_y2[n-2]+1.0);
			for (k=(int)n-2; k >=0; k--)
			{
				m_y2[k]= m_y2[k]*m_y2[k+1]+m_u[k];
			}
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxSpline::Configure								2002-10-7*/
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

bool nxSpline::Configure( const nx1dArray<double>& x,  const nx1dArray<double>& y, double yp0, double ypn)
{
	const double*	xptr =x.ArrayBasePtr();
	const double*	yptr =y.ArrayBasePtr();
	bool			ok;

	ok = x.ArrayRankSpecs()->IsContiguous() && (y.ArrayRankSpecs()->IsContiguous());
	ok = ok && (x.size() == y.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxSpline::Configure for nx1dArray only works for contiguous nxArrays of the same size");
	}
	else
	{
		ok = Configure( xptr,  yptr, x.size(), yp0, ypn);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxSpline::Configure		2013-11-1*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxSpline::Configure( const std::vector<double>& userx,  const std::vector<double>& usery, double yp0, double ypn)
{
	bool				ok;

	ok =  (userx.size() == usery.size());
	ok = ok && Configure( &userx.front(),&usery.front(), userx.size(), yp0, ypn );
	return ok;
}
/*-----------------------------------------------------------------------------
 *					nxSpline::Interpolate							2002-10-7*/
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

double nxSpline::Interpolate( double x, double badvalue ) const
{
	double*		begin = m_x;
	double*		end   = m_x + m_numcoefs;
	double*		khi;
	double*		klo;
	intptr_t	ilo;
	intptr_t		ihi;
	bool		ok;
	double		y = 0;
	double		h,a,b;
	
	ok = (m_numcoefs > 0);
	if (ok)
	{
		khi = std::upper_bound( begin, end, x );			//find  location where  (begin < x < end )
		ok = ((khi > begin) && (khi < end));				// see if we are out of range
		if (ok)
		{
			klo = khi-1;
			ilo = (klo-begin);
			ihi = ilo+1;
			h   =  *khi - *klo;
 			a   = (*khi - x)/h;
			b   = (x - *klo)/h;

			y   =   a*m_y[ilo]
				  + b*m_y[ihi]
				  +( (a*a*a-a)*m_y2[ilo] + (b*b*b-b)*m_y2[ihi] )*(h*h)/6.0;
		}
	}
	if (!ok)
	{
		y = badvalue;
	}
	return y;
}

/*-----------------------------------------------------------------------------
 *					nxSpline::Integrate								2002-10-7*/
/**	Integrates the spline from the first point to the last.  It returns the summed
 *	integral at each control point in the spline. This allows a quick integral
 *	to be performed on a set of data with pretty good precision.
 *
 *	\param integral
 *		The array returns an integral value for each control point in the
 *		spline. May be NULL in which case only the integrated sum is
 *		returned. 
 *
 *	\return 
 *		The total integrated value from the beginning to the end
 *
 **/
/*---------------------------------------------------------------------------*/

double nxSpline::Integrate( nx1dArray<double>* integral ) const
{
	double		x1,x2;
	double		dx;
	double		y1,y2;
	double		yp1,yp2;
	int			i;
	double		sum;

	sum = 0.0;
	if (m_numcoefs > 1)
	{
		if (integral != NULL)
		{
			integral->SetSize(m_numcoefs);
			integral->At(0) = 0;
		}
		x2  = m_x[0];
		y2  = m_y[0];
		yp2 = m_y2[0];
		for (i = 1; i < m_numcoefs; i++)
		{
			x1   = x2;
			x2   = m_x[i];
			y1   = y2;
			y2   = m_y[i];
			yp1  = yp2;
			yp2  = m_y2[i];
			dx   = (x2-x1);
			sum += dx*0.5*( (y1+y2) - (yp1+yp2)*(dx*dx)/6.0 );
			if (integral != NULL) integral->At(i) = sum;
		}
	}
	else
	{
		if (integral != NULL) integral->erase();
	}
	return sum;
}


