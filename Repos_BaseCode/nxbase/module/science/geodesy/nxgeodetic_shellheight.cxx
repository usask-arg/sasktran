#include "nxbase_geodesy.h"


/*---------------------------------------------------------------------------
 *'					zbrent                                          2002-1-15
 *	Finds the zero crossing of function func(x).
 *-------------------------------------------------------------------------*/

template <class T>
double zbrent( T func, double x1, double x2, double tol, HRESULT* status)
{
	const	int    ITMAX = 500;
	const	double EPS   = 3.0e-8;
	int		iter;
	double	a=x1;
	double	b=x2;
	double	c = 0;
	double  d = 0;
	double  e = 0;
	double  min1,min2;
	double	fa= func(a);
	double	fb= func(b);
	double	fc,p,q,r,s,tol1,xm;

	*status = S_OK;
	if (fb*fa > 0.0) *status = E_FAIL;
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if (fb*fc > 0.0) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)  q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
		fb= func(b);
	}
	*status = E_FAIL;
	return -9999.00; //nrerror("Maximum number of iterations exceeded in ZBRENT");
}


/*---------------------------------------------------------------------------
 *'					class HeightOffsetEvaluator                     2002-1-15
 *	A class used whilk locating a shell heigth location.  The templated zbrent
 *	function calls the operator() to evaluate the offset in height of a
 *	certain point x
 *-------------------------------------------------------------------------*/

class HeightOffsetEvaluator
{
	nxGeodetic*	m_geoid;					// The geoid used for calculations
	nxVector	m_tangentpoint;				// The location of the tangent point of this look vector
	nxVector	m_look;						// The look direction as a unit vector (from the satellite).
	double		m_H;

	public:
				HeightOffsetEvaluator( nxGeodetic* geoid, const nxVector& tanpoint, const nxVector& look, double H );
		nxBOOL	FindCrossingPoint( double lmin, double lmax, nxVector* entrypoint );
		double	operator() ( double x );
};

/*---------------------------------------------------------------------------
 *'					HeightOffsetEvaluator::HeightOffsetEvaluator    2002-1-15
 *-------------------------------------------------------------------------*/

HeightOffsetEvaluator::HeightOffsetEvaluator( nxGeodetic* geoid, const nxVector& tanpoint, const nxVector& look, double H )
{
	m_geoid = geoid;
	m_tangentpoint = tanpoint;
	m_look         = look;
	m_H            = H;
}

/*---------------------------------------------------------------------------
 *'					HeightOffsetEvaluator::operator()               2002-1-15
 *-------------------------------------------------------------------------*/

double HeightOffsetEvaluator::operator()(double x)
{
	nxVector	location;

	location = m_tangentpoint - m_look*x;
	m_geoid->FromGeocentricVector( location );
//	printf ("Location  = %f, height = %f\n", (double)x, (double)m_geoid->Height() );
	return (m_geoid->Height()-m_H);
}


/*---------------------------------------------------------------------------
 *'					FindCrossingPoint                               2002-1-15
 *-------------------------------------------------------------------------*/

nxBOOL HeightOffsetEvaluator::FindCrossingPoint( double lmin, double lmax, nxVector* entrypoint )
{
	double					delta = 0.05*lmin;
	double					answer;
	HRESULT					status;
	HeightOffsetEvaluator&	evaluator = *this;

	answer = evaluator( lmin );									// Get the height offset closest to tangent point
	while (answer >= 0.0 )										// if
	{
		lmin  -= delta;
		answer = evaluator( lmin );
	}

	answer = evaluator( lmax );
	while (answer <= 0.0 )
	{
		lmax  += 0.05*lmax;
		answer = evaluator( lmax );
	}
	answer = zbrent( evaluator, lmin, lmax, 0.1, &status );		/* An estimate to the min location*/
	if (status == S_OK)
	{
		evaluator(answer);
		*entrypoint = m_geoid->Location();
	}
	else
	{
		entrypoint->SetCoords(0,0,0);
	}
	return (status== S_OK);
}

/*---------------------------------------------------------------------------
 *'					nxGeodetic::GetShellHeightLocation              2002-1-15
 *	Calculate the entrance and exit points of a line of sight with geodetic
 *	height  H.
 *
 *	The entrance and exit points are approximately symmetric about the tangent
 *	point.  The entrance point is the point located from the tangent point
 *	in the opposite direction to the look direction.
 *
 *	The exit point is the point from the tangent point in the direction of the
 *	look vector.
 *-------------------------------------------------------------------------*/

nxBOOL nxGeodetic::GetShellHeightLocation(  double H, const nxVector& observerposition, const nxVector& look, nxVector* entrypoint, nxVector* exitpoint,  nxVector* lostangentpoint, double tangenth)
{
	nxVector	tangentpoint;				// The location of the tangent point of this look vector
	double		Ht;
	double		Rmin;
	double		Rmax;
	nxBOOL		ok;
	nxBOOL		ok1,ok2;
	double		lmin,lmax;

	if (!look.IsValid() || look.IsZero() )
	{
		nxLog::Record(NXLOG_WARNING,"nxGeodetic::GetShellHeightLocation, the look vector is undefined");
	}

	if (lostangentpoint != NULL)
	{
		tangentpoint = *lostangentpoint;
		Ht           = tangenth;
	}
	else
	{
		FromTangentPointLocation( observerposition, look );					// Get tthe tangent point
		tangentpoint = Location();											// as this is in between the exit and entrance
		Ht = Height();														// get the tangent height
	}
	ok = (H > Ht);																// and make sure their is a valid solution
	if (ok)																		// if there is then
	{
		if (IsPureSphere())															// If we have a pure sphere
		{																			// then the solution
			lmax = sqrt( (2*m_ReA+H+Ht)*(H-Ht) );									// is considerably simplified
			nxVector	d = look.UnitVector();
			*entrypoint  = tangentpoint - lmax*d;
			if (exitpoint != NULL)
			{
				*exitpoint  = tangentpoint + lmax*d;
			}
		}
		else
		{
			HeightOffsetEvaluator	evaluator( this, tangentpoint, look, H);		// Create the height offset evaluator object

			Rmin = SemiMinorAxis();													// Get the Semi Minor axis
			Rmax = SemiMajorAxis();													// Get the Semi major axis
			lmin = 0.9*sqrt( (H - Ht)*(H + Ht +2*Rmin ) );							// Minimum distance from tangent point to shell height
			lmax = 1.1*sqrt( (H - Ht)*(H + Ht +2*Rmax ) );							// Maximum distance from tangent point to intersection

			ok1 = evaluator.FindCrossingPoint(  lmin, lmax,  entrypoint );			// Find the shell from tangent point towards observer
			ok2 = (exitpoint == NULL);
			if (!ok2)
			{
				ok2 = evaluator.FindCrossingPoint( -lmin, -lmax, exitpoint  );		// Find the shell from tangent point away from observer
			}
			ok = ok1 && ok2;
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"nxGeodetic::GetShellHeightLocation, Error retrieving entry and exit points for for height (%f)", (double)H );
			}
		}
	}
	if (!ok)
	{
		entrypoint->SetInvalid();
		if (exitpoint != NULL) exitpoint->SetInvalid();
	}
	return ok;
}

