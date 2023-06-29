#include "nxbase_math.h"
#include <float.h>


/*---------------------------------------------------------------------------
 *'					nxQuaternion::nxQuaternion                     2001-11-28
 *-------------------------------------------------------------------------*/

nxQuaternion::nxQuaternion()
{
	m_q[0] = 0;
	m_q[1] = 0;
	m_q[2] = 0;
 	m_q[3] = 0;
}


/*---------------------------------------------------------------------------
 *'					nxQuaternion::nxQuaternion                     2001-11-28
 *-------------------------------------------------------------------------*/

nxQuaternion::nxQuaternion( double a, double b, double c, double d)
{
	FromScalars( a,b,c,d );
}

/*---------------------------------------------------------------------------
 *'					nxQuaternion::FromScalars                      2001-11-28
 *-------------------------------------------------------------------------*/

void nxQuaternion::FromScalars( double a, double b, double c, double d)
{
	m_q[0] = a;
	m_q[1] = b;
	m_q[2] = c;
 	m_q[3] = d;
}

/*---------------------------------------------------------------------------
 *'					nxQuaternion::FromQUATERNION                   2001-11-28
 *-------------------------------------------------------------------------*/

void nxQuaternion::FromQUATERNION( const QUATERNION& Q )
{
	m_q[0] = Q[0];
	m_q[1] = Q[1];
	m_q[2] = Q[2];
 	m_q[3] = Q[3];
}

//---------------------------------------------------------------------------
//						nxQuaternion::nxQuaternion
//---------------------------------------------------------------------------

nxQuaternion::nxQuaternion( double a, const nxVector& qv )
{
	FromScalarAndVector( a, qv );
}

//---------------------------------------------------------------------------
//						nxQuaternion::Conjugate
//---------------------------------------------------------------------------

nxQuaternion nxQuaternion::Conjugate()
{
	return nxQuaternion( m_q[0], -m_q[1], -m_q[2], -m_q[3] );
}

/*---------------------------------------------------------------------------
 *'					nxQuaternion::Norm                             2001-11-28
 *-------------------------------------------------------------------------*/

double nxQuaternion::Norm() const
{
	return sqrt (m_q[0]*m_q[0] + m_q[1]*m_q[1]  + m_q[2]*m_q[2]	 + m_q[3]*m_q[3]  );
}

/*---------------------------------------------------------------------------
 *'					nxQuaternion::Norm                             2001-11-28
 *-------------------------------------------------------------------------*/

void nxQuaternion::Normalize()
{
	double norm = Norm();
	double factor;

	if (norm > DBL_MIN)
	{
		factor = 1.0/norm;
		m_q[0] *= factor;
		m_q[1] *= factor;
		m_q[2] *= factor;
		m_q[3] *= factor;
	}
}

/*---------------------------------------------------------------------------
 *'					nxQuaternion::DotProduct                       2001-11-28
 *-------------------------------------------------------------------------*/

double nxQuaternion::DotProduct( const nxQuaternion& other ) const
{
	return m_q[1]*other.m_q[1] + m_q[2]*other.m_q[2] + m_q[3]*other.m_q[3] + m_q[0]*other.m_q[0];
}


/*<--------------------------------------------------------------------------
 *'					nxQuaternion::UnitSlerp
 *	Provides spherical linear interpolation between this quaternion and another.
 *	The interpolation is controlled by the parameter alpha which is assumed
 *	to be in the range [0,1].
 *
 *	If alpha is 0 then it should return "this" quaternion.  If alpha is 1.0 it
 *	will return the "endquat" quaternion.  For alpha between 0 and 1 it will return
 *	a linear interpolation of the two.
 *
 *	This code was lifted from:
 *  http://www.acm.org/pubs/tog/GraphicsGems/gemsiii/quatspin.c
 *	so it has reasonable heritage.
 *
 *'	PARAMETERS:
 *`		const nxQuaternion& endquat
 *			The endpoint quaternion (assumed to be normalized to unity)
 *
 *`		double alpha
 *			The interpolation parameter [0..1].  ) returns "this" quaternion
 *			1 returns endquat quaternion
 *
 *'	RETURNS:
 *	nxQuaternion
 *
 *'	HISTORY:
 * 2001-11-28
 *>------------------------------------------------------------------------*/

nxQuaternion nxQuaternion::UnitSlerp( const nxQuaternion& endquat, double alpha) const
{
	double			beta;							// complementary interp parameter
	double			theta;							// angle between A and B
	double			sin_t, cos_t;					// sine, cosine of theta
	double			phi;							// theta plus spins
	nxBOOL			bflip;							// use negation of B?

	cos_t = DotProduct(endquat);						// cosine theta = dot product of A and B
 	if (cos_t < 0.0)								// if B is on opposite hemisphere from A, use -B instead
	{
		cos_t = -cos_t;
		bflip = nxTRUE;
	}
	else
	{
		bflip = nxFALSE;
	}

	if ((1.0 - cos_t) < DBL_EPSILON)				   	// if B is (within precision limits) the same as A,
	{												// just linear interpolate between A and B.
		beta = 1.0 - alpha;							//
 	}
	else											// otherwise we have
	{												// the normal case
 		theta = acos(cos_t);
 		phi   = theta;
 		sin_t = sin(theta);
 		beta  = sin(theta - alpha*phi) / sin_t;
 		alpha = sin(alpha*phi) / sin_t;
 	}
	if (bflip) alpha = -alpha;

	/* interpolate */
 	double w = beta*m_q[0]  + alpha*endquat.m_q[0];
 	double x = beta*m_q[1]  + alpha*endquat.m_q[1];
  	double y = beta*m_q[2]  + alpha*endquat.m_q[2];
	double z = beta*m_q[3]  + alpha*endquat.m_q[3];

	return nxQuaternion( w, x, y, z );
}


//---------------------------------------------------------------------------
//						nxQuaternion::VectorComponent
//---------------------------------------------------------------------------

nxVector nxQuaternion::VectorComponent() const
{
	return nxVector( m_q[1], m_q[2], m_q[3] );
}

//---------------------------------------------------------------------------
//						nxQuaternion::ScalarComponent
//---------------------------------------------------------------------------

double nxQuaternion::ScalarComponent() const
{
	return m_q[0];
}

//---------------------------------------------------------------------------
//						nxQuaternion::FromScalarAndVector
//---------------------------------------------------------------------------

void nxQuaternion::FromScalarAndVector(double a, const nxVector& qv)
{
	m_q[0] = a;
	m_q[1] = qv.X();
	m_q[2] = qv.Y();
 	m_q[3] = qv.Z();
}

/*---------------------------------------------------------------------------
 *'					nxQuaternion::FromEuler                        2001-11-27
 *	Converts a 3 element array of Euler angle (x,y,z) degrees to a quaternion
 *	assuming the Euler angles are applied in the sequence zyx.
 *-------------------------------------------------------------------------*/

void nxQuaternion::FromEuler( const nxVector& eulerdegrees)
{
	nxVector t;
	double sintx;
	double sinty;
	double sintz;
	double costx;
	double costy;
	double costz;

	t     = eulerdegrees*(0.5*nxmath::ONE_DEGREE);
	sintx = sin(t.X());
	sinty = sin(t.Y());
	sintz = sin(t.Z());
	costx = cos(t.X());
	costy = cos(t.Y());
	costz = cos(t.Z());

	m_q[0]  =  sintx*sinty*sintz + costx*costy*costz;
	m_q[1]  = -costx*sinty*sintz + sintx*costy*costz;
	m_q[2]  =  costx*sinty*costz + sintx*costy*sintz;
	m_q[3]  = -sintx*sinty*costz + costx*costy*sintz;
}

//---------------------------------------------------------------------------
//						nxQuaternion::operator*
//---------------------------------------------------------------------------

nxQuaternion nxQuaternion::operator*( const nxQuaternion& other )
{
	nxQuaternion	answer = *this;
	return (answer *= other);
}


//---------------------------------------------------------------------------
//						nxQuaternion::operator*
//---------------------------------------------------------------------------

nxQuaternion& nxQuaternion::operator*=( const nxQuaternion& other )
{
	double a,b,c,d;
	double x,y,z,w;

	a = m_q[0];
	b = m_q[1];
	c = m_q[2];
	d = m_q[3];

	x = other.m_q[0];
	y = other.m_q[1];
	z = other.m_q[2];
	w = other.m_q[3];

	double tmp_00 = (d-c)*(z-w);
	double tmp_01 = (a+b)*(x+y);
	double tmp_02 = (a-b)*(z+w);
	double tmp_03 = (c+d)*(x-y);
	double tmp_04 = (d-b)*(y-z);
	double tmp_05 = (d+b)*(y+z);
	double tmp_06 = (a+c)*(x-w);
	double tmp_07 = (a-c)*(x+w);
	double tmp_08 = tmp_05 + tmp_06 + tmp_07;
	double tmp_09 = 0.5*(tmp_04 + tmp_08);
	FromScalars( tmp_00 + tmp_09 - tmp_05,  tmp_01 + tmp_09 - tmp_08, tmp_02 + tmp_09 - tmp_07, tmp_03 + tmp_09 - tmp_06 );
	return *this;
}

//---------------------------------------------------------------------------
//						nxQuaternion::RotateVector
//---------------------------------------------------------------------------

nxVector nxQuaternion::RotateVector( const nxVector& original )
{
	nxQuaternion	p( 0, original );
	return (((*this)*p)*Conjugate()).VectorComponent();
}

/*---------------------------------------------------------------------------
 *'					nxQuaternion::FromAxisToAxisRotation             2003-2-2
 *-------------------------------------------------------------------------*/

void nxQuaternion::FromAxisToAxisRotation  ( const nxVector& from, const nxVector& to )
{
    nxVector	vhalf;
//	nxVector	cross;
	nxVector	fromu;
	nxVector	axis;

	fromu = from.UnitVector();
	vhalf = (fromu + to.UnitVector()).UnitVector();		// Get the vector halfway between the two
	axis  = fromu^vhalf;								// Get the rotation axis  (sin (theta/2) factor)
	FromScalarAndVector( (fromu & vhalf), axis);		//
}

/*---------------------------------------------------------------------------
 *'					nxQuaternion::SetToIdentity                      2003-2-2
 *-------------------------------------------------------------------------*/

void nxQuaternion::SetToIdentity()
{
	m_q[0] = 1;
	m_q[1] = 0;
	m_q[2] = 0;
	m_q[3] = 0;
}


