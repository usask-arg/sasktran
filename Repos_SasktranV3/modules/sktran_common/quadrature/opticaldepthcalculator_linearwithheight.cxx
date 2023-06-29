#include "../sktran_common.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalDepthCalculator_LinearWithHeight::SetEndPoints		 2014- 4- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_OpticalDepthCalculator_LinearWithHeight::ConfigureQuadratureCoefficients( const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint )
{
	HELIODETIC_UNITVECTOR	lookvector;
	double					costheta0;
	double					costheta1;
	double					r0;
	double					t0;
	double					r1;
	double					t1;
	double					rt;
	double					dr;


	r0 = startpoint.Radius();							// Get the radius of the first point
	r1 = endpoint.Radius();								// Get the radius of the second point
	dr = r1 - r0;										// Get the change in radius
	if (dr >= 0 )										// get the unit look vector in the upward direction
	{													// If the last point is above the first
		lookvector = (endpoint.Vector() - startpoint.Vector()).UnitVector();
	}
	else								// or if first point is above the first
	{
		lookvector = (startpoint.Vector() - endpoint.Vector()).UnitVector();
	}
	costheta0 = startpoint.CosZenithAngle( lookvector);
	costheta1 = endpoint.CosZenithAngle( lookvector );
	t0 = r0*costheta0;									// Distance of start point from straight line tangent point
	t1 = r1*costheta1;									// Distance of start point from straight line tangent point
	rt = r0*sqrt(1.0 - costheta0*costheta0);			// Radius of straight line tangent point = r0.sin(theta0)
	ConfigureQuadratureCoefficients( r0, r1, t0, t1, rt );
	return true;
}
			

/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalDepthCalculator_LinearWithHeight::ConfigureQuadratureCoefficients		 2014- 12- 4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_OpticalDepthCalculator_LinearWithHeight::ConfigureQuadratureCoefficients(double r0, double r1, double t0, double t1, double rt )
{
	m_r0 = r0;
	m_r1 = r1;
	m_dr = r1 - r0;											// Get the change in radius
	if( fabs(m_dr) <= 0.001 )								// If the change in radius is less than 1 millimeter
	{														// Then we probably have a zero path length integral and we are close to tangential
		m_k0   = 100.0*fabs(t1-t0);							// so the quadrature coefficients are simply set to limiting values, for tangential conditions (t1-t0) approximates the path length in the shell
		m_k1   = 0.0;
		m_dt1  = 0.0;
		m_dt2  = 0.0;
		m_dr   = 0.001;
	}
	else
	{
		if ( t1 >= t0 )										// Case 1, Upward rays, or very short cell segments at the tangent point
		{
			m_dt1 = t1 - t0;
			if (std::abs(rt) < 10)
			{
				m_dt2 = 0.5*((r1*t1 - r0 * t0));
			}
			else
			{
				m_dt2 = 0.5*((r1*t1 - r0 * t0) + rt * rt*log((r1 + t1) / (r0 + t0)));
			}
		}
		else 												// case 2. downward ray, passes from upper shell to lower shell
		{
			m_dt1 = t0 - t1;
			if (std::abs(rt) < 10)
			{
				m_dt2 = 0.5*((r0*t0 - r1 * t1));
			}
			else
			{
				m_dt2 = 0.5*((r0*t0 - r1 * t1) + rt * rt*log((r0 + t0) / (r1 + t1)));
			}
		}
		m_k0 =  100.0/m_dr*( r1*m_dt1 - m_dt2);
		m_k1 = -100.0/m_dr*( r0*m_dt1 - m_dt2);
	}
	return true;
}

double SKTRAN_OpticalDepthCalculator_LinearWithHeight::OpticalDepthFromStartToEnd( double sigma0_percm, double sigma1_percm ) const
{  
	double	od;
	
	od = m_k0*sigma0_percm + m_k1*sigma1_percm;			// OD should equal this
	NXASSERT(( od >= 0.0));
	return od;
}


