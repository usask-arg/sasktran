/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalDepthCalculator_LinearWithHeight		 2014- 4- 24*/
/** @ingroup odintegrate 
 *	A class to calculate the optical depth between two points in the atmosphere.
 *	The class is given the location of the end points and the extinction per cm
 *	at bioth ends. It then assumes extinctions varies linearly with radius
 *	ie. simga(r) = sigmaK + r*sigmaF. This calculation is only strictly valid
 *	for atmospheres that are radially symmetric. 
 *
 *	The code calculates the total extinction between two points joined by a
 *	straight line, ie. no ray curvature calculations are made. The code 
 *	handles most situations but the user must make sure that the two end 
 *	points do not straddle the tangent point as this messes up the integral
 *	(no checks are made in this code for speed sake).
 *
 *	The class was written in a two stage process so it could be used more efficiently
 *	in the occultation engine. The optimization involves doing most of the calculation
 *	once the end points are defined in #SetEndPoints, as this only involves geometry but does require
 *	transcendental functions like log. OpticalDepths are then quickly calculated for each
 *	wavelength/wavenumber using #OpticalDepthFromStartToEnd.
 *
 *	The class stores internal geometry information and cannot be used in a
 *	multi-threaded environment where different threads make calls to SetEndPoints
 *	on the same instance. Method OpticalDepthFromStartToEnd is completely thread safe.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_OpticalDepthCalculator_LinearWithHeight
{
	private:
		double		m_k0;
		double		m_k1;
		double		m_dt1;
		double		m_dt2;
		double		m_r0;
		double		m_r1;
		double		m_dr;

	public:
		bool		ConfigureQuadratureCoefficients	( const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint );
		bool		ConfigureQuadratureCoefficients	( double r0, double r1, double t0, double t1, double rt );
		double		SigmaAtRadius					( double r, double sigma0, double sigma1   ) const { return ((m_r1-r)*sigma0 + (r-m_r0)*sigma1)/m_dr;}
		double		OpticalDepthFromStartToEnd		( double sigma0_percm, double sigma1_percm ) const;// {  return m_k0*sigma0_percm + m_k1*sigma1_percm;}
};



