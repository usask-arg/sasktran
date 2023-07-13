//#pragma once

//#include "sktran_common_internals.h"

//#include <boost/numeric/ublas/vector.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/triangular.hpp>
//#include <boost/numeric/ublas/lu.hpp>
//#include <boost/numeric/ublas/io.hpp>

//using namespace std::numeric::ublas;
//using namespace std;


/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_Profile		 2016- 7- 12*/
/** 
 *	@ingroup refindex
 * This class calculates the index of refraction as a function of altitude, wavelength, temperature, and pressure.
 * That is, it calculates n vs. h given the current optical state.  It makes use of the std::numeric::ublas
 * linear algebra library to create a cubic spline interpolation matrix for speedy and accurate interpolation
 * of the values as a function of altitude when implemented in the ray-tracing class SKTRAN_RayBaseCurvedGeometry.
 *
 * The index of refraction itself is calculated by the worker optical class skRTRefractiveIndex_MoistAir, a class
 * that makes use of Ciddor & Edlen's equation for calculting the refractivity of air.  As of yet I haven't been
 * able to incorporate humidity into the working model, but the option remains open to do so.
 *
 * Lorne Jensen, Aug 2012
 *
**/
/*---------------------------------------------------------------------------*/

class skRTRefractiveIndex_Profile
{
	protected:
		skRTRefractiveIndex_MoistAir*			m_calculator;		//!< Index of refraction for moist air - performs the actual calculations for this class
		std::vector<double>						m_refractiveindex;  //!< Array that actually holds the values of n vs. h
		std::vector<double>						m_pressure;			//!< The total pressure profile
		std::vector<double>						m_temperature;		//!< The temperature profile
		GEODETIC_INSTANT						m_referencepoint;	//!< The reference point for the optical properties call
		std::vector<double>						m_heights;			//!< The heights of the midpoints of each cell
		std::vector<double>						m_M;				//!< Spline interpolation matrix

	protected:
		void									UpdateRefractiveIndex( SKTRAN_AtmosphericOpticalState_V21 *opticalstate, double wavenum );
		bool									InitializeCubicSplineInterpolation( );

//		bool									MatrixInverse( matrix<double>& input, matrix<double>& inverse );   // This has become a static function within the module to hide the matrix and ublas interfaces

	public:
												skRTRefractiveIndex_Profile ();
		virtual								   ~skRTRefractiveIndex_Profile ();
		bool									CalculateProfile( SKTRAN_AtmosphericOpticalState_V21 *opticalstate,
																  const SKTRAN_GridDefRayTracingShells_V21 *raytracingspecs,
																  double wavelen_nm,
																  GEODETIC_INSTANT referencepoint );
		double									ExponentialLinearInterp(double altitude);
		size_t									NumPoints	() const { return m_heights.size(); }
		std::vector<double>						GetProfile	() const { return m_refractiveindex; }
		std::vector<double>						GetInterpolationMatrix() const { std::vector<double> copy; copy = m_M; return copy; }
		std::vector<double>						GetHeightProfile () const { return m_heights; }
		double									At( size_t idx ) const { if ( idx < 0 || idx > m_heights.size() - 1 ) return -1.0; else return m_refractiveindex[idx]; }
		double									HeightAt( size_t idx ) const { if ( idx < 0 || idx > m_heights.size() - 1 ) return -1.0; else return m_heights[idx]; }
};
