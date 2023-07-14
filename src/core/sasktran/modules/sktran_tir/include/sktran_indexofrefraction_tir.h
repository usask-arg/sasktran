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

class skRTRefractiveIndex_Profile_TIR : public skRTRefractiveIndex_Profile
{

	private:
		void									UpdateRefractiveIndexTIR( SKTRAN_TIR_AtmosphericOpticalState *opticalstate, double wavenum );

//		bool									MatrixInverse( matrix<double>& input, matrix<double>& inverse );   // This has become a static function within the module to hide the matrix and ublas interfaces

	public:
												skRTRefractiveIndex_Profile_TIR ();
		virtual								   ~skRTRefractiveIndex_Profile_TIR ();
		bool									CalculateProfileTIR( SKTRAN_TIR_AtmosphericOpticalState* opticalstate,
																     const SKTRAN_GridDefRayTracingShells_V21* raytracingspecs,
																     double wavelen_nm,
																     GEODETIC_INSTANT referencepoint );
};
