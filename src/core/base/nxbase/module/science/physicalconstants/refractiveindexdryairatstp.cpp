#include "nxbase_core.h"
#include "nxbase_math.h"



/*-----------------------------------------------------------------------------
 *					VacuumToAir		 2015- 1- 6*/
/** Use equation 1 from Ciddor 1996  Applied Optics , 35, 1566 to calculate refractivity
 *	of air at STP.
 *
 *  refractivity of air at STP = 10^8( n_{as}-1) = k_1/(k_0-\sigma^2) + k_3/( k_2-\sigma^2)
 *
 **/
/*---------------------------------------------------------------------------*/

double RefractiveIndexDryAirSTP::RefractivityAtSTP( double wavenum_cm1_invacuum )
{
	double			refractivity = 0.0;
	double			sigma;
	double			sigma2;
	const double	k0 = 238.0185;
	const double	k1 = 5792105.0E-08;
	const double    k2 = 57.362;
	const double    k3 = 167917.0E-08;

	if (wavenum_cm1_invacuum <= 50000.0)						// This formula only applies to wavelengths greater than 200 nm
	{															// then
		sigma  = wavenum_cm1_invacuum*1.0e-4;					// Convert wavenumber per cm to wavenumer per micrometer
		sigma2 = sigma*sigma;									// get wavenumber squared
		refractivity = k1/(k0 - sigma2) + k3/(k2 - sigma2);		// Compute conversion factorsupplied by Ciddor
	}															// and we are donr
	return refractivity;
}


/*-----------------------------------------------------------------------------
 *					RefractiveIndexDryAirSTP::VacuumWavenumToAir		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double RefractiveIndexDryAirSTP::VacuumWavenumberToAir( double wavenum_cm1_invacuum )
{
	double factor;

	factor = 1.0 +  RefractivityAtSTP(wavenum_cm1_invacuum) ;
	wavenum_cm1_invacuum *= factor;																	// Apply the conversion
	return wavenum_cm1_invacuum;
}



/*-----------------------------------------------------------------------------
 *					RefractiveIndexDryAirSTP::VacuumWavelengthToAir		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double RefractiveIndexDryAirSTP::VacuumWavelengthToAir( double wavelen_nm_invacuum )
{
	double factor;

	factor = 1.0 +  RefractivityAtSTP( 1.0E7/wavelen_nm_invacuum) ;
	wavelen_nm_invacuum /= factor;																	// Apply the conversion
	return wavelen_nm_invacuum;
}

/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndexDryAirSTP::AirWavenumToVacuum		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double RefractiveIndexDryAirSTP::AirWavenumberToVacuum( double wavenum_cm1_inair )
{
	double factor;
	double	wavenumapprox;

	factor = 1.0 +  RefractivityAtSTP(wavenum_cm1_inair) ;						// Do an approximate conversion
	wavenumapprox = wavenum_cm1_inair/factor;									// Where we assume air is the vacuum wavenumber
	factor = 1.0 +  RefractivityAtSTP(wavenumapprox) ;							// and then do it again using the new value
	wavenum_cm1_inair /= factor;												// and apply the conversion, should be pretty close
	return wavenum_cm1_inair;
}

/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndexDryAirSTP::AirWavenumToVacuum		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double RefractiveIndexDryAirSTP::AirWavelengthToVacuum( double wavelength_cm1_inair )
{
	return 1.0E7/AirWavenumberToVacuum( 1.0E7/wavelength_cm1_inair);
}

