/*-----------------------------------------------------------------------------
 *					RefractiveIndexDryAirSTP		 2015- 1- 6*/
/** Calculates the refractive index of air at STP  used mostlky for conversion
 *	of wavelength/wavenumber from vacuum to air at STP.
 **/
/*---------------------------------------------------------------------------*/

class RefractiveIndexDryAirSTP
{
	public:
		static double					RefractivityAtSTP					( double wavenum_cm1_in_vacuum );
		static double					AirWavenumberToVacuum				( double wavenum_cm1_inair );
		static double					AirWavelengthToVacuum				( double wavelen_nm_inair );
		static double					VacuumWavenumberToAir				( double wavenum_cm1_invacuum );
		static double					VacuumWavelengthToAir				( double wavelen_nm_invacuum );
};
