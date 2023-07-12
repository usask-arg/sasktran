#include <skopticalproperties21.h>

/*-----------------------------------------------------------------------------
 *					skEmission_Thermal::skEmission_Thermal			2017-8-23*/
/** **/
/*---------------------------------------------------------------------------*/

skEmission_Thermal::skEmission_Thermal()
{
	m_atmospheric_state = nullptr;
	m_groundemissivity = 1.0;
	m_temperatureK = 0.0;
	m_absorptionperm = 0.0;
	m_isground = false;
	m_absorptionisset = false;
}


/*---------------------------------------------------------------------------
 *             skEmission_Thermal::skEmission_Thermal             2020-08-20 */
/** **/
/*---------------------------------------------------------------------------*/

skEmission_Thermal::skEmission_Thermal(double groundemissivity)
{
	m_atmospheric_state = nullptr;
	m_groundemissivity = groundemissivity;
	m_temperatureK = 0.0;
	m_absorptionperm = 0.0;
	m_isground = false;
	m_absorptionisset = false;
}


/*---------------------------------------------------------------------------
 *            skEmission_Thermal::~skEmission_Thermal             2020-08-20 */
/** **/
/*---------------------------------------------------------------------------*/
skEmission_Thermal::~skEmission_Thermal()
{
	if (m_atmospheric_state != nullptr ) m_atmospheric_state->Release();
}


/*---------------------------------------------------------------------------
 *            skEmission_Thermal::SetAtmosphericState             2020-08-20 */
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_Thermal::SetAtmosphericState( skClimatology* neutralatmosphere ) 
{
	if (neutralatmosphere != nullptr) neutralatmosphere->AddRef();
	if (m_atmospheric_state != nullptr) m_atmospheric_state->Release();
	m_atmospheric_state = neutralatmosphere;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skEmission_Thermal::SetAbsorptionPerM			2018-1-23*/
/** **/
/*---------------------------------------------------------------------------*/

void skEmission_Thermal::SetAbsorptionPerM( double value )
{
	m_absorptionperm = value;
	m_absorptionisset = true;
}



/*---------------------------------------------------------------------------
 *                skEmission_Thermal::UpdateCache                 2020-08-20 */
/** **/
/*---------------------------------------------------------------------------*/

bool  skEmission_Thermal::UpdateCache( const GEODETIC_INSTANT& pt)
{
	bool ok = (m_atmospheric_state != nullptr);

	ok = ok && m_atmospheric_state->UpdateCache( pt );
	if (!ok)
	{
		if ( m_atmospheric_state == nullptr) nxLog::Record(NXLOG_WARNING,"skEmission_Thermal::UpdateLocation, Error updating temperature climatology cache as the object is not set. Call SetAtmosphericState with a valid climatology");
		else                                 nxLog::Record(NXLOG_WARNING,"skEmission_Thermal::UpdateLocation, Error updating temperature climatology cache from the climatology provided");
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *               skEmission_Thermal::UpdateLocation               2020-08-20 */
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_Thermal::UpdateLocation(const GEODETIC_INSTANT& pt, bool isground)
{
	bool ok = (m_atmospheric_state != nullptr);				// should not need to update the climatology cache since it is already checked in 'SKTRAN_AtmosphericEmission::CalculateEmissions()'
	bool updatecache = false;
	double temperature;

	// Get the temperature at current location
	ok = ok && m_atmospheric_state->GetParameter(SKCLIMATOLOGY_TEMPERATURE_K, pt, &temperature, updatecache);
	m_temperatureK = temperature;
	m_isground = isground;
	if (!ok)
	{
		if ( m_atmospheric_state == nullptr) nxLog::Record(NXLOG_WARNING,"skEmission_Thermal::UpdateLocation, Error fetching temperature as the climatology is not set. Call SetAtmosphericState with a valid climatology");
		else                                 nxLog::Record(NXLOG_WARNING,"skEmission_Thermal::UpdateLocation, Error fetching temperature from the climatology provided");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skEmission_Thermal::IsotropicEmission			2017-8-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_Thermal::IsotropicEmission(double wavenumber, double *isotropicradiance)
{
	bool ok = true;
	double wavelen_nm = 1.0e7 / wavenumber;
	double radiance;

	ok = ok && m_temperatureK >= 0;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skEmission_Thermal::IsotropicEmission, You must successfully call SetAtmosphericState before calling IsotropicEmission");
	}

	radiance = PlanckBlackbody(wavelen_nm, m_temperatureK);
	if (m_isground) {
		radiance = m_groundemissivity * radiance;
	}
	else {
		if (!m_absorptionisset) {
			nxLog::Record(NXLOG_WARNING, "skEmission_Thermal::IsotropicEmission, You must successfully call SetAbsorption before calling IsotropicEmission");
		}
		radiance = m_absorptionperm * radiance;
		m_absorptionisset = false;
	}

	ok = ok && radiance >= 0.0;

	*isotropicradiance = radiance;

	return ok;
}


/*-----------------------------------------------------------------------------
 *					skEmission_Thermal::PlanckBlackbody				2017-8-23*/
/** **/
/*---------------------------------------------------------------------------*/

double skEmission_Thermal::PlanckBlackbody(double wavelen_nm, double temperature_K)
{
	const double c = 299792458.0;			// speed of light       [m/s]
	const double h = 6.6260693e-34;			// Planck's constant    [(m^2 kg) / (s)]
	const double k = 1.3806485e-23;			// Boltzmann's constant [(m^2 kg) / (s^2 K)]
	double       wavelen_m = wavelen_nm / 1.0e9;
	double	     blackbodyradiance;

	// Compute radiance using Planck's Law. Division by 1e13 converts from units of
	// [photons / (m^2 m sr)] to units of [photons / (cm^2 nm sr)].
	blackbodyradiance = (2 * c / pow(wavelen_m, 4)) / (exp(h * c / wavelen_m / k / temperature_K) - 1) / 1e13;

	return blackbodyradiance;
}


