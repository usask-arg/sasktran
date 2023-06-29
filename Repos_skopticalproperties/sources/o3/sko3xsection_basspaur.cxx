#include <skopticalproperties21.h>

#include "bass_paur/bass_paur_o3_quadratic.hpp"

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_O3_OSIRISRes::skOpticalProperties_O3_OSIRISRes		2003-11-28
 *-------------------------------------------------------------------------*/

skOpticalProperties_O3_BassPaurQuadratic::skOpticalProperties_O3_BassPaurQuadratic()
{
	m_isdirty     = true;
	m_backgroundatmosphere = nullptr;
	m_temperature = -99999.0;
	SetInstrumentPSF_FWHM( 0.025 );
	SetInstrumentPointSpacing( 0.05  );		// Spacing of data points in nm
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_BassPaurQuadratic::~skOpticalProperties_O3_BassPaurQuadratic		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_O3_BassPaurQuadratic::~skOpticalProperties_O3_BassPaurQuadratic()
{
	if ( m_backgroundatmosphere != nullptr) m_backgroundatmosphere->Release();
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::SetAtmosphericState		2008-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O3_BassPaurQuadratic::SetAtmosphericState( skClimatology* neutralatmosphere )
{
	if ( neutralatmosphere != nullptr) neutralatmosphere->AddRef();
	if ( m_backgroundatmosphere != nullptr) m_backgroundatmosphere->Release();
	m_backgroundatmosphere = neutralatmosphere;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::SetAtmosphericState		2008-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O3_BassPaurQuadratic::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )
{
	double		kelvin = 0.0;
	bool		ok;

	ok =  ( m_backgroundatmosphere != nullptr);
	ok  = ok && m_backgroundatmosphere->GetParameter( SKCLIMATOLOGY_TEMPERATURE_K, pt, &kelvin, false );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_BassPaurQuadratic::SetLocation, Error fetching temperature from background atmosphere. Have you called SetAtmosphericState properly?");
	}
	else
	{
		Set_Temperature(kelvin);
	}
	if (crosssectionschanged != NULL) *crosssectionschanged = m_isdirty;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_BassPaurQuadratic::Set_Temperature		2009-11-5*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O3_BassPaurQuadratic::Set_Temperature( double kelvin)
{
	m_isdirty = m_isdirty || (m_temperature != kelvin);
	m_temperature = kelvin;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_BassPaurQuadratic::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_O3_BassPaurQuadratic::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_O3_BassPaurQuadratic*	clone;
	bool						ok;

	clone = new skOpticalProperties_O3_BassPaurQuadratic;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_BassPaurQuadratic::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->DeepCopy(*this);
		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_BassPaurQuadratic::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}
*/


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_BassPaurQuadratic::DeepCopy		2009-11-5*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O3_BassPaurQuadratic::DeepCopy( const skOpticalProperties_O3_BassPaurQuadratic& other )
{
	bool	ok;

	ok            =       skOpticalProperties::DeepCopy(other);
	ok            = ok && skWavelengthToPSF_TableConstant::DeepCopy(other);
	m_isdirty     = other.m_isdirty;
	m_temperature = other.m_temperature;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					basspaur_entry_lessthan2		2009-11-5*/
/** **/
/*---------------------------------------------------------------------------*/

static bool basspaur_entry_lessthan2( const skOpticalProperties_O3_BassPaurQuadratic::sk_dummy_basspaur03xsectentry& a, const skOpticalProperties_O3_BassPaurQuadratic::sk_dummy_basspaur03xsectentry& b)
{
	return (a.nm < b.nm);
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_BassPaurQuadratic::BassPaurCrossSection		2009-11-5*/
/** **/
/*---------------------------------------------------------------------------*/

double skOpticalProperties_O3_BassPaurQuadratic::BassPaurCrossSection( const sk_dummy_basspaur03xsectentry* entry ) const
{
	double	T;
	double	sigma;

	T     = m_temperature - 273.0;								// Convert Kelvin to Celsius. Note the approx value of -273.0 is used on purpose to get agreement with other people.
	if (T < -70) T = -70;										// Keep the temperature bounded to the
	if (T > +25) T = +25;										// Range specified by Bass Paur as teh fit goes goofy outside this range.
	sigma = (entry->c0 + (entry->c1 + entry->c2*T)*T)*1.0E-20;
	return sigma;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_BassPaurQuadratic::CalculateCrossSections		2009-11-5*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O3_BassPaurQuadratic::CalculateCrossSections( double wavenum, double* absxs, double* extxs, double* scattxs ) const
{
	bool									ok;
	const sk_dummy_basspaur03xsectentry*	first;
	const sk_dummy_basspaur03xsectentry*	last;
	const sk_dummy_basspaur03xsectentry*	upperentry;
	const sk_dummy_basspaur03xsectentry*	lowerentry;
	sk_dummy_basspaur03xsectentry			dummy;

	double									y1,y2,y;
	double									t1,t2,dt;
	double									nm = 1.0E7/wavenum;						// Get the nanometer wavelength
//	double									t  = Temperature();
	size_t									numentries;

	ok = (m_temperature > 50) && (m_temperature < 1000);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_BassPaurQuadratic::CalculateCrossSections, The temperature is not yet set to be in range. Make sure you call SetAtmosphericState at least once. Forcing temperature to 273.00 K");
//		m_temperature = 273.00;		// Commenting this line out makes it able to run as const, which makes it thread safe.
	}
	y          = 0.0;																	// Default to zero cross-section
	numentries = N_ELEMENTS(m_o3coeffs);
	first      = &(m_o3coeffs[0]);
	last       = first + numentries;
	dummy.nm   = nm;
	upperentry = std::upper_bound( first, last, dummy, basspaur_entry_lessthan2 );			// Find the first entry greater than our wavelength
	ok = (upperentry > first) && (upperentry < last);									// Make sure our wavelength is properly bounded by the BAss-Paur cross-sections
	if (!ok)
	{																	// then
		nxLog::Record( NXLOG_WARNING, "skOpticalProperties_O3_BassPaurQuadratic::CalculateCrossSections, The specified wavelength %e is outside the range of the BAss-Paur cross-sections (%8.4f to %8.4f)", (double)nm, (double)m_o3coeffs[0].nm, (double)m_o3coeffs[numentries-1].nm);
	}
	else
	{
		lowerentry = upperentry - 1;
		y1 = BassPaurCrossSection( lowerentry );
		y2 = BassPaurCrossSection (upperentry );
		t1 = lowerentry->nm;
		t2 = upperentry->nm;
		dt = t2-t1;
		if (dt == 0) y = y1;
		else         y = y1 + (y2-y1)*(nm-t1)/dt;
	}
	*absxs   = y;
	*extxs   = y;
	*scattxs = 0.0;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_BassPaurQuadratic::CalculateCrossSections		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O3_BassPaurQuadratic::CalculateCrossSections( double wavenum, double* absxs, double* extxs, double* scattxs )
{
	return CalculateCrossSections( wavenum, absxs, extxs, scattxs );
}



