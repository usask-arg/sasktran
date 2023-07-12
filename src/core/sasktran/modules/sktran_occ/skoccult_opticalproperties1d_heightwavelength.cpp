#include "skoccultation.h"

/*-----------------------------------------------------------------------------
 *			SKOCCULT_OpticalProperties1D_HeightWavelength::SKOCCULT_OpticalProperties1D_HeightWavelength		 2014- 4- 24*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_OpticalProperties1D_HeightWavelength::SKOCCULT_OpticalProperties1D_HeightWavelength()
{
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_OpticalProperties1D_HeightWavelength::~SKOCCULT_OpticalProperties1D_HeightWavelength		 2014- 4- 24*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_OpticalProperties1D_HeightWavelength::~SKOCCULT_OpticalProperties1D_HeightWavelength()
{
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_OpticalProperties1D_HeightWavelength::ConfigureOptical		 2014- 4- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_OpticalProperties1D_HeightWavelength::ConfigureOptical( const std::vector<double>& wavenumber, SKTRAN_AtmosphericOpticalState_V21& opticalstate )
{
	bool				ok;
	bool				ok1;
	size_t				numh;
	size_t				numw;
	GEODETIC_INSTANT	pt;
	std::vector<double>	absxs;
	std::vector<double>	extxs;
	std::vector<double>	scattxs;

	m_wavenumber = wavenumber;
	numh = m_heights.size();
	numw = m_wavenumber.size();
	ok = (numh > 0) && (numw > 0) &&  (m_coords != nullptr);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKOCCULT_OpticalProperties1D_HeightWavelength::ConfigureOptical, the optical properties table is not yet properly configured");
	}
	else
	{ 
		m_extinctionpercm.assign( numw, m_heights );					// Create the array of height profiles, 1 for each wavenumber

		pt.latitude  = m_coords->ReferencePtLatitude();
		pt.longitude = m_coords->ReferencePtLongitude();
		pt.mjd       = m_coords->ReferencePointMJD();
		for (size_t ih = 0; ih < numh; ih++)
		{
			pt.heightm   = m_heights.at(ih);
			ok1 =        opticalstate.SetTimeAndLocation( pt, false );
			ok1 = ok1 && opticalstate.CalculateMultiWaveCrossSections( wavenumber, &absxs, &extxs, &scattxs );
			for (size_t iw = 0; iw < numw; iw++)
			{
				m_extinctionpercm.at(iw).at(ih) = extxs.at(iw);
			}
			ok = ok && ok1;
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKOCCULT_OpticalProperties1D_HeightWavelength::ConfigureOptical, There were errors configuring the optical properties table. Clearing the table");
			m_extinctionpercm.resize(0);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_OpticalProperties1D_HeightWavelength::TotalExtinctionPerCM		 2014- 4- 24*/
/** **/
/*---------------------------------------------------------------------------*/

double SKOCCULT_OpticalProperties1D_HeightWavelength::TotalExtinctionPerCM( const HELIODETIC_POINT& point, size_t wavenumberindex  ) const
{
	return nxLinearInterpolate::LogInterpolateYatX( point.Altitude(), m_heights, m_extinctionpercm.at(wavenumberindex), nxLinearInterpolate::ENUM_TRUNCATE, 0);
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_OpticalProperties1D_HeightWavelength::ConfigureGeometry		 2014- 4- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_OpticalProperties1D_HeightWavelength::ConfigureGeometry( const SKOCCULT_Specs_Internal* specs )
{

	m_coords  = specs->CoordinateSystemObject();

	m_heights = specs->OpticalTableSpecs()->OpticalPropertiesGrid()->Altitudes();
	return (m_heights.size() > 1);
}

