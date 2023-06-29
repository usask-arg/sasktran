#include "include/sktran_hr_internals.h"

void SKTRAN_HR_WF_Extinction_Table::FillTable( const SKTRAN_HR_WF_Store& pertlist,
									   skOpticalProperties& optprop,
									   const SKTRAN_CoordinateTransform_V2& coords,
									   skClimatology& neutral,
									   const std::vector<double>& wavelen,
									   double mjd )
{
	m_table.SetSize( wavelen.size(), pertlist.StoreSize() ); 
	std::vector<double> extxs;
	std::vector<double> absxs;
	std::vector<double> scattxs;

	m_wavenumber.resize( wavelen.size() );
	extxs.resize( wavelen.size() );
	absxs.resize( wavelen.size() );
	scattxs.resize( wavelen.size() );

	for( int idx = (int) wavelen.size()-1; idx >=0 ; idx-- )	{
		m_wavenumber[idx] = 1e7 / wavelen[ wavelen.size() - idx - 1 ];
	}

	for( size_t pertidx = 0; pertidx < pertlist.StoreSize(); pertidx++ )
	{
		HELIODETIC_POINT pertpoint = pertlist.RawAccess(pertidx)->PerturbationLocation( coords );
		GEODETIC_INSTANT pertgeo = coords.PointToGeodetic( pertpoint, mjd );
		
		bool change = true;
		optprop.SetAtmosphericState( &neutral);
		optprop.SetLocation( pertgeo, &change );
		optprop.CalculateCrossSectionsArray( &m_wavenumber.front(), (int)m_wavenumber.size(), &absxs.front(), &extxs.front(), &scattxs.front() );
		for(size_t idx = 0; idx < extxs.size(); idx++ )
		{
			m_table.At( idx, pertidx ) = extxs[idx];
		}
	}
}

double SKTRAN_HR_WF_Extinction_Table::ExtinctionAt( double wavelen, size_t pertidx )
{
	size_t idx = std::find_if( std::begin(m_wavenumber), std::end(m_wavenumber), [&] (double a) {return abs(a- 1e7/wavelen) < 1e-10;}) - std::begin(m_wavenumber);

	return m_table.At(idx, pertidx);
}