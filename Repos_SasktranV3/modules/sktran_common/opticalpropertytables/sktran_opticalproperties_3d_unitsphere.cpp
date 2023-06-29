#include "../sktran_common.h"
#include <iostream>
#include <omp.h>


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::SKTRAN_OpticalPropertiesTable_UnitSphere		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableOpticalProperties_3D_UnitSphere::SKTRAN_TableOpticalProperties_3D_UnitSphere()
{
	m_unitsphere = NULL;
	m_scatteranglegrid = NULL;
	m_wavelengthgrid = NULL;
	m_albedo = NULL;
	m_alts = NULL;
	m_speedhelper.resize( omp_get_max_threads() );
	m_firsttime = true;
	m_forcecacheupdates = false;
	m_inonedimmode = false;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::~SKTRAN_TableOpticalProperties_3D_UnitSphere		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableOpticalProperties_3D_UnitSphere::~SKTRAN_TableOpticalProperties_3D_UnitSphere()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::ConfigureGeometry		2013-09-09*/
/**  
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::ConfigureGeometry( const SKTRAN_Specifications_Base* specs )
{
	bool ok = true;

	size_t		numloc = m_unitsphere->NumUnitVectors();

	if (numloc == 1)
	{
		m_inonedimmode = true;
	}

	size_t		numalt = m_alts->NumAltitudes();
	size_t		numscat = m_scatteranglegrid->NumAngles();
	//size_t      numwavel = m_wavelengtharray.size();
	size_t		numwavel = m_wavelengthgrid == NULL ? 1 : m_wavelengthgrid->NumWavelengths();

	if( 0 == numwavel )
	{
		numwavel = 1;
	}

	m_totalextinction.resize( numwavel );
	m_scatextinction.resize( numwavel );
	for( size_t wavelidx = 0; wavelidx < numwavel; wavelidx++ )
	{
		m_totalextinction[wavelidx].resize( numloc );
		m_scatextinction[wavelidx].resize( numloc );
		for( size_t locidx = 0; locidx < numloc; locidx++ )
		{
			m_totalextinction[wavelidx][locidx].resize( numalt );
			m_scatextinction[wavelidx][locidx].resize( numalt );
		}
	}

	m_scatprops->Allocate( numwavel*numloc*numalt*numscat );

	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::ConfigureOptical		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::ConfigureOptical( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate )
{
	bool ok = true;

	//return ConfigureOpticalMT( wavelen, opticalstate );

	opticalstate.PhaseGridHint(m_scatanglevector);

	m_wavelindex = std::find_if( std::begin(m_wavenumber), std::end(m_wavenumber), [&] (double a) {return abs(a- 1e7/wavelen) < 1e-10;}) - std::begin(m_wavenumber);

	skBRDF*					albedo;
	GEODETIC_INSTANT		geopoint;
	HELIODETIC_VECTOR		heliovector;
	nxVector				geovector;

	ok = ok && opticalstate.GetAlbedoObject( &albedo );
	if( ok )
	{
		if( NULL != m_albedo ) m_albedo->Release();
		m_albedo = albedo;
		m_albedo->AddRef();
		m_mjd = opticalstate.GetTimeAndLocation().mjd;
		m_wavelen = wavelen;
		opticalstate.SetWavelength( m_wavelen );
		
		geopoint.mjd = m_mjd;
		for( size_t locidx = 0; locidx < m_unitsphere->NumUnitVectors(); locidx++ )
		{	// unit sphere is in heliodetic coordinates, so convert to geographic
			heliovector.SetCoords( m_unitsphere->UnitVectorAt( locidx ).X(),
						   m_unitsphere->UnitVectorAt( locidx ).Y(),
						   m_unitsphere->UnitVectorAt( locidx ).Z() );

			geovector = CoordinatesPtr()->HelioVectorToGeographic( heliovector );
			geopoint.latitude  = geovector.Latitude();
			geopoint.longitude = geovector.Longitude();
			geopoint.heightm = 0.0;
			ok = ok && opticalstate.SetTimeAndLocation( geopoint, m_forcecacheupdates );
			for( size_t altidx = 0; altidx < m_alts->NumAltitudes(); altidx++ )
			{
				geopoint.heightm   = m_alts->At( altidx );
				ok = ok && opticalstate.SetTimeAndLocation( geopoint, false );
				if( 0 == m_wavelengtharray.size() )
					ok = ok && FillTablesAtIndex( altidx, locidx, opticalstate );
				else
					ok = ok && FillTablesAtIndexMultiWavel( altidx, locidx, opticalstate );
			}
		}
	}
	//DumpTotalExtinction();
	m_firsttime = false;
	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::FillTablesAtIndex( size_t altidx, size_t locidx, SKTRAN_AtmosphericOpticalState_V21& opticalstate )
{
	bool ok = true;

	double					cosangle;
	double					p11;
	skRTPhaseMatrix			pmatrix;

	bool polarized = m_scatprops->NeedsFullMueller();

	m_totalextinction[m_wavelindex][locidx][altidx] = opticalstate.ExtinctionPercm();
	m_scatextinction[m_wavelindex][locidx][altidx]  = opticalstate.ScatteringPercm();

	for( SKTRAN_GridIndex scatidx = 0; scatidx < m_scatteranglegrid->NumAngles(); scatidx++ )
	{
		cosangle = m_scatteranglegrid->At( scatidx );
		if (polarized)
		{
			ok = ok && opticalstate.VectorPhaseMatrix(cosangle, &pmatrix);
		}
		else
		{
			ok = ok && opticalstate.ScalarPhaseMatrix(std::pair<double, size_t>{cosangle, (size_t)scatidx}, p11);
			pmatrix.At(1, 1) = p11;
		}
		m_scatprops->StorePolarizationPropsCM2( TableSubToInd(m_wavelindex,locidx,altidx,scatidx), pmatrix, opticalstate );
	}
	
	if( !ok )
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableOpticalProperties_3D_UnitSphere::FillTablesAtIndex, Error configuring the Optical State at altidx[%i] locidx[%i]", (int)altidx, (int)locidx);
	}


	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::FillTablesAtIndexMultiWavel( size_t altidx, size_t locidx, SKTRAN_AtmosphericOpticalState_V21& opticalstate )
{
	bool             ok = true;
	skRTPhaseMatrix	 pmatrix;

//	double					cosangle;
//	skRTPhaseMatrix			pmatrix;

	std::vector<double> absxs;
	std::vector<double> extxs;
	std::vector<double> scattxs;

	if( m_firsttime )
	{
		nx2dArray<skRTPhaseMatrix> Ptemp;
		Ptemp.SetSize(m_wavenumber.size(), m_scatanglevector.size());
		opticalstate.CalculateMultiWaveCrossSectionsAndPhaseMatrix( m_wavenumber, &absxs, &extxs, &scattxs, m_scatanglevector, &Ptemp );
		for( size_t wavelidx = 0; wavelidx < m_wavenumber.size(); wavelidx++ )
		{
			NXASSERT((extxs[wavelidx] < 10000));
			m_totalextinction[wavelidx][locidx][altidx] = extxs[wavelidx];
			m_scatextinction[wavelidx][locidx][altidx] = scattxs[wavelidx];
			for( SKTRAN_GridIndex scatidx = 0; scatidx < m_scatteranglegrid->NumAngles(); scatidx++ )
			{
                m_scatprops->StorePolarizationPropsCM2( TableSubToInd(wavelidx,locidx,altidx,scatidx), Ptemp.At( wavelidx, scatidx), opticalstate );
			}
		}
	}

	if( !ok )
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableOpticalProperties_3D_UnitSphere::FillTablesAtIndex, Error configuring the Optical State at altidx[%i] locidx[%i]", (int)altidx, (int)locidx);
	}


	return ok;
}


inline SKTRAN_GridIndex SKTRAN_TableOpticalProperties_3D_UnitSphere::TableSubToInd ( size_t wavidx, SKTRAN_GridIndex locidx, SKTRAN_GridIndex altidx, SKTRAN_GridIndex angidx ) const
{
	return (( wavidx 
		*m_unitsphere->NumUnitVectors() + locidx)
		*m_alts->NumAltitudes() + altidx)
		*m_scatteranglegrid->NumAngles() + angidx;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::SetAltitudes		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::SetAltitudes(  const SKTRAN_GridDefOpticalPropertiesRadii_V21& alts )
{
	bool ok = true;

	alts.AddRef();
	if (m_alts != NULL) m_alts->Release();
	m_alts = &alts;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::SetUnitSphere		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::SetUnitSphere( const SKTRAN_UnitSphere_V2& unitsphere )
{
	bool ok = true;
	
	unitsphere.AddRef();
	if (m_unitsphere != NULL) m_unitsphere->Release();
	m_unitsphere = &unitsphere;

	m_inonedimmode = (m_unitsphere->NumUnitVectors() == 1);

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::SetScatterGrid( const SKTRAN_GridDefScatterAngle_V21& scatgrid )
{
	bool ok = true;

	scatgrid.AddRef();
	if (m_scatteranglegrid != NULL) m_scatteranglegrid->Release();
	m_scatteranglegrid = &scatgrid;
	m_scatanglevector.resize( m_scatteranglegrid->size() );
	for( size_t idx = 0; idx < m_scatteranglegrid->size(); idx++ )
	{
		m_scatanglevector[idx] = m_scatteranglegrid->At( idx );
	}

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::SetWavelengthGrid(const SKTRAN_GridDefWavelength& wavgrid)
{
	bool ok = true;

	wavgrid.AddRef();
	if (m_wavelengthgrid != NULL) m_wavelengthgrid->Release();
	m_wavelengthgrid = &wavgrid;
	m_wavelengtharray.resize( m_wavelengthgrid->size() );
	m_wavenumber.resize( m_wavelengthgrid->size() );
	for( size_t idx = 0; idx < m_wavelengthgrid->size(); idx++ )
	{
		m_wavelengtharray[idx] = m_wavelengthgrid->At( idx );
		m_wavenumber[idx] = 1e7 / m_wavelengthgrid->At( idx );
	}

	return ok;
}



bool SKTRAN_TableOpticalProperties_3D_UnitSphere::SetWavelengths( std::vector<double>& wavelengths )
{
	bool ok = true;
	m_wavelengtharray = wavelengths;
	m_wavenumber.resize(m_wavelengtharray.size());

	for (size_t idx = 0; idx < m_wavelengtharray.size(); idx++)
	{
		m_wavenumber[m_wavenumber.size() - idx - 1] = 1e7 / m_wavelengtharray[idx];
	}

	// insist that the wavelengths are sorted in ascending order for interpolation (there may be a better way, as this may not be necessary if there is no wavelength interpolation)
	if (wavelengths.size() > 0) ok = ok && std::is_sorted(wavelengths.begin(), wavelengths.end());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableOpticalProperties_3D_UnitSphere::SetWavelengths, Wavelength array must be sorted in ascending order.");
	}

	if (m_wavelengthgrid != NULL) m_wavelengthgrid->Release();
	SKTRAN_GridDefWavelength* grid = new SKTRAN_GridDefWavelength;
	if (wavelengths.size() > 0)
	{
		grid->ConfigureWavelengths(&wavelengths.front(), wavelengths.size());
	}
	m_wavelengthgrid = grid;
	
	return ok;
}


void SKTRAN_TableOpticalProperties_3D_UnitSphere::ReleaseResources()
{
	if( nullptr != m_unitsphere ) m_unitsphere->Release();
	if( nullptr != m_scatteranglegrid ) m_scatteranglegrid->Release();
	if( nullptr != m_wavelengthgrid ) m_wavelengthgrid->Release();
	if( nullptr != m_albedo ) m_albedo->Release();
	if( nullptr != m_alts ) m_alts->Release();
	m_unitsphere = nullptr;
	m_scatteranglegrid = nullptr;
	m_albedo = nullptr;
	m_alts = nullptr;

	m_totalextinction.clear();
	m_scatextinction.clear();
	m_wavelengtharray.clear();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::TotalExtinctionPerCM		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_TableOpticalProperties_3D_UnitSphere::TotalExtinctionPerCM( const HELIODETIC_POINT& point ) const
{
	return InterpTable( m_totalextinction[m_wavelindex], point );
}

double SKTRAN_TableOpticalProperties_3D_UnitSphere::TotalExtinctionPerCM(double wavelength, const HELIODETIC_POINT& point) const
{
	return InterpTable(m_totalextinction, wavelength, point);
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::DumpTotalExtinction() const
{
	double x,y,z;
	nxVector geopoint;
	HELIODETIC_VECTOR heliovector;
	double extinctionsum;
	nxLog::Record(NXLOG_WARNING, "Outputing Extinction\n");
	std::ofstream outfile;
	std::stringstream stream;
	stream << "c:/total_extinction.txt";
	outfile.open(stream.str().c_str());

	for( size_t locidx = 0; locidx < m_unitsphere->NumUnitVectors(); locidx++ )
	{
		heliovector.SetCoords( m_unitsphere->UnitVectorAt( locidx ).X(),
							   m_unitsphere->UnitVectorAt( locidx ).Y(),
							   m_unitsphere->UnitVectorAt( locidx ).Z() );
		geopoint = CoordinatesPtr()->HelioVectorToGeographic( heliovector );
		x = geopoint.X();
		y = geopoint.Y();
		z = geopoint.Z();
		extinctionsum = 0;
		std::cout << m_totalextinction[0][locidx][40] << std::endl;
	}
	outfile.close();
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::InterpTable		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_TableOpticalProperties_3D_UnitSphere::InterpTable( const std::vector< std::vector<double> >& table, const HELIODETIC_POINT& loc ) const
{
	double result = 0;
	
	double					altweights[2];			// assume only two altitude indices
	double					locweights[3];			// assumone three location indicies, for triangular interp
	size_t					altindices[2];
	size_t					locindices[3];

	size_t						numloc;
	size_t						numalt;
	
	CalcAltIndices( loc, altweights, altindices, numalt );
	CalcSphereIndices( loc, locweights, locindices, numloc );

	for( size_t locidx = 0; locidx < numloc; locidx++ )
	{
		for( size_t altidx = 0; altidx < numalt; altidx++ )
		{
			result += locweights[locidx] * altweights[altidx] * table[locindices[locidx]][altindices[altidx]];
		}
	}
	return result;
}

double SKTRAN_TableOpticalProperties_3D_UnitSphere::InterpTable(const std::vector< std::vector< std::vector<double> > >& table, double wavelength, const HELIODETIC_POINT& loc) const
{
	double result = 0;

	double					wlweights[2];			// assume two wavelength indices
	double					altweights[2];			// assume only two altitude indices
	double					locweights[3];			// assumone three location indicies, for triangular interp
	size_t					wlindices[2];
	size_t					altindices[2];
	size_t					locindices[3];
	size_t						numwl;
	size_t						numloc;
	size_t						numalt;

	CalcAltIndices(loc, altweights, altindices, numalt);
	CalcSphereIndices(loc, locweights, locindices, numloc);
	CalcWavelengthIndices(wavelength, wlweights, wlindices, numwl);

	for (size_t wlidx = 0; wlidx < numwl; wlidx++)
	{
		for (size_t locidx = 0; locidx < numloc; locidx++)
		{
			for (size_t altidx = 0; altidx < numalt; altidx++)
			{
				result += wlweights[wlidx] * locweights[locidx] * altweights[altidx] * table[wlindices[wlidx]][locindices[locidx]][altindices[altidx]];
			}
		}
	}
	return result;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::ScatteringExtinctionPerCM		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_TableOpticalProperties_3D_UnitSphere::ScatteringExtinctionPerCM( const HELIODETIC_POINT& point ) const
{
	return InterpTable( m_scatextinction[m_wavelindex], point );
}


double SKTRAN_TableOpticalProperties_3D_UnitSphere::ScatteringExtinctionPerCM(double wavelength, const HELIODETIC_POINT& point) const
{
	return InterpTable(m_scatextinction, wavelength, point);
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetScatteringCoefficientCM2( const HELIODETIC_POINT& point, double cosangle,  SKTRAN_PhaseMatrixScalar* scatcoeff ) const
{
	bool ok = true;

	SKTRAN_GridIndex gridindices[12];
	double           gridweights[12];
	size_t           numNonZero;

	ok = ok && GetUniquePointWeights( point, cosangle, gridindices, gridweights, numNonZero );
	ok = ok && GetScatteringCoefficientCM2_FromWeights( gridindices, gridweights, numNonZero, *scatcoeff );

	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetScatteringCoefficientCM2(double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff) const
{
	bool ok = true;

	SKTRAN_GridIndex gridindices[24];
	double           gridweights[24];
	size_t           numNonZero;

	ok = ok && GetUniquePointWeights(wavelength, point, cosangle, gridindices, gridweights, numNonZero);
	ok = ok && GetScatteringCoefficientCM2_FromWeights(gridindices, gridweights, numNonZero, *scatcoeff);

	return ok;
}


inline bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetScatteringCoefficientCM2_FromWeights( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_PhaseMatrixScalar& scatcoeff ) const
{
    return m_scatprops->GetPhaseFunctionCM2 ( gridindex, gridweight, numels, scatcoeff );
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::GetAlbedo		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetBRDF( const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf ) const

{
	bool ok = true;
	GEODETIC_INSTANT		geopoint;

	if( NULL == m_albedo )
	{
		*brdf = 0.0;
		ok = true;
	}
	else
	{
		geopoint = CoordinatesPtr()->PointToGeodetic( point, m_mjd );
		ok = m_albedo->BRDF( m_wavelen, geopoint,mu_in, mu_out, cosdphi, brdf );
	}

	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetBRDF(double wavelength, const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const
{
	bool ok = true;
	GEODETIC_INSTANT		geopoint;

	if (NULL == m_albedo)
	{
		*brdf = 0.0;
		ok = true;
	}
	else
	{
		geopoint = CoordinatesPtr()->PointToGeodetic(point, m_mjd);
		ok = m_albedo->BRDF(wavelength, geopoint, mu_in, mu_out, cosdphi, brdf);
	}

	return ok;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_TableOpticalProperties_3D_UnitSphere::GetBRDFGeodetic		2017-05-01*/
/** Returns the albedo at the specified point**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetBRDFGeodetic(const GEODETIC_INSTANT& userpoint, double mu_in, double mu_out, double cosdphi, double* brdf) const

{
	bool					ok;

	if (m_albedo == NULL)
	{
		*brdf = 0.0;
		ok = true;
	}
	else
	{
		ok = m_albedo->BRDF(m_wavelen, userpoint, mu_in, mu_out, cosdphi, brdf);
	}
	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetBRDFGeodetic(double wavelength, const GEODETIC_INSTANT& userpoint, double mu_in, double mu_out, double cosdphi, double* brdf) const
{
	bool					ok;

	if (m_albedo == NULL)
	{
		*brdf = 0.0;
		ok = true;
	}
	else
	{
		ok = m_albedo->BRDF(wavelength, userpoint, mu_in, mu_out, cosdphi, brdf);
	}
	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetEffectiveExtinctionPerCMWithHeight1( const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1 ) const
{
	if( r->ExtinctionAtCellStart( startPtIndex ) < 0 )
	{
		HELIODETIC_POINT start;
		r->LocationOfPoint( startPtIndex, &start );
		*sigma0 = TotalExtinctionPerCM( start );
		r->SetExtinction( startPtIndex, *sigma0 );
		NXASSERT( ((*sigma0 >= 0) && (*sigma0 < 1.0E5)) );
	}
	else
	{
		*sigma0 = r->ExtinctionAtCellStart( startPtIndex );
		NXASSERT( ((*sigma0 >= 0) && (*sigma0 < 1.0E5)) );
	}
	if( r->ExtinctionAtCellStart( startPtIndex + 1 ) < 0 )
	{
		HELIODETIC_POINT end;
		r->LocationOfPoint( startPtIndex + 1, &end );
		*sigma1 = TotalExtinctionPerCM( end );
		r->SetExtinction( startPtIndex + 1, *sigma1 );
		NXASSERT( ((*sigma1 >= 0) && (*sigma1 < 1.0E5)) );
	}
	else
	{
		*sigma1 = r->ExtinctionAtCellStart( startPtIndex + 1 );
		NXASSERT( ((*sigma1 >= 0) && (*sigma1 < 1.0E5)) );
	}
	return true;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetEffectiveExtinctionPerCMWithHeight1(double wavelength, const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1) const
{
	if (r->ExtinctionAtCellStart(startPtIndex) < 0)
	{
		HELIODETIC_POINT start;
		r->LocationOfPoint(startPtIndex, &start);
		*sigma0 = TotalExtinctionPerCM(wavelength, start);
		r->SetExtinction(startPtIndex, *sigma0);
		NXASSERT(((*sigma0 >= 0) && (*sigma0 < 1.0E5)));
	}
	else
	{
		*sigma0 = r->ExtinctionAtCellStart(startPtIndex);
		NXASSERT(((*sigma0 >= 0) && (*sigma0 < 1.0E5)));
	}
	if (r->ExtinctionAtCellStart(startPtIndex + 1) < 0)
	{
		HELIODETIC_POINT end;
		r->LocationOfPoint(startPtIndex + 1, &end);
		*sigma1 = TotalExtinctionPerCM(wavelength, end);
		r->SetExtinction(startPtIndex + 1, *sigma1);
		NXASSERT(((*sigma1 >= 0) && (*sigma1 < 1.0E5)));
	}
	else
	{
		*sigma1 = r->ExtinctionAtCellStart(startPtIndex + 1);
		NXASSERT(((*sigma1 >= 0) && (*sigma1 < 1.0E5)));
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::GetEffectiveExtinctionPerCMWithHeight		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetEffectiveExtinctionPerCMWithHeight1( const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1 ) const
{
	*sigma0 = TotalExtinctionPerCM( startpoint );
	*sigma1 = TotalExtinctionPerCM( endpoint   );
	NXASSERT( *sigma0 >= 0.0 );
	NXASSERT( *sigma1 >= 0.0 );
	return true;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetEffectiveExtinctionPerCMWithHeight1(double wavelength, const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1) const
{
	*sigma0 = TotalExtinctionPerCM(wavelength, startpoint);
	*sigma1 = TotalExtinctionPerCM(wavelength, endpoint);
	NXASSERT(*sigma0 >= 0.0);
	NXASSERT(*sigma1 >= 0.0);
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::IsOptionTrue		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::IsOptionTrue( SKTRAN_TableOpticalProperties_Base::OPTIONSENUM options) const
{
	return false;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::CalcAltIndices		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::CalcAltIndices( const HELIODETIC_POINT& loc, double* weights, size_t* indices, size_t& numindex ) const
{
	bool ok = true;

	double alt = loc.Altitude();
	numindex = 2;

	ok = ok && m_alts->FindBoundingIndices( alt, SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, &indices[0], &weights[0], &indices[1], &weights[1] );


	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_3D_UnitSphere::CalcSphereIndices		2013-09-09*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::CalcSphereIndices( const HELIODETIC_POINT& loc, double* weights, size_t* indices, size_t& numindex ) const
{
	bool ok = true;
	
	numindex = m_unitsphere->MaxNumInterpIndex();
	if (numindex == 1)
	{
		weights[0] = 1;
		indices[0] = 0;
		return ok;
	}
	HELIODETIC_UNITVECTOR locunit = loc.UnitVector();
	nxVector	locvector;
	locvector.SetCoords( locunit.X(),
						 locunit.Y(),
						 locunit.Z());


	ok = ok && m_unitsphere->Triangulate(locvector, indices, weights, 3, m_speedhelper[omp_get_thread_num()] );

	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::CalcScatterIndices( double cosangle, double* weights, SKTRAN_GridIndex* indices, size_t& numindex ) const
{
	bool ok = true;

	numindex = 2;

	ok = ok && m_scatteranglegrid->FindBoundingIndices( cosangle, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &indices[0], &weights[0], &indices[1], &weights[1] );

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::CalcWavelengthIndices(double wavelength, double* weights, size_t* indices, size_t& numindex) const
{
	bool ok = true;

	if (m_wavelengthgrid->NumWavelengths() > 0)
	{
		ok = ok && m_wavelengthgrid->FindBoundingIndices(wavelength, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &indices[0], &weights[0], &indices[1], &weights[1]);
		numindex = 2;
	}
	else
	{
		weights[0] = 1.0;
		indices[0] = 0;
		numindex = 1;
	}
	return ok;


}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetUniquePointWeights( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_GridIndex gridindices[12], double gridweights[12], size_t& numNonZero ) const 
{
	bool ok = true;
		
	double					altweights[2];			// assume only two altitude indices
	double					scatweights[2];
	double					locweights[3];			// assumone three location indicies, for triangular interp
	size_t					altindex[2];
	size_t					locindex[3];
	SKTRAN_GridIndex		scatindex[2];
	size_t	                numalt, numloc, numscat;
	const double            minval = 0.0;

	ok = ok && CalcAltIndices( point, altweights, altindex, numalt );
	ok = ok && CalcScatterIndices( cosAngle, scatweights, scatindex, numscat );
	ok = ok && CalcSphereIndices( point, locweights, locindex, numloc );
	//if(NULL!=pmatrix) *pmatrix = (lowerAltZeroAngle[loIndex]*(loAltWeight*loWeight) + lowerAltZeroAngle[hiIndex]*(loAltWeight*hiWeight)) + (upperAltZeroAngle[loIndex]*(hiAltWeight*loWeight) + upperAltZeroAngle[hiIndex]*(hiAltWeight*hiWeight));
    
	numNonZero = 0;
	for(int locidx=0; locidx<numloc; ++locidx)
	{
		for(int altidx=0; altidx<numalt; ++altidx)
		{
			for(int scatidx=0; scatidx<numscat; ++scatidx)
			{
				if( minval < locweights[locidx]*altweights[altidx]*scatweights[scatidx] )
				{
					gridindices[numNonZero] = TableSubToInd( m_wavelindex, locindex[locidx], altindex[altidx], scatindex[scatidx] );
					gridweights[numNonZero] = locweights[locidx]*altweights[altidx]*scatweights[scatidx];
					++numNonZero;
				}
			}
		}
	}
	

	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetUniquePointWeights(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_GridIndex gridindices[24], double gridweights[24], size_t& numNonZero) const
{
	bool ok = true;

	double					wlweights[2];			// assume two wavelength indices
	double					altweights[2];			// assume only two altitude indices
	double					scatweights[2];
	double					locweights[3];			// assumone three location indicies, for triangular interp
	size_t					wlindex[2];
	size_t					altindex[2];
	size_t					locindex[3];
	SKTRAN_GridIndex		scatindex[2];
	size_t	                numwl, numalt, numloc, numscat;
	const double            minval = 0.0;

	ok = ok && CalcWavelengthIndices(wavelength, wlweights, wlindex, numwl);
	ok = ok && CalcAltIndices(point, altweights, altindex, numalt);
	ok = ok && CalcScatterIndices(cosAngle, scatweights, scatindex, numscat);
	ok = ok && CalcSphereIndices(point, locweights, locindex, numloc);
	//if(NULL!=pmatrix) *pmatrix = (lowerAltZeroAngle[loIndex]*(loAltWeight*loWeight) + lowerAltZeroAngle[hiIndex]*(loAltWeight*hiWeight)) + (upperAltZeroAngle[loIndex]*(hiAltWeight*loWeight) + upperAltZeroAngle[hiIndex]*(hiAltWeight*hiWeight));

	numNonZero = 0;
	for (int wlidx = 0; wlidx < numwl; ++wlidx)
	{
		for (int locidx = 0; locidx < numloc; ++locidx)
		{
			for (int altidx = 0; altidx < numalt; ++altidx)
			{
				for (int scatidx = 0; scatidx < numscat; ++scatidx)
				{
					if (minval < locweights[locidx] * altweights[altidx] * scatweights[scatidx])
					{
						gridindices[numNonZero] = TableSubToInd(wlindex[wlidx], locindex[locidx], altindex[altidx], scatindex[scatidx]);
						gridweights[numNonZero] = wlweights[wlidx] * locweights[locidx] * altweights[altidx] * scatweights[scatidx];
						++numNonZero;
					}
				}
			}
		}
	}


	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetScatteringMatrixCM2 ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix   ) const 
{
	bool ok = true;

	SKTRAN_GridIndex gridindices[12];
	double           gridweights[12];
	size_t           numNonZero;

	ok = ok && GetUniquePointWeights( point, cosAngle, gridindices, gridweights, numNonZero );
	ok = ok && GetScatteringMatrixCM2_FromWeights( gridindices, gridweights, numNonZero, pmatrix );

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetScatteringMatrixCM2(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix) const
{
	bool ok = true;

	SKTRAN_GridIndex gridindices[12];
	double           gridweights[12];
	size_t           numNonZero;

	ok = ok && GetUniquePointWeights(wavelength, point, cosAngle, gridindices, gridweights, numNonZero);
	ok = ok && GetScatteringMatrixCM2_FromWeights(gridindices, gridweights, numNonZero, pmatrix);

	return ok;
}

inline bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetScatteringMatrixCM2_FromWeights( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_ScatMat_MIMSNC& pmatrix ) const
{
	bool ok = true;
	ok = ok && m_scatprops->GetScatteringMatrixCM2( gridindex, gridweight, numels, pmatrix);
	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetResultOfUnpolarizedScatterCM2 ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec ) const
{
	bool ok = true;

	SKTRAN_GridIndex gridindices[12];
	double           gridweights[12];
	size_t           numNonZero;

	ok = ok && GetUniquePointWeights( point, cosAngle, gridindices, gridweights, numNonZero );
	ok && m_scatprops->GetResultOfUnpolarizedScatterCM2( gridindices, gridweights, numNonZero, stokesvec);

	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::GetResultOfUnpolarizedScatterCM2(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec) const
{
	bool ok = true;

	SKTRAN_GridIndex gridindices[24];
	double           gridweights[24];
	size_t           numNonZero;

	ok = ok && GetUniquePointWeights(wavelength, point, cosAngle, gridindices, gridweights, numNonZero);
	ok && m_scatprops->GetResultOfUnpolarizedScatterCM2(gridindices, gridweights, numNonZero, stokesvec);

	return ok;
}


bool SKTRAN_TableOpticalProperties_3D_UnitSphere::CreateInterpolationForPoint( const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator ) const
{

    bool ok = true;


    SKTRAN_GridIndex gridindices[12];
    double gridweights[12];

    double					altweights[2];			// assume only two altitude indices
    double					locweights[3];			// assumone three location indicies, for triangular interp
    size_t					altindex[2];
    size_t					locindex[3];
    size_t	                numalt, numloc;
    const double            minval = 0.0;

    const size_t numAngles = m_scatteranglegrid->NumAngles();
    SKTRAN_ScatMat_MIMSNC pmatrix;
    std::vector<SKTRAN_ScatMat_MIMSNC> pmatrices;
    pmatrices.reserve( numAngles );

    ok = ok && CalcAltIndices( point, altweights, altindex, numalt );
    ok = ok && CalcSphereIndices( point, locweights, locindex, numloc );
    const size_t step2 = m_scatteranglegrid->NumAngles();
    const size_t step1 = m_alts->NumAltitudes()*step2;
    size_t numNonZero;

    // Weights never change 
    numNonZero = 0;
    for(int locidx=0; locidx<numloc; ++locidx)
    {
        for(int altidx=0; altidx<numalt; ++altidx)
        {
            gridweights[numNonZero] = locweights[locidx]*altweights[altidx];
            ++numNonZero;
        }
    }


    for (size_t aidx = 0; aidx < numAngles; ++aidx) {
        
        //ok = ok && CalcScatterIndices( m_scatteranglegrid->At(aidx), scatweights, scatindex, numscat ); // Always scatindex=aidx, weights[0]=1

        numNonZero = 0;
        for(int locidx=0; locidx<numloc; ++locidx)
        {
            for(int altidx=0; altidx<numalt; ++altidx)
            {
                gridindices[numNonZero] = TableSubToInd( m_wavelindex, locindex[locidx], altindex[altidx], aidx );
                ++numNonZero;
            }
        }

        ok = ok && m_scatprops->GetScatteringMatrixCM2( gridindices, gridweights, numNonZero, pmatrix );
        pmatrices.push_back( pmatrix*100.0 );
    }


    double scatteringExtinction = ScatteringExtinctionPerCM ( point );
    double totalExtinction      = TotalExtinctionPerCM      ( point );
    interpolator.InjectTable( pmatrices, m_scatteranglegrid, scatteringExtinction/totalExtinction, point );

    return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere::CreateInterpolationForPoint(double wavelength, const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator) const
{

	bool ok = true;


	SKTRAN_GridIndex gridindices[24];
	double gridweights[24];

	double					wlweights[2];
	double					altweights[2];			// assume only two altitude indices
	double					locweights[3];			// assumone three location indicies, for triangular interp
	size_t					wlindex[2];
	size_t					altindex[2];
	size_t					locindex[3];
	size_t	                numwl, numalt, numloc;
	const double            minval = 0.0;

	const size_t numAngles = m_scatteranglegrid->NumAngles();
	SKTRAN_ScatMat_MIMSNC pmatrix;
	std::vector<SKTRAN_ScatMat_MIMSNC> pmatrices;
	pmatrices.reserve(numAngles);

	ok = ok && CalcWavelengthIndices(wavelength, wlweights, wlindex, numwl);
	ok = ok && CalcAltIndices(point, altweights, altindex, numalt);
	ok = ok && CalcSphereIndices(point, locweights, locindex, numloc);
	const size_t step2 = m_scatteranglegrid->NumAngles();
	const size_t step1 = m_alts->NumAltitudes()*step2;
	size_t numNonZero;

	// Weights never change 
	numNonZero = 0;
	for (int wlidx = 0; wlidx < numwl; ++wlidx)
	{
		for (int locidx = 0; locidx < numloc; ++locidx)
		{
			for (int altidx = 0; altidx < numalt; ++altidx)
			{
				gridweights[numNonZero] = wlweights[wlidx] * locweights[locidx] * altweights[altidx];
				++numNonZero;
			}
		}
	}

	for (size_t aidx = 0; aidx < numAngles; ++aidx) {

		//ok = ok && CalcScatterIndices( m_scatteranglegrid->At(aidx), scatweights, scatindex, numscat ); // Always scatindex=aidx, weights[0]=1

		numNonZero = 0;
		for (int wlidx = 0; wlidx < numwl; ++wlidx)
		{
			for (int locidx = 0; locidx < numloc; ++locidx)
			{
				for (int altidx = 0; altidx < numalt; ++altidx)
				{
					gridindices[numNonZero] = TableSubToInd(wlindex[wlidx], locindex[locidx], altindex[altidx], aidx);
					++numNonZero;
				}
			}
		}
		ok = ok && m_scatprops->GetScatteringMatrixCM2(gridindices, gridweights, numNonZero, pmatrix);
		pmatrices.push_back(pmatrix*100.0);
	}
	
	double scatteringExtinction = ScatteringExtinctionPerCM(wavelength, point);
	double totalExtinction = TotalExtinctionPerCM(wavelength, point);
	interpolator.InjectTable(pmatrices, m_scatteranglegrid, scatteringExtinction / totalExtinction, point);

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant::CalcAltIndices(const HELIODETIC_POINT & loc, double * weights, size_t * indices, size_t & numindex) const
{
	bool ok = true;

	double alt = loc.Altitude();
	numindex = 1;

	m_alts->IndexOfPointBelowOrEqual(alt, &indices[0]);
	weights[0] = 1.0;

	return ok;
}

bool SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant::ConfigureGeometry(const SKTRAN_Specifications_Base * specs)
{
	bool ok = true;
	ok = ok && SKTRAN_TableOpticalProperties_3D_UnitSphere::ConfigureGeometry(specs);
	if (!m_inonedimmode)
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant::ConfigureGeometry, 3D mode not supported.");
		ok = false;
	}
	return ok;
}

