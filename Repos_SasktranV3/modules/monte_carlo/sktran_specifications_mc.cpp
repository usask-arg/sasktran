#include "include/sktran_montecarlo_internals.h"

SKTRAN_Specifications_MC::SKTRAN_Specifications_MC()
{
	m_raytracingshells         = NULL;
	m_amfshells				   = NULL;
	m_scatteranglegrid         = NULL;
	m_opticalpropradii         = NULL;
    m_opticalunitsphere        = nullptr;
	m_wavelengthgrid		   = NULL;
    m_sunposition.SetInvalid();
	
	ConfigureDefaults( );
}

SKTRAN_Specifications_MC::~SKTRAN_Specifications_MC()
{
	ReleaseResources();
}

void SKTRAN_Specifications_MC::ReleaseResources()
{
	m_hasBeenFinalized = false;
	ReleaseGrids();
	if(nullptr!=m_opticalunitsphere) m_opticalunitsphere->Release(); m_opticalunitsphere=nullptr;
}


void SKTRAN_Specifications_MC::ReleaseGrids()
{
	m_hasBeenFinalized = false;
	if(NULL!=m_scatteranglegrid) m_scatteranglegrid->Release(); m_scatteranglegrid=NULL;
//	if(NULL!=m_raytracingshells) m_raytracingshells->Release(); m_raytracingshells=NULL;
	if(NULL!=m_opticalpropradii) m_opticalpropradii->Release(); m_opticalpropradii=NULL;
	if (NULL != m_wavelengthgrid) m_wavelengthgrid->Release(); m_wavelengthgrid = NULL;
}

bool SKTRAN_Specifications_MC::SetSunGeographicPosition( const nxVector& sun )
{
	bool ok = true;

	m_sunposition = sun.UnitVector();

    return ok;
}

bool SKTRAN_Specifications_MC::GetSun( nxVector* target ) const
{
	bool ok = true;

	*target = m_sunposition;

	return ok;
}

bool SKTRAN_Specifications_MC::Allocate()
{
	bool ok = true;

	m_raytracingshells.reset( new SKTRAN_GridDefRayTracingShells_V21);
	m_amfshells.reset( new SKTRAN_GridDefAirMassFactorShells);
	m_scatteranglegrid = new SKTRAN_GridDefScatterAngle_V21;
	m_opticalpropradii = new SKTRAN_GridDefOpticalPropertiesRadii_V21;
	m_wavelengthgrid = new SKTRAN_GridDefWavelength;
	if(NULL!=m_raytracingshells) m_raytracingshells->SetStatic();
	if(NULL!= m_amfshells) m_amfshells->SetStatic();
	if(NULL!=m_scatteranglegrid) m_scatteranglegrid->AddRef();
	if(NULL!=m_opticalpropradii) m_opticalpropradii->AddRef();
	if (NULL != m_wavelengthgrid) m_wavelengthgrid->AddRef();
	
//	ok = ok && NULL!=m_raytracingshells;
	ok = ok && NULL != m_amfshells;
	ok = ok && NULL!=m_scatteranglegrid;
	ok = ok && NULL!=m_opticalpropradii;
	ok = ok && NULL != m_wavelengthgrid;

	if(!ok){
		ReleaseGrids();
		nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::Allocate, Could not allocate grids.");
	}

	return ok;
}

bool SKTRAN_Specifications_MC::ConfigureDefaults()
{
	bool ok = true;
	

	m_mcEngineType             = EngineType::mc;
	m_sunType                  = SunType::point;
    m_emissionTableType        = EmissionTableType::doNothing;
	m_solarTableType           = SolarTableType::dim2; // TODO: Should change this to 3d table?
	//m_thermalTableType         = ThermalType::none;
	m_solarRayTracerType       = RayTracerType::generic;
	m_LOSRayTracerType         = RayTracerType::generic;
	m_MSRayTracerType          = RayTracerType::generic;
	m_optPropsIntType          = OptPropIntType::straight;
	m_optTableType             = OptTableType::dim1;
	m_photonLogType            = LogType::none;
	m_symmetryType             = SymmetryType::none;
	m_polType                  = PolType::none;
	m_atmosphereHasDelta       = AtmoHasDelta::no;
	m_scatterType			   = ScatterType::elastic;
	m_wavelengthType		   = WavelengthType::single;
	m_secondary				   = SecondaryOutput::none;

	m_numPhotonsPerLOS         .resize(1,10000);		// Default to 1% precision
	m_precisionMC              .resize(1, 0.01);
	m_numRayTracingAlts        = 200 + 1;
	m_numOptPropAlts           = 200 + 1;
	m_numSolTableCosSza        = 1801;
	m_thread_chunkSize         = 1;			// Allow threads to receive very small chunks
	m_rngSeed                  = 0;			// Set to zero to seed RNG_i with seed(clock,i)
	m_scattAngleResolution     = 0.5;
	m_TOAHeight                = 100000.0;
	m_curvedRayStepSize        = m_TOAHeight / (m_numRayTracingAlts-1); //1000.0;
	m_adaptivemaxopticaldepth  = 0.1;
	m_adaptiveminratio         = 0.9;
	m_minRelScatterWeight      = 0.0;		// Default to no truncation
	m_scatterPosRes            = 50.0;
	m_sineSunApexAngle         = 0.0;	// Start with point sun, and make user give value for non-zero point sun
	m_minFractionHigherOrder   = std::vector<double>(1, 0.1);
    m_solarTableAltDelta       = 500.0;
    m_groundShiftAltitude      = 0.0;
	m_threads_allowDynamic = true;
	m_chunkSize = 1;
	m_nadirreferencepointonground = false;
	m_referencepoint = std::vector<double>();

	m_amfSpeciesHandle		   = SKCLIMATOLOGY_UNDEFINED;
	m_manualAMFHeights.resize(101);
	for (size_t amfidx = 0; amfidx < 101; amfidx++)
		m_manualAMFHeights[amfidx] = amfidx * 1e3;

	m_hasBeenFinalized         = false;		// FinalizeSpecs() needs to be called after ConfigureDefaults

	return ok;
}

bool SKTRAN_Specifications_MC::FinalizeSpecs()
{
	bool ok = true;
	bool ok_nonessential = true;

	//double deltaShell;
	//double *rayShellAlts, *optShellAlts, *amfShellAlts;

	std::vector<double> heights;
	bool uniform;
	
	ReleaseGrids();
	ok = ok && Allocate();

	// configure the shells for ray tracing
	ok = ok && GetProfileAlts(ProfileType::rayTracing, heights, uniform);
	ok = ok && m_raytracingshells->ConfigureHeights(heights);
	if (uniform) ok_nonessential = ok_nonessential && m_raytracingshells->SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Problem setting ray tracing shell altitude grid.");

	// check outputQuantity settings and configure the AMF shells
	switch (m_secondary)
	{
		case SecondaryOutput::none:
			m_amfshells = nullptr;
			break;
		case SecondaryOutput::lengthAMF:
			ok = ok && m_scatterType == ScatterType::elastic;
			ok = ok && GetProfileAlts(ProfileType::airMassFactor, heights, uniform);
			ok = ok && m_amfshells->ConfigureHeights(heights, GetGroundShiftAlt(), GetTOAHeight());
			if (uniform) ok_nonessential = ok_nonessential && m_amfshells->SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Problem setting air mass factor shell altitude grid.");
			if (m_minFractionHigherOrder[0] < 1) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::FinalizeSpecs, minFractionHigherOrder should be set to 1 when calculating AMFs.");
			break;
		case SecondaryOutput::opticalDepthAMF:
			ok = ok && m_scatterType == ScatterType::elastic;
			ok = ok && GetProfileAlts(ProfileType::airMassFactor, heights, uniform);
			ok = ok && m_amfshells->ConfigureHeights(heights, GetGroundShiftAlt(), GetTOAHeight());
			if (uniform) ok_nonessential = ok_nonessential && m_amfshells->SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Problem setting air mass factor shell altitude grid.");
			if (m_minFractionHigherOrder[0] < 1) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::FinalizeSpecs, minFractionHigherOrder should be set to 1 when calculating AMFs.");
			ok = ok && !(m_amfSpeciesHandle == SKCLIMATOLOGY_UNDEFINED);
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, no AMF species was specified.");
			break;
		case SecondaryOutput::ringSpectrum:
			m_amfshells = nullptr;
			ok = ok && m_scatterType == ScatterType::both || m_scatterType == ScatterType::manualBoth;
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, the SecondaryOutput 'ringSpectrum' requires the ScatterType to be 'both' or 'manualBoth'.");
			break;
		case SecondaryOutput::fillingInParameter:
			m_amfshells = nullptr;
			ok = ok && m_scatterType == ScatterType::both || m_scatterType == ScatterType::manualBoth;
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, the SecondaryOutput 'fillingInParameter' requires the ScatterType to be 'both' or 'manualBoth'.");
			break;
		case SecondaryOutput::elasticRaman:
			m_amfshells = nullptr;
			ok = ok && m_scatterType == ScatterType::both || m_scatterType == ScatterType::manualBoth;
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, the SecondaryOutput 'elasticRaman' requires the ScatterType to be 'both' or 'manualBoth'.");
			break;
		default:
			ok = false;
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, unknown SecondaryOutput.");
			break;
	}

	// configure the altitudes for optical properties
	ok = ok && GetProfileAlts(ProfileType::opticalProperties, heights, uniform);
	ok = ok && m_opticalpropradii->ConfigureAltitudes(&heights[0], heights.size());
	if (uniform) ok_nonessential = ok_nonessential && m_opticalpropradii->SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Problem setting optical properties altitude grid.");

	ok = ok && m_scatteranglegrid->Configure(m_scattAngleResolution);
	ok_nonessential = ok_nonessential && m_scatteranglegrid->SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);

	// configure the wavelength grid
	switch (m_scatterType)
	{
		case ScatterType::both:				// same configuration as inelastic, intentional fallthrough
		case ScatterType::manualBoth:		// same configuration as inelastic, intentional fallthrough
		case ScatterType::manualInelastic:	// same configuration as inelastic, intentional fallthrough
		case ScatterType::inelastic:
			ok = ok && m_opticalpropertieswavelengths.size() > 0;
			ok = ok && m_wavelengthgrid->ConfigureWavelengths(&m_opticalpropertieswavelengths[0], m_opticalpropertieswavelengths.size());
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Inelastic scattering mode requires you to specify an optical wavelength grid");

			ok = ok && m_solarSpectrum.size() == m_wavelengthgrid->NumWavelengths();
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Inelastic scattering mode requires you to specify an incident solar spectrum on the optical wavelength grid");
			break;
		case ScatterType::elastic:
			//if (m_wavelengtharray.size() > 0) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::FinalizeSpecs, The specified wavelength array is ignored in elastic scattering mode.");
			if (m_opticalpropertieswavelengths.size() > 0)
			{
				ok = ok && m_wavelengthgrid->ConfigureWavelengths(&m_opticalpropertieswavelengths[0], m_opticalpropertieswavelengths.size()); // This will just make it slightly less efficient, so I should remove it later, but it's a good test for now
				if (m_solarSpectrum.size() > 0)
				{
					ok = ok && m_solarSpectrum.size() == m_wavelengthgrid->NumWavelengths();
					if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Incident solar spectrum does not match the optical wavelength grid");
				}
				else
				{
					m_solarSpectrum.resize(m_wavelengthgrid->NumWavelengths());
					std::fill(m_solarSpectrum.begin(), m_solarSpectrum.end(), 1.0);
					nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::FinalizeSpecs, No incident solar spectrum was provided so values of 1.0 were assumed");
				}
			}
			break;
		default:
			ok = false;
			nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Error configuring the wavelength grid.");
			break;
	}
	
	// check for consistency in the elastic/inelastic specs
	switch (m_scatterType)
	{
		case ScatterType::elastic:
			ok = ok && m_minFractionHigherOrder.size() == 1;
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, minFractionHigherOrder must be a single value in elastic scatter mode");
			break;
		case ScatterType::both: // intentional fallthrough
		case ScatterType::inelastic:
			ok = ok && m_minFractionHigherOrder.size() == 1;
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, minFractionHigherOrder must be a single value in non-optimized inelastic scatter mode");
			ok_nonessential = m_minFractionHigherOrder.front() == 1.0;
			if (!ok_nonessential) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::FinalizeSpecs, minFractionHigherOrder should be set to 1 in non-optimized inelastic scatter mode");
			break;
		case ScatterType::manualBoth: // intentional fallthrough
		case ScatterType::manualInelastic:
			ok = ok && m_minFractionHigherOrder.size() > 1;
			if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, in optimized inelastic mode minFractionHigherOrder must have size > 1");
		default:
			break;
	}

	// check for consistency in the simultaneous-wavelength specs
	if (m_wavelengthType == WavelengthType::simultaneous)
	{
		ok = ok && m_secondary != SecondaryOutput::lengthAMF && m_secondary != SecondaryOutput::opticalDepthAMF;
		if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Simultaneous wavelength mode is not compatible with AMFs.");
		ok = ok && m_radiancewavelengths.size() > 1;
		if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, More than one wavelength must be specified for simultaneous wavelength mode.");
		auto it = std::find(m_radiancewavelengths.begin(), m_radiancewavelengths.end(), m_primaryWavelength);
		ok = ok && it != m_radiancewavelengths.end();
		if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Invalid primary wavelength for simultaneous wavelength mode.");
		if (ok) m_primaryWavelengthIndex = it - m_radiancewavelengths.begin();
		ok = ok && m_opticalpropertieswavelengths.size() > 1;
		if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::FinalizeSpecs, Optical property wavelengths must be set and its length must be at least 2 for simultaenous wavelength mode.");
		ok_nonessential = *std::min_element(m_opticalpropertieswavelengths.begin(), m_opticalpropertieswavelengths.end()) <= *std::min_element(m_radiancewavelengths.begin(), m_radiancewavelengths.end()) && *std::max_element(m_opticalpropertieswavelengths.begin(), m_opticalpropertieswavelengths.end()) >= *std::max_element(m_radiancewavelengths.begin(), m_radiancewavelengths.end());
		if (!ok_nonessential) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::FinalizeSpeces, Radiance will be calculated at wavelengths outside the specified optical property wavelengths; extraoplation will occur.");
	}

	m_hasBeenFinalized = ok;
	if( !m_hasBeenFinalized ) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::FinalizeSpecs, Could not finalize specs.");

	return ok;

}

bool SKTRAN_Specifications_MC::CreateSun_Point( std::unique_ptr<SKTRAN_Sun_Base>& target ) const
{
	bool ok = true;
	
	std::unique_ptr<SKTRAN_Sun_Point> sunp( new SKTRAN_Sun_Point );
	ok = ok && nullptr!=sunp;
	if(ok){
        target = std::move(sunp);
    }else {
        nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateSun_Point, Error creating sun.");
    }

	return ok;
}

bool SKTRAN_Specifications_MC::CreateSun_RandomDisc( std::unique_ptr<SKTRAN_Sun_Base>& target, const std::vector<SKTRAN_RNG>& rngs, size_t numthreads ) const
{
	bool ok = true;

	std::unique_ptr<SKTRAN_Sun_RandomDisc> sunp ( new SKTRAN_Sun_RandomDisc );
	ok = ok && nullptr!=sunp;
	ok = ok && sunp->Initialize(GetSineSunApexAngle(), rngs, numthreads);
	if(ok){
		target = move(sunp);
    } else{
        nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateSun_RandomDisc, Error creating sun.");
    }

	return ok;
}

bool SKTRAN_Specifications_MC::CreateSun( std::unique_ptr<SKTRAN_Sun_Base>& target, const std::vector<SKTRAN_RNG>& rngs, size_t numthreads ) const
{
	bool ok = true;
    
	switch( m_sunType )
	{
		case SunType::point:
			ok = ok && CreateSun_Point( target );
			break;
		case SunType::randomDisc:
			ok = ok && CreateSun_RandomDisc( target, rngs, numthreads );
			break;
		default:
			ok = false;
			break;
	}

	return ok;
}

bool SKTRAN_Specifications_MC::CreateConfigurationManager( std::unique_ptr<SKTRAN_ConfigurationManager_MC>& target) const 
{
	bool ok = true;

	std::unique_ptr<SKTRAN_ConfigurationManager_MC> cman( new SKTRAN_ConfigurationManager_MC );
    
	ok = ok && nullptr!=cman;

	if(ok) target = std::move(cman);

	return ok;

}


bool SKTRAN_Specifications_MC::SetPrecisionMC	(double d)
{
	return SetPrecisionMC( std::vector<double>( 1, d ) );
}

bool SKTRAN_Specifications_MC::SetPrecisionMC ( const std::vector<double>& dvec )
{
	bool ok = true;
	ok = ok && std::all_of(dvec.begin(), dvec.end(), [](double d){ return 0.0<=d;});
	if(ok){
		m_precisionMC.resize(dvec.size());
		std::copy(dvec.begin(), dvec.end(), m_precisionMC.begin());
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetPrecisionMC, Must have non-negative precision.");
	}

	return ok;
}

bool SKTRAN_Specifications_MC::SetScattAngleResolution	(double d)
{
	bool ok = true;
	ok = d > 0.0;
	if(ok){
		m_scattAngleResolution = d;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetScattAngleResolution, Must have positive scatter angle resolution.");
	}
	return ok;
}

bool SKTRAN_Specifications_MC::SetCurvedRayStepSize		(double d)
{
	bool ok = true;
	ok = d > 0.0;
	if(ok){
		m_curvedRayStepSize = d;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetCurvedRayStepSize, Cannot have negative curved ray step size.");
	}
	return ok;
}


bool SKTRAN_Specifications_MC::SetNumPhotonsPerLOS (size_t n)
{
	return SetNumPhotonsPerLOS( std::vector<size_t>(1,n) );
}

bool SKTRAN_Specifications_MC::SetNumPhotonsPerLOS (const std::vector<size_t>& nvec)
{
	bool ok = true;
	m_numPhotonsPerLOS.resize(nvec.size());
	std::copy(nvec.begin(), nvec.end(), m_numPhotonsPerLOS.begin());
	return ok;
}

bool SKTRAN_Specifications_MC::SetPrimaryWavelength(double d)
{
	m_primaryWavelength = d;
	return true;
}

bool SKTRAN_Specifications_MC::SetAdaptOptDepthMax     (double d)
{
	bool ok = true;
	ok = d > 0.0;
	if(ok){
		m_adaptivemaxopticaldepth = d;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetNumPhotonsPerLOS, Cannot have negative adaptive otpical depth.");
	}
	return ok;
}


bool SKTRAN_Specifications_MC::SetAdaptOptDepthMinRatio (double d)
{
	bool ok = true;
	ok = d > 0.0;
	if(ok){
		m_adaptiveminratio = d;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetAdaptOptDepthMinRatio, Cannot have negative adaptive otpical depth.");
	}
	return ok;
}

bool SKTRAN_Specifications_MC::SetTOAHeight         (double d)
{
	bool ok = true;
	ok = d > 0.0;
	if(ok){
		m_TOAHeight = d;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetTOAHeight, Cannot have negative top of atmosphere height.");
	}
	return ok;
}


bool SKTRAN_Specifications_MC::SetMinimumRelPathWeight(double d)
{
	bool ok = true;
	ok = 0.0 <= d && d <= 1.0;
	if(ok){
		m_minRelScatterWeight = d;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetMinimumRelPathWeight, Cannot have relative path weight outside [0,1], use 0.0 for no truncation.");
	}

	return ok;
}

bool SKTRAN_Specifications_MC::SetScatterPositionRes (double d)
{
	bool ok = true;
	ok = d >= 0.0;
	if(ok){
		m_scatterPosRes = d;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetScatterPositionRes, Cannot have a negative scatter position resolution.");
	}

	return ok;
}

bool SKTRAN_Specifications_MC::SetSineSunApexAngle (double d)
{
	bool ok = true;
	ok = d >= 0.0;
	if(ok){
		m_sineSunApexAngle = d;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetSineSunApexAngle, Cannot have a negative sun aperture angle.");
	}

	return ok;
}

bool SKTRAN_Specifications_MC::SetMinFractionHigherOrder(double d)
{
	bool ok = true;
	ok = d >= 0.0;
	if (ok) {
		m_minFractionHigherOrder = std::vector<double>(1, d);
	}
	else {
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetMinFractionHigherOrder, Cannot have a negative fraction higher order (can use zero though).");
	}

	return ok;
}

bool SKTRAN_Specifications_MC::SetMinFractionHigherOrder (const std::vector<double>& d)
{
	bool ok = true;
	double s = 0;
	for (auto&& f : d)
	{
		ok = f >= 0.0;
		s += f;
	}
	ok = ok && s <= 1.0;

	if(ok){
		m_minFractionHigherOrder = d;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetMinFractionHigherOrder, Cannot have a negative fraction higher order (can use zero though).");
	}

	return ok;
}

bool SKTRAN_Specifications_MC::SetSolarTableAltDelta( double d)
{
    bool ok = true; 

    ok = ok && d > 0.0;
    if(ok){
        m_solarTableAltDelta      = d; 
    } else{
        nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetSolarTableAltDelta, Cannot have negative solar table altitude step size.");
    }
 
   return ok;
}


bool SKTRAN_Specifications_MC::SetGroundShiftAlt(double d)
{
    bool ok = true; 

    ok = ok && d >= 0.0;
    if(ok){
        m_groundShiftAltitude = d; 
    } else{
        nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetGroundShiftAlt, Cannot have negative ground shift altitude.");
    }

    return ok;
}


bool SKTRAN_Specifications_MC::CreateEmissionTable( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_EmissionTable_Base** target ) const 
{
	bool ok = true;

	ok = ok && NULL!=target;
	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateEmissionTable, Must give valid pointer reference.");
	ok = ok && NULL==*target;
	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateEmissionTable, User should clear mem space, don't want to lose a reference.");

	switch(m_emissionTableType)
	{
    case EmissionTableType::doNothing:
        ok = ok && CreateEmissionTable_DoNothing( target );
        break;
	case EmissionTableType::dim1:
		ok = ok && CreateEmissionTable_1DTable( coords, target );
		break;
	default:
		ok = false;
		break;
	}

	return ok;
}


bool SKTRAN_Specifications_MC::CreateEmissionTable_DoNothing ( SKTRAN_EmissionTable_Base** target ) const
{
    bool ok = true;

    ok = ok && NULL!=target;
	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateEmissionTable_DoNothing, Must give valid pointer reference.");
	ok = ok && NULL==*target;
	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateEmissionTable_DoNothing, User should clear mem space, don't want to lose a reference.");

    *target = new SKTRAN_EmissionTable_DoNothing;
	ok = ok && NULL!=target;
	if(ok) (*target)->AddRef();

    return ok;
}


bool SKTRAN_Specifications_MC::CreateEmissionTable_1DTable   ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_EmissionTable_Base** target ) const
{
    bool ok = true;
	bool ok_nonessential = true;

    ok = ok && NULL!=target;
	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateEmissionTable_1DTable, Must give valid pointer reference.");
	ok = ok && NULL==*target;
	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateEmissionTable_1DTable, User should clear mem space, don't want to lose a reference.");

    if(ok){
        SKTRAN_EmissionTable_1D* table = new SKTRAN_EmissionTable_1D;
   	    ok = ok && NULL!=table;

		// make some heights
		SKTRAN_GridDefRayTracingShells_V21 altgrid; 
		std::vector<double> heights;
		bool uniform;
		ok = ok && GetProfileAlts(ProfileType::solarTable, heights, uniform);
		ok = ok && altgrid.ConfigureHeights(heights);
		if (uniform) ok_nonessential = ok_nonessential && altgrid.SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
		if (!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateEmissionTable_1DTable, There was a problem creating the altitude grid.");

        GEODETIC_INSTANT point;
        point.latitude  = coords->ReferencePtLatitude();
        point.longitude = coords->ReferencePtLongitude();
        point.heightm   = 0.0;
        point.mjd       = coords->ReferencePointMJD();
        //ok = ok && table->SetGeometry( altgrid, point );
        nxLog::Record( NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateEmissionTable_1DTable, New implementation of geometry initialization not yet implemented." ); ok = false; 

	    if(ok){
            *target = table;
            (*target)->AddRef();
        } else{
            delete table;
        }
    }

    return ok;
}

bool SKTRAN_Specifications_MC::CreateSolarTransmissionTable(std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_SolarTransmission_Base** target, size_t numthreads ) const {
	bool ok = true;

	ok = ok && NULL!=target;
	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateSolarTransmissionTable, Must give valid pointer reference.");
	ok = ok && NULL==*target;
	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateSolarTransmissionTable, User should clear mem space, don't want to lose a reference.");

//	SKTRAN_Sun_Base* sunp = NULL;

	switch (m_scatterType)
	{
	case ScatterType::elastic:
		switch (m_solarTableType)
		{
		case SolarTableType::doNothing:
			ok = ok && CreateSolarTransmissionTable_DoNothing(target);
			break;
		case SolarTableType::noTable:
			ok = ok && CreateSolarTransmissionTable_NoTable(target);
			break;
		case SolarTableType::dim2:
			ok = ok && CreateSolarTransmissionTable_2DTable(coords, target);
			break;
		case SolarTableType::dim3:
			ok = ok && CreateSolarTransmissionTable_3DTable(coords, target);
			break;
		default:
			ok = false;
			break;
		}
		break;
	case ScatterType::manualInelastic: // intentional fallthrough
	case ScatterType::inelastic:
		switch (m_solarTableType)
		{
		case SolarTableType::doNothing:
			ok = ok && CreateSolarTransmissionTable_DoNothing(target);
			break;
		case SolarTableType::noTable:
			ok = ok && CreateSolarTransmissionTable_Inelastic_NoTable(target);
			break;
		case SolarTableType::dim2:
			ok = false;
			break;
		case SolarTableType::dim3:
			ok = false;
			break;
		default:
			ok = false;
			break;
		}
		break;
	case ScatterType::manualBoth: // intentional fallthrough
	case ScatterType::both:
		switch (m_solarTableType)
		{
		case SolarTableType::doNothing:
			ok = ok && CreateSolarTransmissionTable_DoNothing(target);
			break;
		case SolarTableType::noTable:
			ok = ok && CreateSolarTransmissionTable_Ring_NoTable(target);
			break;
		case SolarTableType::dim2:
			ok = false;
			break;
		case SolarTableType::dim3:
			ok = false;
			break;
		default:
			ok = false;
			break;
		}
		break;
	default:
		ok = false;
		break;
	}

	return ok;
}


bool SKTRAN_Specifications_MC::CreateSolarTransmissionTable_DoNothing(SKTRAN_SolarTransmission_Base** target) const {
	bool ok = true;

	*target = new SKTRAN_SolarTransmission_DoNothing;
	ok = ok && NULL!=target;
	if(ok) (*target)->AddRef();

	return ok;
}


bool SKTRAN_Specifications_MC::CreateSolarTransmissionTable_NoTable(SKTRAN_SolarTransmission_Base** target) const {
	bool ok = true;

	if (m_solarSpectrum.size() > 0)
	{
		SKTRAN_SolarTransmission_NoTable_reuseRays_SolarSpectrum_MC* solartable = new SKTRAN_SolarTransmission_NoTable_reuseRays_SolarSpectrum_MC;
		ok = ok && solartable->SetSolarSpectrum(m_solarSpectrum);
		ok = ok && solartable->SetWavelengthGrid(m_wavelengthgrid);

		if (ok)
		{
			*target = solartable;
			(*target)->AddRef();
		}
		else
		{
			delete solartable;
		}
	}
	else
	{
		SKTRAN_SolarTransmission_NoTable_reuseRays_MC* solartable = new SKTRAN_SolarTransmission_NoTable_reuseRays_MC;

		if (ok)
		{
			*target = solartable;
			(*target)->AddRef();
		}
		else
		{
			delete solartable;
		}
	}
	
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Specifications_MC::CreateSolarTransmissionTable_2DTable		 2014- 11- 21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Specifications_MC::CreateSolarTransmissionTable_2DTable( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_SolarTransmission_Base** target) const 
{
    bool ok = true;
	bool ok_nonessential = true;
	// make some heights
	SKTRAN_GridDefRayTracingShells_V21 altgrid;
	std::vector<double> heights;
	bool uniform;
	ok = ok && GetProfileAlts(ProfileType::solarTable, heights, uniform);
	ok = ok && altgrid.ConfigureHeights(heights);
	if (uniform) ok_nonessential = ok_nonessential && altgrid.SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
	if (!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateSolarTransmissionTable_2DTable, There was a problem creating the altitude grid.");

    // and some cossza
    SKTRAN_GridDefCosSZA_V21	cosszagrid;
    cosszagrid.AllocateGridArray( m_numSolTableCosSza );
    ok = ok && 0 < cosszagrid.NumAngles();
    if(ok){
        double delta = 2.0 / double(cosszagrid.NumAngles() - 1);
        for( size_t szaidx = 0; szaidx < cosszagrid.NumAngles(); szaidx++ )
        {
            cosszagrid.AtVar( szaidx ) = -1.0 + szaidx*delta;
        }
    }

    SKTRAN_SolarTransmission_2D* st = new SKTRAN_SolarTransmission_2D();
    ok = ok && NULL!=st;
    ok = ok && st->SetGeometry( coords, altgrid, cosszagrid );
    if(ok){
        *target = st;
        (*target)->AddRef();
    } else{
        delete st;
    }

    return ok;
}


bool SKTRAN_Specifications_MC::CreateSolarTransmissionTable_3DTable( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_SolarTransmission_Base** target) const 
{
    bool ok = true;
	bool ok_nonessential = true;

    SKTRAN_Sun_Point* sun = new SKTRAN_Sun_Point;
	
	// make some heights
	SKTRAN_GridDefRayTracingShells_V21 altgrid;
	std::vector<double> heights;
	bool uniform;
	ok = ok && GetProfileAlts(ProfileType::solarTable, heights, uniform);
	ok = ok && altgrid.ConfigureHeights(heights);
	if (uniform) ok_nonessential = ok_nonessential && altgrid.SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
	if (!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateSolarTransmissionTable_3DTable, There was a problem creating the altitude grid.");

    // and some cossza
    SKTRAN_GridDefCosSZA_V21	cosszagrid;

    double refcossza = coords->ReferencePoint(0).CosSZA();
    double cosszadiff = 0.001;
    
    double maxcossza = 1;
    double mincossza = -1;
 
    size_t		numcossza = 1801;
    cosszagrid.AllocateGridArray( numcossza );
    for( size_t szaidx = 0; szaidx < numcossza; szaidx++ )
    {
        cosszagrid.AtVar( szaidx ) = mincossza + ( (maxcossza - mincossza) * double(szaidx) / (double(numcossza) -1 ));
    }
    cosszagrid.SetGridSearchMode( SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM );

    // Create solar table longitude grid
    SKTRAN_GridDefSLON_V21		slongrid;
    ok = ok && ConfigureSLonGrid( slongrid, coords );
 

    //std::unique_ptr<SKTRAN_SolarTransmission_3D> solartable ( new SKTRAN_SolarTransmission_3D );
    SKTRAN_SolarTransmission_3D* solartable  =  new SKTRAN_SolarTransmission_3D(false);
    ok = ok && solartable->InitializeGeometry( coords, altgrid, cosszagrid, slongrid );
    //ok = ok && solartable->ConfigureOptical( m_solarrayfactory, m_optintegrator.get() );
    //ok = ok && solartable->SetSun( sun );
    //ok = ok && solartable->FillTable();

    if(ok){
        //*target = solartable.get();
        *target = solartable;
        (*target)->AddRef();
    } else{
        delete solartable;
    }

    return ok;
}

bool SKTRAN_Specifications_MC::CreateSolarTransmissionTable_Inelastic_NoTable(SKTRAN_SolarTransmission_Base** target) const
{
	bool ok = true;

	SKTRAN_SolarTransmission_Inelastic_MC* solartable = new SKTRAN_SolarTransmission_Inelastic_MC;
	ok = ok && solartable->SetSolarSpectrum(m_solarSpectrum);
	ok = ok && solartable->SetWavelengthGrid(m_wavelengthgrid);

	if (ok)
	{
		*target = solartable;
		(*target)->AddRef();
	}
	else
	{
		delete solartable;
	}

	return ok;
}

bool SKTRAN_Specifications_MC::CreateSolarTransmissionTable_Ring_NoTable(SKTRAN_SolarTransmission_Base ** target) const
{
	bool ok = true;

	SKTRAN_SolarTransmission_Ring_MC* solartable = new SKTRAN_SolarTransmission_Ring_MC;
	ok = ok && solartable->SetSolarSpectrum(m_solarSpectrum);
	ok = ok && solartable->SetWavelengthGrid(m_wavelengthgrid);

	if (ok)
	{
		*target = solartable;
		(*target)->AddRef();
	}
	else
	{
		delete solartable;
	}

	return ok;
}

bool SKTRAN_Specifications_MC::ConfigureSLonGrid( SKTRAN_GridDefSLON_V21& slongrid, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords ) const		
{		
    bool ok = true; 		
    // Note that if we are in SZA table mode, In3dMode() returns false		
    if( m_optTableType == OptTableType::dim1 || m_optTableType == OptTableType::dim1_constant)		
    {		
        // just need one slon at the TP		
        ok = ok && slongrid.AllocateGridArray( 1 );		
        slongrid.AtVar( 0 ) = 0.0;		
    } else{		
        //if( m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere )		
        //{		
        //    // slon is going to vary in each cone		
        //    ret.AllocateGridArray( m_numcones*2 + 1 );		
        //    ret.AtVar( 0 ) = m_numcones * m_coneanglesep * (-1) * nxmath::Pi / 180;		
        //    for( size_t idx = 1; idx < m_numcones*2 + 1; idx++ )		
        //    {		
        //        ret.AtVar( idx ) = ret.AtVar( idx-1 ) + m_coneanglesep * nxmath::Pi / 180;		
        //    }		
        //}		
        //if( m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_LOSPlane || m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_SZA )		
        //{		
        //    // TODO: calculate actual slon grid		
        //    // worst case is when the angle grid is the slon direction		
        //    ret.AllocateGridArray( m_anglegrid.size() );		
        //    for( size_t idx = 0; idx < m_anglegrid.size(); idx++ )		
        //    {		
        //        ret.AtVar( idx ) = m_anglegrid[idx] * nxmath::Pi / 180;		
        //    }		
        //}		
        ok = false;		
        nxLog::Record(NXLOG_ERROR,"SKTRAN_Specifications_MC::ConfigureSLonGrid, Haven't implemented the 'cones', etc for 3d solar transmission table in MC. _noTable should work.");		
    }		
    if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::ConfigureSLonGrid, This needs to be common among all engines, and I hate the way it's implemented now. Rework this later.");		
    		
    return ok;		
}		


//bool SKTRAN_Specifications_MC::CreateThermalEmissionTable( SKTRAN_ThermalEmission_2D** target ) const
//{
//	bool ok = true;
//
//	
//	ok = ok && NULL!=target;
//
//	switch(GetThermalEmissionType()){
//		case ThermalType::none:
//			ok = ok && CreateThermalEmissionTable_none   ( target ); 
//			break;
//		case ThermalType::dim2:
//			ok = ok && CreateThermalEmissionTable_2DTable( target );
//			break;
//		default:
//			ok = false;
//			break;
//	}
//
//	return ok;
//}

//bool SKTRAN_Specifications_MC::CreateThermalEmissionTable_none( SKTRAN_ThermalEmission_2D** target ) const
//{
//	*target = NULL;
//	return true;
//}

//bool SKTRAN_Specifications_MC::CreateThermalEmissionTable_2DTable( SKTRAN_ThermalEmission_2D** target ) const
//{
//	bool ok = true;
//	double delta;
//
//	SKTRAN_ThermalEmission_2D* table;
//	if(ok){
//		table = new SKTRAN_ThermalEmission_2D();
//	} else{
//        table = nullptr;
//		nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateThermalEmissionTable, Must have valid target.");
//	}
//
//    			
//    size_t numSolarTableAlts = GetNumSolarTableAlts();		
//    double groundShift = GetGroundShiftAlt();
//	std::vector<SKTRAN_Distance> heights;
//	heights.resize( numSolarTableAlts );
//	ok = ok && 1<heights.size();
//	if(ok){
//		delta = GetSolarTableAltDelta();
//		for( size_t altidx = 0; altidx < heights.size(); altidx++ )
//		{
//			heights[altidx] = groundShift + double(altidx*delta);
//		}
//	}
//
//	SKTRAN_GridDefCosSZA_V21	cosszagrid;
//	SKTRAN_GridDefSLON_V21      slongrid;
//	cosszagrid.AllocateGridArray( m_numSolTableCosSza );
//	slongrid. AllocateGridArray( m_numSolTableCosSza );
//	ok = ok && 0 < cosszagrid.NumAngles();
//	if(ok){
//		delta = 2.0 / double(cosszagrid.NumAngles() - 1);
//		for( size_t szaidx = 0; szaidx < cosszagrid.NumAngles(); szaidx++ )
//		{
//			cosszagrid.AtVar( szaidx ) = -1.0 + szaidx*delta;
//			slongrid.AtVar(szaidx) = 0.0;
//		}
//	}
//
//	ok = ok && table->SetGeometry( heights, cosszagrid, slongrid );
//	if(ok){
//		*target = table;
//		(*target)->AddRef();
//	} else{
//		delete table;
//	}
//
//	return ok;
//}

bool SKTRAN_Specifications_MC::CreateOpticalPropertiesUnitSphere( const SKTRAN_CoordinateTransform_V2& coords )
{
	nxVector referenceHELIO;
	switch( m_optTableType )
	{
		case OptTableType::dim1:
		case OptTableType::dim1_constant:
			referenceHELIO.SetCoords( coords.ReferencePoint(0.0).UnitVector().X(),coords.ReferencePoint(0.0).UnitVector().Y(),coords.ReferencePoint(0.0).UnitVector().Z() );
			m_opticalunitsphere = new SKTRAN_UnitSphere_Dummy( referenceHELIO );
			m_opticalunitsphere->AddRef();
			return true;
			break;
		case OptTableType::dim3_delaunay:
			return CreateDelaunaySphere( coords );
			break;
		default:
			nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateOPticalPropertyTables, invalid type" );
			break;
	}
	return false;
}

bool SKTRAN_Specifications_MC::CreateDelaunaySphere( const SKTRAN_CoordinateTransform_V2& coords )
{
	size_t numcircles = 100;				// number of circles around the tangent point
	size_t pointspercircle = 10;
	double circlesep = 0.1;
	double zenith;						// the angle off of the tangent
	double azimuth;						// angle around the cone
	HELIODETIC_POINT temp;
	HELIODETIC_VECTOR tempv;
	nxVector		tangentpoint;

	size_t numpoints = numcircles * pointspercircle + 1;	// extra 1 is for tangent point

	std::vector< nxVector > unitvecs;
	unitvecs.resize( numpoints );

	bool hascreatedsphere = false;
	SKTRAN_UnitSphere_Delaunay_nonTabledLookup* spherelocal;

	temp = coords.ReferencePoint( 0.0 );
	tempv.SetCoords( temp.Vector().X(), 0, temp.Vector().Z() );
	tangentpoint.SetCoords( tempv.UnitVector().X(),
						   tempv.UnitVector().Y(),
						   tempv.UnitVector().Z() ); // tangent point;
	unitvecs[0] = tangentpoint;

	// The delaunay triangulation code fails frequently, if it does fail nudge the tangent
	// point very slightly off and retry
	size_t numnudges = 0;
	while ( !hascreatedsphere && numnudges < 100 )
	{
		for( size_t circleidx = 0; circleidx < numcircles; circleidx++ )
		{
			for( size_t pointidx = 0; pointidx < pointspercircle; pointidx++ )
			{ 
				zenith = (circleidx+1) * circlesep;
				azimuth = 2*nxmath::Pi * (pointidx + (circleidx % 2)*0.5) / (pointspercircle);
				unitvecs[pointidx + circleidx * pointspercircle + 1] = CalcRotatedVector( unitvecs[0], zenith, azimuth );
			}
		}

		nxVector tangentopposite;
		tangentopposite.SetCoords( -1.0*unitvecs[0].X(),
								   -1.0*unitvecs[0].Y(),
								   -1.0*unitvecs[0].Z() );
		spherelocal = new SKTRAN_UnitSphere_Delaunay_nonTabledLookup;
		hascreatedsphere = spherelocal->CreateTriangulation( &unitvecs[0], numpoints, &tangentopposite );
		if( !hascreatedsphere )
		{
			delete spherelocal;
			// nudge the reference point off slightly
			unitvecs[0] = CalcRotatedVector( tangentpoint, 0.001, numnudges*2*nxmath::Pi/100 );
			++numnudges;
		}
	}
	if( hascreatedsphere )
	{
		m_opticalunitsphere = spherelocal;
		m_opticalunitsphere->AddRef();
	}
	else
	{
		nxLog::Record( NXLOG_WARNING, "Delaunay triangulation failed even after 100 nudges" );
	}

	return true;

}

nxVector SKTRAN_Specifications_MC::CalcRotatedVector( const nxVector& tangent, double zenith, double azimuth ) const
{
	// note that tangent_y is 0 since we force the tangent to be at slon 0
	// so rotate the zenith in 2 dimensions first

	nxVector result;

	result.SetCoords( nxmath::cosd(zenith) * tangent.X() - nxmath::sind(zenith) * tangent.Z(),
					  0,
					  nxmath::sind(zenith) * tangent.X() + nxmath::cosd(zenith) * tangent.Z() );

	// now we need to rotate about the tangent point, but still remembering
	// tangent_y = result_y = 0
	double cosazi = cos(azimuth);
	double sinazi = sin(azimuth);
	double ux = tangent.X();
	double uz = tangent.Z();
	double rx = result.X();
	double rz = result.Z();

	result.SetCoords( ( cosazi + ux*ux*(1-cosazi) ) * rx + (ux*uz*(1-cosazi)) * rz,
		              uz*sinazi*rx + -1.0*ux*sinazi*rz,
					  uz*ux*(1-cosazi)*rx + (cosazi + uz*uz*(1-cosazi))*rz );

	return result;
}

bool SKTRAN_Specifications_MC::CreateOpticalPropertyTables( SKTRAN_TableOpticalProperties_Base** optprop, SKTRAN_TableOpticalProperties_MCBase** mcoptprop, const SKTRAN_CoordinateTransform_V2& coords)
{
	bool	ok = true;
	 
	ok = ok && nullptr!=optprop && nullptr!=mcoptprop;
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateOpticalPropertyTables, User must provide valid table references.");
	if(ok){
		if (m_optTableType == OptTableType::dim1_constant)
		{
			*mcoptprop = new SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant_MC;
		}
		else
		{
			*mcoptprop = new SKTRAN_TableOpticalProperties_3D_UnitSphere_MC;
		}
		CreateOpticalPropertiesUnitSphere( coords );
	}

	ok = ok && nullptr!=*mcoptprop;
	if(!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_Engine_MC_V21::CreateOpticalPropertyTables, Error creating optical properties table");

	// Table contains both SKTRAN and MC parts -- make pointer to the same object to access SKTRAN-type functionality
	if ( m_optTableType == OptTableType::dim1_constant )
	{
		*optprop = dynamic_cast<SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant*>(*mcoptprop);
	}
	else
	{
		*optprop = dynamic_cast<SKTRAN_TableOpticalProperties_3D_UnitSphere*>(*mcoptprop);
	}
	ok = ok && nullptr!=*optprop;
	ok = ok && nullptr!=optprop;
	if (!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_Engine_MC_V21::CreateOpticalPropertyTables, Table must contain basic SKTRAN-type functionality");

	std::unique_ptr< SKTRAN_PolarizationProperties_Base > polobj;
	ok = ok && CreatePolarizationObject ( polobj );
	if(ok) (*optprop)->SetPolarizationProperties( polobj );

	// create inelastic optical property table
	std::shared_ptr<SKTRAN_TableOpticalProperties_Inelastic_Base> inelasticOptTable;
	ok = ok && CreateInelasticTable(inelasticOptTable);

	ok = ok && CreatePolarizationObject(polobj);
	if (ok) inelasticOptTable->SetPolarizationProperties(polobj);

	ok = ok && (*mcoptprop)->SetInelasticProperties(inelasticOptTable);

	if(ok){
		ok = ok && (*optprop)->AddRef();
		if (m_optTableType == OptTableType::dim1_constant)
		{
			ok = ok && dynamic_cast<SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant*>(*optprop)->ConfigureGeometry(this);	// TODO this absolutely needs to be fixed
		}
		else
		{
			ok = ok && dynamic_cast<SKTRAN_TableOpticalProperties_3D_UnitSphere*>(*optprop)->ConfigureGeometry(this);	// TODO this absolutely needs to be fixed
		}
	}

	if (!ok) nxLog::Record(NXLOG_WARNING,"SKTRAN_Engine_MC_V21::CreateOpticalPropertyTables, Error configuring optical properties table");

	//// create amf optical table
	//*amfoptprop = new SKTRAN_TableOpticalProperties_1D_Height_V3;
	//std::unique_ptr< SKTRAN_PolarizationProperties_Base > amfpolobj;
	//ok = ok && CreatePolarizationObject(amfpolobj);
	//if (ok) (*amfoptprop)->SetPolarizationProperties(amfpolobj);
	//ok = ok && (*amfoptprop)->ConfigureGeometry(*GetScatterAngleGrid(), *GetOpticalPropRadii());



	return ok;

}

bool SKTRAN_Specifications_MC::CreateAveragingKernel(SKTRAN_PhotonLog_Base** target) const 
{
	bool ok = true;

	ok = ok && nullptr!=target;

	switch(m_photonLogType)
	{
		case LogType::none:
			(*target) = new SKTRAN_PhotonLog_Null;
			break;
		case LogType::obsPlane:
			(*target) = new SKTRAN_PhotonLog_AveKernel;
			ok = ok && NULL!=target;
			break;
		case LogType::stDev:
			{
			SKTRAN_PhotonLog_StDev* pl = new SKTRAN_PhotonLog_StDev;
			if(nullptr!=pl) pl->ConfigureHistoryIntervals ( 100, 100 );
			(*target) = pl;
			}
			break;
		case LogType::radAlongLOS:
			(*target) = new SKTRAN_PhotonLog_RadianceOnLos;
			ok = ok && nullptr!=target;
			break;
		case LogType::photAlongLOS:
			(*target) = new SKTRAN_PhotonLog_PhotonsOnLos;
			ok = ok && nullptr!=target;
			break;
		case LogType::scatPtOnLOS:
			(*target) = new SKTRAN_PhotonLog_ScatterPtOnLos;
			ok = ok && nullptr!=target;
			break;
		default:
			nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreateAveragingKernel, Averaging kernel type not recognized.");
			ok = false;
			break;
	}
	ok = ok && nullptr!=*target;

	return ok;
}

bool SKTRAN_Specifications_MC::CreateScatterOperator( std::shared_ptr<SKTRAN_MCScatterOperator_Base>& target ) const
{
	bool ok = true;

	SKTRAN_HPFOSet* hpfos;

    std::shared_ptr<SKTRAN_MCScatterOperator_Base> scatterOp;

	switch (m_scatterType)
	{
	case ScatterType::elastic:
		switch (m_polType)
		{
		case PolType::none:
			scatterOp = std::shared_ptr<SKTRAN_MCScatterOperator_Base>(new SKTRAN_MCScatterOperator_Scalar);
			break;
		case PolType::pol:
			scatterOp = std::shared_ptr<SKTRAN_MCScatterOperator_Base>(new SKTRAN_MCScatterOperator_Polarized);
			break;
		case PolType::pv1:
			scatterOp = std::shared_ptr<SKTRAN_MCScatterOperator_Base>(new SKTRAN_MCScatterOperator_PseudoPolarized);
			break;
		default:
			nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateScatterOperator, Optical properties table type not recognized.");
			break;
		} break;
	case ScatterType::both:				// intentional fallthrough
	case ScatterType::manualBoth:		// intentional fallthrough
	case ScatterType::manualInelastic:	// intentional fallthrough
	case ScatterType::inelastic:
		switch (m_polType)
		{
		case PolType::none:
			scatterOp = std::shared_ptr<SKTRAN_MCScatterOperator_Base>(new SKTRAN_MCScatterOperator_ScalarInelastic);
			break;
		case PolType::pol:
			scatterOp = nullptr;
			nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateScatterOperator, Polarized inelastic scattering has not been implemented.");
			break;
		case PolType::pv1:
			scatterOp = nullptr;
			nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateScatterOperator, Polarized inelastic scattering has not been implemented.");
			break;
		default:
			nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateScatterOperator, Optical properties table type not recognized.");
			break;
		} break;
	}
	ok = ok && nullptr!=scatterOp;

	switch(m_symmetryType)
	{
	case SymmetryType::none:
		hpfos = new SKTRAN_HPFOSet_NoSymmetry;
		break;
	case SymmetryType::horiz:
		hpfos = new SKTRAN_HPFOSet_HorizSymmetry;
		break;
	default:
        hpfos = nullptr;
		nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateScatterOperator, Atmospheric symmetry type not recognized.");
		break;
	}
	ok = ok && nullptr!=hpfos;

	ok = ok && scatterOp->SetHPFlipOpSymmetry( hpfos );
	if(ok){
        target = scatterOp;
	} else{
		delete hpfos;
	}

	return ok;
}

bool SKTRAN_Specifications_MC::AddInfoToGenericRayTracer( SKTRAN_RayTracer_Straight_Generic& raytracer, const SKTRAN_CoordinateTransform_V2& coords, const SKTRAN_GridDefRayTracingShells_V21* shells) const
{
	bool ok = true;

	// first add the raytracing heights
	// TODO: different for solar/los

	for( size_t idx = 0; idx < shells->NumShells(); idx++ )
	{
		raytracer.AddGeometryObject( std::unique_ptr<SKTRAN_GeometryObject_Sphere>(new SKTRAN_GeometryObject_Sphere( shells->At( idx ) + coords.AltitudeToRadius(0.0) ) ) );
	}

	// now add the ground/atmo
	raytracer.SetEarthRadius( coords.AltitudeToRadius(0.0) + GetGroundShiftAlt() );
	raytracer.SetUpperAtmoRadius( coords.AltitudeToRadius(0.0) + GetTOAHeight() );
	if( m_optTableType == OptTableType::dim3_delaunay )
	{
		// TODO: make these member varibles
		size_t numcones = 100;
		size_t profilepercone = 10;
		double coneanglesep = 0.1;
		// First add the cones
		// series of cones centered on the reference point, spaced evenly with m_coneanglesep 
		HELIODETIC_VECTOR refhelio = coords.ReferencePoint(0.0).Vector();
		nxVector coneunitvec = nxVector(refhelio.X(), refhelio.Y(), refhelio.Z()).UnitVector();

		for( size_t idx = 1; idx < numcones; idx++ )
		{
			raytracer.AddGeometryObject( std::unique_ptr<SKTRAN_GeometryObject_Cone> ( new SKTRAN_GeometryObject_Cone( coneunitvec, idx*coneanglesep ) ) );
		}

		// Now add the planes
		// Define each plane by two vectors, the reference, and a point on a circle 10 degrees off
		// since we stagger the spacing every circle we need twice as many planes
		double zenith = 10 * nxmath::Pi/180;
		double azimuth;
		nxVector secondvec;
		nxVector normal;
		for( size_t idx = 0; idx < profilepercone*2; idx++ )
		{
			azimuth = 2*nxmath::Pi * (idx) / (profilepercone*2);
			secondvec = CalcRotatedVector( coneunitvec, zenith, azimuth );
			normal = coneunitvec.Cross( secondvec ).UnitVector();
			raytracer.AddGeometryObject( std::unique_ptr<SKTRAN_GeometryObject_Plane> ( new SKTRAN_GeometryObject_Plane( normal ) ) );
		}
	}

	return ok;
}

			
bool SKTRAN_Specifications_MC::GetProfileAlts(SKTRAN_Specifications_MC::ProfileType profileType, std::vector<double>& profile, bool& uniform) const
{	
	bool ok = true;
	bool ok_nonessential = true;

	bool manual;
	bool mustspan = true, span = true;
	size_t numShells;
	double shellSpacing;

	const double toaHeight = GetTOAHeight();
	const double groundShiftHeight = GetGroundShiftAlt();

	switch (profileType)
	{
		case ProfileType::rayTracing:
			if (m_manualRayTracingHeights.size() > 0) {
				manual = true;
				profile = m_manualRayTracingHeights;
				numShells = profile.size();
			}
			else {
				manual = false;
				numShells = m_numRayTracingAlts;
				shellSpacing = (toaHeight - groundShiftHeight) / (numShells - 1);
			}
			break;
		case ProfileType::opticalProperties:
			if (m_manualopticalheights.size() > 0)
			{
				manual = true;
				profile = m_manualopticalheights;
				numShells = profile.size();
			}
			else {
				manual = false;
				numShells = m_numOptPropAlts;
				shellSpacing = (toaHeight - groundShiftHeight) / (numShells - 1);
			}
			break;
		case ProfileType::airMassFactor:
			manual = true;
			profile = m_manualAMFHeights;
			numShells = profile.size();
			mustspan = false;
			break;
		case ProfileType::solarTable:
			if (m_manualSolarTableHeights.size() > 0) {
				manual = true;
				profile = m_manualSolarTableHeights;
				numShells = profile.size();
			}
			else {
				manual = false;
				shellSpacing = m_solarTableAltDelta;
				numShells = (size_t)std::ceil((toaHeight - groundShiftHeight) / shellSpacing) + 1;
			}
			break;
		default:
			nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::GetProfileAlts, unknown profile type.");
	}

	if (manual) {
		ok = ok && numShells > 1;
		if (!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::GetProfileAlts, manual grid has length less than 2.");

		// check for increasing order
		for (size_t shidx = 0; shidx < numShells - 1; shidx++) ok = ok && profile[shidx + 1] - profile[shidx] > 0.0;
		if (!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::GetProfileAlts, manual grid is not in increasing order.");

		// check that the manual grid spans the whole atmosphere
		span = span && *std::max_element(profile.begin(), profile.end()) == toaHeight;
		span = span && *std::min_element(profile.begin(), profile.end()) == groundShiftHeight;
		ok = ok && !(mustspan && !span);
		if (!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::GetProfileAlts, minimum and maximum altitudes in manual grid do not match surface and TOA altitudes.");
		

		// check for uniform spacing
		uniform = true;
		double delta = profile[1] - profile[0];
		for (size_t shidx = 0; shidx < numShells - 1; shidx++) uniform = uniform && profile[shidx + 1] - profile[shidx] == delta;

		if (profileType == ProfileType::airMassFactor && !span) uniform = false;  // if AMF shells do not span the atmosphere, extra shells will be added which will make them nonuniform
	}
	else {
		ok = ok && numShells > 1;
		if (!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::GetProfileAlts, there are fewer than 2 profile alts specified.");

		profile.resize(numShells);
		for (size_t shidx = 0; shidx < numShells; shidx++) profile[shidx] = groundShiftHeight + shidx * shellSpacing;

		uniform = true;
	}    	
	
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Specifications_MC::SetRayTracers		 2014- 11- 19*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Specifications_MC::SetRayTracers(SKTRAN_Engine_MC_V21* engine, 	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords ) const
{
	bool ok = true;
	std::unique_ptr< SKTRAN_RayFactory<	SKTRAN_RayOptical_Straight,
									    SKTRAN_RayTracer_Shells,
										SKTRAN_RayStorage_Straight_MC> >  rayfactorystraight;

	std::unique_ptr< SKTRAN_RayFactory< SKTRAN_RayOptical_Curved,
										SKTRAN_RayTracer_Curved_Shells,
										SKTRAN_RayStorage_CurvedPiecewise_MC> >	rayfactorycurved;

	std::unique_ptr< SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
										SKTRAN_RayTracer_Straight_Generic,
										SKTRAN_RayStorage_Straight_MC> > rayfactorygeneric;

	switch(m_solarRayTracerType)
	{
		case RayTracerType::shell:
			rayfactorystraight.reset( new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Shells, SKTRAN_RayStorage_Straight_MC>(coords));
			rayfactorystraight->RayTracer()->Initialize( m_raytracingshells );
			engine->SetRayFactory_SOLAR( std::move(rayfactorystraight));
			break;
		case RayTracerType::curved:
			nxLog::Record(NXLOG_WARNING,"SKTRAN_Specifications_MC::SetRayTracers, **** TODO **** Curved solar rays need proper refractive index definition");
			rayfactorycurved.reset ( new SKTRAN_RayFactory< SKTRAN_RayOptical_Curved, SKTRAN_RayTracer_Curved_Shells, SKTRAN_RayStorage_CurvedPiecewise_MC>(coords));
			rayfactorycurved->RayTracer()->Initialize( m_raytracingshells, std::move(std::unique_ptr<skRTRefractiveIndex_Profile>(new skRTRefractiveIndex_Profile()))  );
			engine->SetRayFactory_SOLAR( std::move(rayfactorycurved));
			break;
		case RayTracerType::generic:
			rayfactorygeneric.reset( new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Straight_Generic, SKTRAN_RayStorage_Straight_MC>(coords));
			ok = ok && AddInfoToGenericRayTracer( *rayfactorygeneric->RayTracer(), *coords, m_raytracingshells.get() );
			engine->SetRayFactory_SOLAR( std::move(rayfactorygeneric));
			break;
		default:
			ok = false;
	}
	


	switch(m_LOSRayTracerType)
	{
		case RayTracerType::shell:
			rayfactorystraight.reset( new SKTRAN_RayFactory<	SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Shells, SKTRAN_RayStorage_Straight_MC>(coords));
			rayfactorystraight->RayTracer()->Initialize( m_raytracingshells);
			engine->SetRayFactory_LOS( std::move(rayfactorystraight));
			break;
		case RayTracerType::curved:
			nxLog::Record(NXLOG_WARNING,"SKTRAN_Specifications_MC::SetRayTracers, **** TODO **** Curved LOS rays need proper refractive index definition");
			rayfactorycurved.reset( new SKTRAN_RayFactory<	SKTRAN_RayOptical_Curved, SKTRAN_RayTracer_Curved_Shells, SKTRAN_RayStorage_CurvedPiecewise_MC>(coords));
			rayfactorycurved->RayTracer()->Initialize( m_raytracingshells, std::move(std::unique_ptr<skRTRefractiveIndex_Profile>(new skRTRefractiveIndex_Profile()))  );
			engine->SetRayFactory_LOS( std::move(rayfactorycurved));
			break;
		case RayTracerType::generic:
			rayfactorygeneric.reset( new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Straight_Generic, SKTRAN_RayStorage_Straight_MC>(coords));
			ok = ok && AddInfoToGenericRayTracer( *rayfactorygeneric->RayTracer(), *coords, m_raytracingshells.get() );
			engine->SetRayFactory_LOS( std::move(rayfactorygeneric));
			break;
		default:
			ok = false;
	}

	switch(m_MSRayTracerType)
	{
		case RayTracerType::shell:
			rayfactorystraight.reset( new SKTRAN_RayFactory<	SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Shells, SKTRAN_RayStorage_Straight_MC>(coords));
			rayfactorystraight->RayTracer()->Initialize( m_raytracingshells);
			engine->SetRayFactory_SECONDARY( std::move(rayfactorystraight));
			break;
		case RayTracerType::curved:
			nxLog::Record(NXLOG_WARNING,"SKTRAN_Specifications_MC::SetRayTracers, **** TODO **** Curved LOS rays need proper refractive index definition");
			rayfactorycurved.reset( new SKTRAN_RayFactory<	SKTRAN_RayOptical_Curved, SKTRAN_RayTracer_Curved_Shells, SKTRAN_RayStorage_CurvedPiecewise_MC>(coords));
			rayfactorycurved->RayTracer()->Initialize( m_raytracingshells, std::move(std::unique_ptr<skRTRefractiveIndex_Profile>(new skRTRefractiveIndex_Profile()))  );
			engine->SetRayFactory_SECONDARY( std::move(rayfactorycurved));
			break;
		case RayTracerType::generic:
			rayfactorygeneric.reset( new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Straight_Generic, SKTRAN_RayStorage_Straight_MC>(coords));
			ok = ok && AddInfoToGenericRayTracer( *rayfactorygeneric->RayTracer(), *coords, m_raytracingshells.get() );
			engine->SetRayFactory_SECONDARY( std::move(rayfactorygeneric));
			break;
		default:
			ok = false;
	}

	return ok;
}


bool SKTRAN_Specifications_MC::CreateOpticalPropsIntegrator( SKTRAN_OpticalPropertiesIntegrator_Base** target ) const
{
	bool ok = true;

	ok = ok && NULL!=target;

	switch(m_optPropsIntType)
	{
		case OptPropIntType::straight:
			ok = ok && CreateOpticalPropsIntegrator_Straight( target );
			break;
		case OptPropIntType::adaptive:
			ok = ok && CreateOpticalPropsIntegrator_Adaptive( target );
			break;
		case OptPropIntType::constant:
			ok = ok && CreateOpticalPropsIntegrator_Constant( target );
			break;
		default:
			ok = false;
			break;
	}

	if(ok) (*target)->AddRef();

	return ok;
}

bool SKTRAN_Specifications_MC::CreateOpticalPropsIntegrator_Straight( SKTRAN_OpticalPropertiesIntegrator_Base** target ) const
{
	bool ok = true;

	if (m_wavelengthType == WavelengthType::single && m_scatterType == ScatterType::elastic)
	{
		*target = new SKTRAN_OpticalPropertiesIntegrator_Straight(); // single wavelength for whole calculation
	}
	else
	{
		*target = new SKTRAN_OpticalPropertiesIntegrator_Straight_MC(); // interpolation for inelastic scatters and/or wavelength switching for simultaneous wavelengths
	}

	return ok && NULL!= target;
}

bool SKTRAN_Specifications_MC::CreateOpticalPropsIntegrator_Adaptive( SKTRAN_OpticalPropertiesIntegrator_Base** target ) const
{
	bool ok = true;
	
	if (m_wavelengthType == WavelengthType::single && m_scatterType == ScatterType::elastic)
	{
		SKTRAN_OpticalPropertiesIntegrator_Adaptive* opi = new SKTRAN_OpticalPropertiesIntegrator_Adaptive();
		ok = ok && NULL != opi;
		if (ok) {
			opi->SetMaxOpticalDepthOfCell(GetAdaptOptDepthMax());
			opi->SetMaxExtinctionGradientOfCell(GetAdaptOptDepthMinRatio());
			*target = opi;
		}
	}
	else
	{
		SKTRAN_OpticalPropertiesIntegrator_Adaptive_MC* opi = new SKTRAN_OpticalPropertiesIntegrator_Adaptive_MC();
		ok = ok && NULL != opi;
		if (ok) {
			opi->SetMaxOpticalDepthOfCell(GetAdaptOptDepthMax());
			opi->SetMaxExtinctionGradientOfCell(GetAdaptOptDepthMinRatio());
			*target = opi;
		}
	}

	return ok;
}

bool SKTRAN_Specifications_MC::CreateOpticalPropsIntegrator_Constant(SKTRAN_OpticalPropertiesIntegrator_Base** target) const
{
	bool ok = true;

	if (m_wavelengthType == WavelengthType::single && m_scatterType == ScatterType::elastic)
	{
		*target = new SKTRAN_OpticalPropertiesIntegrator_ConstantLayers; // single wavelength for whole calculation
	}
	else
	{
		*target = new SKTRAN_OpticalPropertiesIntegrator_ConstantLayers_MC; // interpolation for inelastic scatters and/or wavelength switching for simultaneous wavelengths
	}

	return ok && NULL != target;
}


bool SKTRAN_Specifications_MC::CreateRayTracers( SKTRAN_Engine_MC_V21* engine ) const
{
	bool ok = true;

	// Create all ray tracers
//	if( rtType_shell==m_solarRayTracerType || rtType_shell==m_LOSRayTracerType || rtType_shell==m_MSRayTracerType ){
//		ok = ok && CreateRayTracer_Shells( &engine->m_raytracer_shells );
//	}
//	if( rtType_curved==m_solarRayTracerType || rtType_curved==m_LOSRayTracerType || rtType_curved==m_MSRayTracerType ){
//		ok = ok && CreateRayTracer_Shells_Curved( &engine->m_raytracer_shells_curved );
//	}

	return ok;
}

//bool SKTRAN_Specifications_MC::CreateRayTracer_Shells( SKTRAN_RayTracer_Shells** rt ) const
//{
//	bool ok						= true;
//	SKTRAN_RayTracer_Shells*	raytracer_shells;
//	raytracer_shells			= new SKTRAN_RayTracer_Shells();
//
//	ok = ok && rt!=NULL;
//	ok = ok && NULL!=raytracer_shells;
//	ok = ok && raytracer_shells->Initialize( GetRayTracingShells() );
//	if(ok){
//		*rt = raytracer_shells;
//		(*rt)->AddRef();
//	}
//	return ok;
//}

//bool SKTRAN_Specifications_MC::CreateRayTracer_Shells_Curved( SKTRAN_RayTracer_Shells_Curved** rt ) const
//{
//
//	bool							ok = true;
//	SKTRAN_RayTracer_Shells_Curved* raytracer_shells_curved;
//	skRTRefractiveIndex_Profile*	n_profile;
//
//	raytracer_shells_curved			= new SKTRAN_RayTracer_Shells_Curved();
//	n_profile						= new skRTRefractiveIndex_Profile();
//
//	ok = ok && NULL!=rt;
//	ok = ok && NULL!=raytracer_shells_curved;
//	ok = ok && NULL!=n_profile;
//
//	ok = ok && raytracer_shells_curved->Initialize( GetRayTracingShells(), n_profile);
//	if(ok){
//		raytracer_shells_curved->SetRayStepSize( GetCurvedRayStepSize() );
//		(*rt) = raytracer_shells_curved;
//		(*rt)->AddRef();
//	} 
//
//	return ok;
//}

bool SKTRAN_Specifications_MC::CreatePolarizationObject ( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& target ) const
{
	bool ok = true;
    
	if(ok){
		switch(m_polType){
			case PolType::none:
				target = std::unique_ptr< SKTRAN_PolarizationProperties_Base > ( new SKTRAN_PolarizationProperties_NoPolarization );
				break;
				
			case PolType::pv1:
			case PolType::pol:
				//if( m_atmosphereHasDelta ){
				//	(*target) = new SKTRAN_PolarizationProperties_Polarized_Eddington;
				//} else{
					//nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::CreatePolarizationObject, atmoHasDelta needs to be used to select acceleration methods (VROOM).");
					target = std::unique_ptr< SKTRAN_PolarizationProperties_Base > ( new SKTRAN_PolarizationProperties_Polarized );
				//}
				break;
			default:
				nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreatePolarizationObject, Polarization type not recognized.");
		}
	}
	ok = ok && nullptr!=target;
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreatePolarizationObject, Could not make polarization object.");
	
	return ok;
}

bool SKTRAN_Specifications_MC::CreateInelasticTable(std::shared_ptr< SKTRAN_TableOpticalProperties_Inelastic_Base >& target) const
{
	bool ok = true;

	switch (m_scatterType)
	{
		case ScatterType::elastic:
			target.reset( new SKTRAN_TableOpticalProperties_Inelastic_DoNothing );
			break;
		case ScatterType::both:				// intentional fallthrough
		case ScatterType::manualBoth:		// intentional fallthrough
		case ScatterType::manualInelastic:	// intentional fallthrough
		case ScatterType::inelastic:
			target.reset( new SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere );
			break;
		default:
			ok = false;
			break;
	}
	ok = ok && target != nullptr;

	return ok;
}
bool SKTRAN_Specifications_MC::CreateAirMassFactorCalculator(std::unique_ptr<SKTRAN_MCAirMassFactorCalculator_Base>& target) const
{
	bool ok = true;

	switch (m_secondary)
	{
	case SecondaryOutput::lengthAMF:
		target.reset(new SKTRAN_MCAirMassFactorCalculator_Length);
		ok = ok && NULL != target;
		break;
	case SecondaryOutput::opticalDepthAMF:
		target.reset(new SKTRAN_MCAirMassFactorCalculator_OpticalDepth);
		ok = ok && NULL != target;
		ok = ok && target->SetSpecies(m_amfSpeciesHandle);
		break;
	default:
		target.reset(new SKTRAN_MCAirMassFactorCalculator_DoNothing);
		ok = ok && NULL != target;
		break;
	}

	//target.reset(new SKTRAN_MCAirMassFactorCalculator_Base);
	//ok = NULL != target;

	//bool calcAMF = GetAMFType() != SKTRAN_Specifications_MC::AMFType::none;
	//bool calcSpeciesAMF = GetAMFType() == SKTRAN_Specifications_MC::AMFType::species;
	//bool calcGeometricAMF = GetAMFType() == SKTRAN_Specifications_MC::AMFType::geometric;
	//if (ok) target->SetType(calcAMF, calcSpeciesAMF, calcGeometricAMF);

	return ok;
}

bool SKTRAN_Specifications_MC::CreateAirMassFactorOpticalPropertiesTable( SKTRAN_TableOpticalProperties_Base* optprop, SKTRAN_TableOpticalProperties_Base** amfoptprop) const
{
	// I'm passing in the regular optical property table because it may be possible to use it as the amf optical property table in the future
	bool ok = true;

	switch (m_secondary)
	{
		case SecondaryOutput::lengthAMF:
		{
			if (*amfoptprop != nullptr) (*amfoptprop)->Release(); *amfoptprop = nullptr;
			break;
		}
		case SecondaryOutput::opticalDepthAMF:
		{
			std::unique_ptr< SKTRAN_PolarizationProperties_Base > amfpolobj;
			SKTRAN_TableOpticalProperties_3D_UnitSphere* optprop3d;

			switch (m_optTableType)
			{
				case OptTableType::dim1:
				case OptTableType::dim2_plane:
				case OptTableType::dim3_delaunay:
				{
					optprop3d = new SKTRAN_TableOpticalProperties_3D_UnitSphere;
					break;
				}
				case OptTableType::dim1_constant:
				{
					optprop3d = new SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant;		
					break;
				}
				default:
				{
					if (*amfoptprop != nullptr) (*amfoptprop)->Release(); *amfoptprop = nullptr;
					optprop3d = nullptr;
					break;
				}
			}

			ok = ok && optprop3d != nullptr;
			ok = ok && CreatePolarizationObject(amfpolobj);
			if (ok) optprop3d->SetPolarizationProperties(amfpolobj);
			ok = ok && optprop3d->SetAltitudes(*GetOpticalPropRadii());
			ok = ok && optprop3d->SetScatterGrid(*GetScatterAngleGrid());
			ok = ok && optprop3d->SetUnitSphere(*GetOpticalUnitSphere());
			ok = ok && optprop3d->SetWavelengthGrid(*GetWavelengthGridOpticalProperties());
			ok = ok && optprop3d->ConfigureGeometry(this);
			if (ok) *amfoptprop = optprop3d;

			if (ok)
			{
				(*amfoptprop)->AddRef();
			}
			else
			{
				if (*amfoptprop != nullptr) (*amfoptprop)->Release(); *amfoptprop = nullptr;
				ok = false;
				nxLog::Record(NXLOG_ERROR, "SKTRAN_Specifications_MC::CreateAirMassFactorOpticalPropertiesTable, problem creating AMF optical properties table.");
			}
			break;
		}
		default:
		{
			if (*amfoptprop != nullptr) (*amfoptprop)->Release(); *amfoptprop = nullptr;
			break;
		}
	}
	return ok;
}

bool SKTRAN_Specifications_MC::CreateAirMassFactorOpticalPropsIntegrator( SKTRAN_OpticalPropertiesIntegrator_Base* integrator, SKTRAN_TableOpticalProperties_Base* amfoptprop, SKTRAN_OpticalPropertiesIntegrator_Base** amfintegrator ) const
{
	// I'm passing in the regular optical property integrator because it may be possible to use it as the amf optical property integrator in the future
	bool ok = true;

	switch(m_secondary)
	{
	case SecondaryOutput::lengthAMF:
		if (*amfintegrator != nullptr) (*amfintegrator)->Release(); *amfintegrator = nullptr;
		break;
	case SecondaryOutput::opticalDepthAMF:
		CreateOpticalPropsIntegrator(amfintegrator);
		(*amfintegrator)->SetOpticalProps(amfoptprop);
		break;
	default:
		if (*amfintegrator != nullptr) (*amfintegrator)->Release(); *amfintegrator = nullptr;
		break;
	}
	return ok;
}

bool SKTRAN_Specifications_MC::SetAirMassFactorRayTracers(std::unique_ptr<SKTRAN_MCAirMassFactorCalculator_Base>& amfcalc, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords) const
{
	// this is essentially copied from SetRayTracers, there is likely a way to do this with less copy and paste but this is it for now

	bool ok = true;

	if (m_secondary != SecondaryOutput::lengthAMF && m_secondary != SecondaryOutput::opticalDepthAMF) return ok;
	
	std::unique_ptr< SKTRAN_RayFactory<	SKTRAN_RayOptical_Straight,
		SKTRAN_RayTracer_Shells,
		SKTRAN_RayStorage_Straight_MC> >  rayfactorystraight;

	std::unique_ptr< SKTRAN_RayFactory< SKTRAN_RayOptical_Curved,
		SKTRAN_RayTracer_Curved_Shells,
		SKTRAN_RayStorage_CurvedPiecewise_MC> >	rayfactorycurved;

	std::unique_ptr< SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
		SKTRAN_RayTracer_Straight_Generic,
		SKTRAN_RayStorage_Straight_MC> > rayfactorygeneric;

	switch (m_solarRayTracerType)
	{
	case RayTracerType::shell:
		rayfactorystraight.reset(new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Shells, SKTRAN_RayStorage_Straight_MC>(coords));
		rayfactorystraight->RayTracer()->Initialize(m_amfshells);
		amfcalc->SetRayFactory_SOLAR(std::move(rayfactorystraight), false);
		break;
	case RayTracerType::curved:
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetRayTracers, **** TODO **** Curved solar rays need proper refractive index definition");
		rayfactorycurved.reset(new SKTRAN_RayFactory< SKTRAN_RayOptical_Curved, SKTRAN_RayTracer_Curved_Shells, SKTRAN_RayStorage_CurvedPiecewise_MC>(coords));
		rayfactorycurved->RayTracer()->Initialize(m_amfshells, std::move(std::unique_ptr<skRTRefractiveIndex_Profile>(new skRTRefractiveIndex_Profile())));
		amfcalc->SetRayFactory_SOLAR(std::move(rayfactorycurved), true);
		break;
	case RayTracerType::generic:
		rayfactorygeneric.reset(new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Straight_Generic, SKTRAN_RayStorage_Straight_MC>(coords));
		ok = ok && AddInfoToGenericRayTracer(*rayfactorygeneric->RayTracer(), *coords, m_amfshells.get() );
		amfcalc->SetRayFactory_SOLAR(std::move(rayfactorygeneric), false);
		break;
	default:
		ok = false;
	}

	switch (m_LOSRayTracerType)
	{
	case RayTracerType::shell:
		rayfactorystraight.reset(new SKTRAN_RayFactory<	SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Shells, SKTRAN_RayStorage_Straight_MC>(coords));
		rayfactorystraight->RayTracer()->Initialize(m_amfshells);
		amfcalc->SetRayFactory_LOS(std::move(rayfactorystraight), false);
		break;
	case RayTracerType::curved:
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetRayTracers, **** TODO **** Curved LOS rays need proper refractive index definition");
		rayfactorycurved.reset(new SKTRAN_RayFactory<	SKTRAN_RayOptical_Curved, SKTRAN_RayTracer_Curved_Shells, SKTRAN_RayStorage_CurvedPiecewise_MC>(coords));
		rayfactorycurved->RayTracer()->Initialize(m_amfshells, std::move(std::unique_ptr<skRTRefractiveIndex_Profile>(new skRTRefractiveIndex_Profile())));
		amfcalc->SetRayFactory_LOS(std::move(rayfactorycurved), true);
		break;
	case RayTracerType::generic:
		rayfactorygeneric.reset(new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Straight_Generic, SKTRAN_RayStorage_Straight_MC>(coords));
		ok = ok && AddInfoToGenericRayTracer(*rayfactorygeneric->RayTracer(), *coords, m_amfshells.get());
		amfcalc->SetRayFactory_LOS(std::move(rayfactorygeneric), false);
		break;
	default:
		ok = false;
	}

	switch (m_MSRayTracerType)
	{
	case RayTracerType::shell:
		rayfactorystraight.reset(new SKTRAN_RayFactory<	SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Shells, SKTRAN_RayStorage_Straight_MC>(coords));
		rayfactorystraight->RayTracer()->Initialize(m_amfshells);
		amfcalc->SetRayFactory_SECONDARY(std::move(rayfactorystraight), false);
		break;
	case RayTracerType::curved:
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Specifications_MC::SetRayTracers, **** TODO **** Curved LOS rays need proper refractive index definition");
		rayfactorycurved.reset(new SKTRAN_RayFactory<	SKTRAN_RayOptical_Curved, SKTRAN_RayTracer_Curved_Shells, SKTRAN_RayStorage_CurvedPiecewise_MC>(coords));
		rayfactorycurved->RayTracer()->Initialize(m_amfshells, std::move(std::unique_ptr<skRTRefractiveIndex_Profile>(new skRTRefractiveIndex_Profile())));
		amfcalc->SetRayFactory_SECONDARY(std::move(rayfactorycurved), true);
		break;
	case RayTracerType::generic:
		rayfactorygeneric.reset(new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight, SKTRAN_RayTracer_Straight_Generic, SKTRAN_RayStorage_Straight_MC>(coords));
		ok = ok && AddInfoToGenericRayTracer(*rayfactorygeneric->RayTracer(), *coords, m_amfshells.get());
		amfcalc->SetRayFactory_SECONDARY(std::move(rayfactorygeneric), false);
		break;
	default:
		ok = false;
	}

	return ok;
	
}


bool SKTRAN_Specifications_MC::ConfigureSimultaneousWavelengths(SKTRAN_SimultaneousWavelengthManager& simwl) const
{
	bool ok = true;
	std::vector<double> empty;

	switch (m_wavelengthType)
	{
		case WavelengthType::single: 
			simwl.Configure(empty, 0.0, 0);
			break;
		case WavelengthType::simultaneous:
			simwl.Configure(m_radiancewavelengths, m_primaryWavelength, m_primaryWavelengthIndex);
			break;
		default:
			ok = false;
	}
	return ok;
}

bool SKTRAN_Specifications_MC::CreateOptimalScatterSequenceManager(std::unique_ptr<SKTRAN_OptimalScatterSequenceManager_Base>& target)
{
	bool ok = true;			
	

	switch (m_scatterType)
	{
		case ScatterType::elastic:
		{
			SKTRAN_OptimalScatterSequenceManager_Uniform* sequenceManager = new SKTRAN_OptimalScatterSequenceManager_Uniform;
			const std::vector<double> minFrac = GetMinFractionHigherOrder();
			ok = ok && sequenceManager->SetNumDistinctOrders(8);
			ok = ok && sequenceManager->SetMinFractionHigherOrder(minFrac.front());
			ok = ok && sequenceManager->SetMinNumSamplesHigherOrder(100);
			ok = ok && sequenceManager->ConfigureStatisticsExport(m_exportFileName);
			target.reset(sequenceManager);
			break;
		}
		case ScatterType::both:
		{
			SKTRAN_OptimalScatterSequenceManager_UniformSecondary* sequenceManager;
			switch (m_secondary)
			{
				case SecondaryOutput::ringSpectrum:
					sequenceManager = new SKTRAN_OptimalScatterSequenceManager_UniformRing;
					break;
				case SecondaryOutput::fillingInParameter:
					sequenceManager = new SKTRAN_OptimalScatterSequenceManager_UniformFillingIn;
					break;
				case SecondaryOutput::elasticRaman:
					sequenceManager = new SKTRAN_OptimalScatterSequenceManager_UniformElastic;
					break;
				default: // no point calculating both types of scatter if one of these secondary outputs is not requested
					ok = false;
					break;
			}
			
			const std::vector<double> minFrac = GetMinFractionHigherOrder();
			ok = ok && sequenceManager->SetNumDistinctOrders(8);
			ok = ok && sequenceManager->SetMinFractionHigherOrder(minFrac.front());
			ok = ok && sequenceManager->SetMinNumSamplesHigherOrder(100);
			ok = ok && sequenceManager->ConfigureStatisticsExport(m_exportFileName);
			target.reset(sequenceManager);
			break;
		}
		case ScatterType::inelastic:
		{
			SKTRAN_OptimalScatterSequenceManager_Uniform* sequenceManager = new SKTRAN_OptimalScatterSequenceManager_Uniform;
			const std::vector<double> minFrac = GetMinFractionHigherOrder();
			ok = ok && sequenceManager->SetNumDistinctOrders(8);
			ok = ok && sequenceManager->SetMinFractionHigherOrder(minFrac.front());
			ok = ok && sequenceManager->SetMinNumSamplesHigherOrder(100);
			ok = ok && sequenceManager->ConfigureStatisticsExport(m_exportFileName);
			target.reset(sequenceManager);
			break;
		}
		case ScatterType::manualBoth:
		{
			SKTRAN_OptimalScatterSequenceManager_OptimizedSecondary* sequenceManager;
			switch (m_secondary)
			{
			case SecondaryOutput::ringSpectrum:
				sequenceManager = new SKTRAN_OptimalScatterSequenceManager_OptimizedRing;
				break;
			case SecondaryOutput::fillingInParameter:
				sequenceManager = new SKTRAN_OptimalScatterSequenceManager_OptimizedFillingIn;
				break;
			case SecondaryOutput::elasticRaman:
				sequenceManager = new SKTRAN_OptimalScatterSequenceManager_OptimizedElastic;
				break;
			default: // no point calculating both types of scatter if one of these secondary outputs is not requested
				ok = false;
				break;
			}

			ok = ok && sequenceManager->SetMinFractionHigherOrder(m_minFractionHigherOrder);
			ok = ok && sequenceManager->SetMaxRamanOrders(m_maxRamanOrders);
			ok = ok && sequenceManager->ConfigureStatisticsExport(m_exportFileName);
			target.reset(sequenceManager);
			break;
		}
		case ScatterType::manualInelastic:
		{
			SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic* sequenceManager = new SKTRAN_OptimalScatterSequenceManager_OptimizedInelastic;
			ok = ok && sequenceManager->SetMinFractionHigherOrder(m_minFractionHigherOrder);
			ok = ok && sequenceManager->SetMaxRamanOrders(m_maxRamanOrders);
			ok = ok && sequenceManager->ConfigureStatisticsExport(m_exportFileName);
			target.reset(sequenceManager);
			break;
		}
		default:
			ok = false;
			break;
	}

	ok = ok && target;

	return ok;
}


bool SKTRAN_Specifications_MC::CreatePhotonTemplate(std::unique_ptr<const SKTRAN_MCPhoton_Base>& target) const
{
	bool ok = true;

	switch (m_wavelengthType)
	{
		case WavelengthType::single:
			switch (m_scatterType)
			{
				case ScatterType::elastic:
				{
					SKTRAN_MCPhoton* tmp = new SKTRAN_MCPhoton;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = false;
					target.reset(tmp);
					break;
				}
				case ScatterType::inelastic:
				{
					SKTRAN_MCPhoton_Inelastic* tmp = new SKTRAN_MCPhoton_Inelastic;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = false;
					target.reset(tmp);
					break;
				}
				case ScatterType::manualInelastic:
				{
					SKTRAN_MCPhoton_Inelastic* tmp = new SKTRAN_MCPhoton_Inelastic;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = true;
					target.reset(tmp);
					break;
				}
				case ScatterType::both:
				{
					SKTRAN_MCPhoton_Ring* tmp = new SKTRAN_MCPhoton_Ring;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = false;
					target.reset(tmp);
					break;
				}
				case ScatterType::manualBoth:
				{
					SKTRAN_MCPhoton_Ring* tmp = new SKTRAN_MCPhoton_Ring;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = true;
					target.reset(tmp);
					break;
				}
				default:
					break;
			}
			break;
		case WavelengthType::simultaneous:
			switch (m_scatterType)
			{
				case ScatterType::elastic:
				{
					SKTRAN_MCPhoton_Simultaneous* tmp = new SKTRAN_MCPhoton_Simultaneous;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = false;
					target.reset(tmp);
					break;
				}
				case ScatterType::inelastic:
				{
					SKTRAN_MCPhoton_SimultaneousInelastic* tmp = new SKTRAN_MCPhoton_SimultaneousInelastic;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = false;
					target.reset(tmp);
					break;
				}
				case ScatterType::manualInelastic:
				{
					SKTRAN_MCPhoton_SimultaneousInelastic* tmp = new SKTRAN_MCPhoton_SimultaneousInelastic;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = true;
					target.reset(tmp);
					break;
				}
				case ScatterType::both:
				{
					SKTRAN_MCPhoton_SimultaneousRing* tmp = new SKTRAN_MCPhoton_SimultaneousRing;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = false;
					target.reset(tmp);
					break;
				}				
				case ScatterType::manualBoth:
				{
					SKTRAN_MCPhoton_SimultaneousRing* tmp = new SKTRAN_MCPhoton_SimultaneousRing;
					tmp->SetWavelengths(m_radiancewavelengths);
					tmp->m_manualScatter = true;
					target.reset(tmp);
					break;
				}
				default:
					ok = false;
					break;
			}
			break;
		default:
			ok = false;
			break;
	}

	ok = ok && target;

	return ok;
}



