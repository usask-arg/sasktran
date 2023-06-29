#include <skopticalproperties21.h>


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::skOpticalProperties_AerosolProfile		2008-2-28*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_AerosolProfile::skOpticalProperties_AerosolProfile()
{
	init();
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::init		2008-2-28*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_AerosolProfile::init()
{
	m_lastpt.latitude  = 0.0;
	m_lastpt.longitude = 0.0;
	m_lastpt.heightm   = 20000.0;
	m_lastpt.mjd       = 56000.0;

	m_mieaerosol            = NULL;
	m_distribution          = NULL;
	m_moderadiusclimatology = NULL;
	m_refractiveindex       = NULL;
	m_isdirty               = true;
	m_heightdependent       = true;

}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::~skOpticalProperties_AerosolProfile		2008-2-28*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_AerosolProfile::~skOpticalProperties_AerosolProfile()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::ReleaseResources		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_AerosolProfile::ReleaseResources()
{
	if (m_refractiveindex       != NULL) m_refractiveindex->Release();
	if (m_moderadiusclimatology != NULL) m_moderadiusclimatology->Release();
	if (m_mieaerosol            != NULL) m_mieaerosol->Release();
	if (m_distribution          != NULL) m_distribution->Release();

	m_moderadiusclimatology  = NULL;
	m_refractiveindex        = NULL;
	m_distribution           = NULL;
	m_mieaerosol             = NULL;
	m_isdirty                = true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::CheckDirtyAndUpdate		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::CheckDirtyAndUpdate()
{
	bool	ok;

	ok = (!m_isdirty);
	if (!ok)
	{
		ok =     (m_mieaerosol      != NULL)
		      && (m_distribution    != NULL)
			  && (m_refractiveindex != NULL);

		ok = ok && m_mieaerosol->Set_ParticleDistribution( m_distribution );
		ok = ok && m_mieaerosol->Set_RefractiveIndex	 ( m_refractiveindex );
		m_isdirty = !ok;
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfile::CheckDirtyAndUpdate, error during update, probably means the object is not properly configured");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::DeepCopy		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_AerosolProfile::DeepCopy( const skOpticalProperties_AerosolProfile& other)
{
	bool	ok = true;
	ReleaseResources();
	if (other.m_mieaerosol            != NULL) ok = ok && other.m_mieaerosol->           CreateClone( (skOpticalProperties**)&m_mieaerosol );
	if (other.m_moderadiusclimatology != NULL) ok = ok && other.m_moderadiusclimatology->CreateClone( &m_moderadiusclimatology );
	if (other.m_distribution          != NULL) ok = ok && other.m_distribution->         CreateClone( &m_distribution );
	if (other.m_refractiveindex       != NULL) ok = ok && other.m_refractiveindex->      CreateClone( &m_refractiveindex );
	m_isdirty = true;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfile::DeepCopy, There were errors deep copying the other object\n");
	}
	return ok;
}
*/


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::SetOpticalProperties		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::SetOpticalProperties( skOpticalProperties_ParticleBase* mieaerosol)
{
	if (mieaerosol != NULL) mieaerosol->AddRef();
	if (m_mieaerosol != NULL) m_mieaerosol->Release();
	m_mieaerosol = mieaerosol;
	SetDirty();
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::SetAerosolChemical		2008-2-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::SetRefractiveIndex( skRTRefractiveIndex* ri )
{
	if (ri                != NULL) ri->AddRef();
	if (m_refractiveindex != NULL) m_refractiveindex->Release();
	m_refractiveindex = ri;
	SetDirty();
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::SetParticleDistribution		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::SetParticleDistribution( skRTParticleDist*		distribution )
{
	if ( distribution != NULL) distribution->AddRef();
	if ( m_distribution != NULL) m_distribution->Release();

	m_distribution = distribution;
	SetDirty();
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::SetParticleSizeClimatlogy		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::SetParticleSizeClimatology( skClimatology*	moderadiusclimatology )
{
	if (  moderadiusclimatology != NULL) moderadiusclimatology->AddRef();
	if (m_moderadiusclimatology != NULL) m_moderadiusclimatology->Release();
	m_moderadiusclimatology = moderadiusclimatology;
	SetDirty();
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::ASA_To_ExtinctionPerCm		 2015- 9- 25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::ASA_To_ExtinctionPerCm( double wavenumber, double asa_um2percm3, double* ext)
{
	double	absxs, extxs, scattxs;
	double	N = 0.0;
	bool	ok;

	ok = CalculateCrossSections( wavenumber, &absxs, &extxs, &scattxs);		// Calculate the extinction of 1 particle
	if (ok)
	{
		N =  m_distribution->ASA_To_N( asa_um2percm3);							// Get the number of particles from the aerosol surface area
		ok = NXFINITE(N);
	}
	*ext = ok ? (N*extxs) : std::numeric_limits<double>::quiet_NaN();															// Calculate the extinction per cm
	return ok;																// and return the value
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::ExtinctionPerCm_to_ASA		 2015- 9- 25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::ExtinctionPerCm_to_ASA( double wavenumber, double extinction, double* asa)
{
	double	absxs, extxs, scattxs;
	double	N;
	bool	ok;

	ok = CalculateCrossSections( wavenumber, &absxs, &extxs, &scattxs);		// Calculate the extinction of 1 particle
	if (ok)
	{
		N    =  extinction/extxs;													// Get the number of particles from the total extinction divided by the extinction of 1 particle
		*asa = m_distribution->N_To_ASA(N) ;										// Calculate the aerosol surface area um2percm3.
		ok = NXFINITE(*asa);
	}
	if (!ok) *asa = std::numeric_limits<double>::quiet_NaN();
	return ok;																// and return the value
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::SetAtmosphericState		2008-2-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::SetAtmosphericState( skClimatology*)
{
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::SetAtmosphericState		2008-2-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged )
{
	bool				ok;
	bool				ok2;
	CLIMATOLOGY_HANDLE	paramids[3];
	double				params[3] = {0,0,0};
	size_t				numparams;
	size_t				i;

	m_lastpt = pt;
	CheckDirtyAndUpdate();
	ok = (m_moderadiusclimatology == NULL);																		// If we have a climatology defined
	if (!ok)																									// then
	{																											// we can ask the distribution
		ok = m_distribution->GetDistributionParameterSpeciesID( paramids, N_ELEMENTS(paramids), &numparams);	// What species it wants from the climatology (in the order A,B,C for the distribution parameters)
		if (ok)																									// If that worked
		{																										// Then 
			for (i=0; i < numparams; i++)																		// ask the climatology for the species the distribution wants
			{																									// for each parameter
				ok2 = m_moderadiusclimatology->GetParameter( paramids[i], pt, &params[i], false );				// query the climatology
				ok  = ok && ok2;																				// check if its ok
			}																									// and do all the parameters
		}
		if (ok)																									// if we are good so far then
		{																										// Set the distribution
			m_distribution->SetDistributionParameters( params[0], params[1], params[2] );						// with the parameters we just fetched from the database
			m_mieaerosol->Set_ParticleDistribution( m_distribution );
			if ( crosssectionschanged != NULL) *crosssectionschanged = true;
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfile::SetAtmosphericState, Error setting teh atmospheric state. That might be an issue thats tricky to notice. Be careful");
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::InternalClimatology_UpdateCache		2011-8-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::InternalClimatology_UpdateCache	( const GEODETIC_INSTANT& pt )
{
	bool ok;

	m_lastpt = pt;
	ok = CheckDirtyAndUpdate();
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR,"skOpticalProperties_AerosolProfile::InternalClimatology_UpdateCache, You must successfully configure the aerosol optical props, distribution and refractive index before this call");
	}
	else
	{
		ok =      (m_moderadiusclimatology != NULL);
		ok = ok && m_moderadiusclimatology->UpdateCache( pt );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_AerosolProfile::InternalClimatology_UpdateCache, Error updating the internal aerosol size distribution climatology ");
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::CalculateCrossSections		2008-2-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs)
{
	bool	ok;

	if (m_isdirty) InternalClimatology_UpdateCache	(m_lastpt );
	ok = m_mieaerosol->CalculateCrossSections( wavenumber, absxs, extxs, scattxs);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::CalculatePhaseMatrix		2008-2-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	return m_mieaerosol->CalculatePhaseMatrix( wavenumber, cosscatterangle, phasematrix);
}

bool skOpticalProperties_AerosolProfile::LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff)
{
	return m_mieaerosol->LegendreCoefficientsP11(wavenumber, coeff, usermaxcoeff, opticalmaxcoeff);
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::SetModeRadiusAndWidth		2008-2-29*/
/** Sets the particle size climatology to a single log normal profile of
 *	mode radius in microns and mode width
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::SetLogNormalProfileClimatology( const double* altmeters, const double* moderadius_microns, const double* modewidth, size_t numalt)
{
	nx2dArray<double>					profile;
	bool								ok;
	size_t								idx;
	skClimatology_UserDefinedTable*		lognormalprofile;							//	The height profile of log normal mode radius and mode width
	CLIMATOLOGY_HANDLE					species[2] = {SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS, SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH};


	lognormalprofile = new skClimatology_UserDefinedTable;							// Create a new climatology
	ok = profile.SetSize( numalt, 3 );												// Copy the height profiles n to an 2D array
	if (ok)	
	{
		for (idx = 0; idx < numalt; idx++ )
		{
			profile.At(idx,0) = altmeters[idx];
			profile.At(idx,1) = moderadius_microns[idx];
			profile.At(idx,2) = modewidth[idx];
		}
		ok =       lognormalprofile->LoadProfileFrom2dArray( species, 2, profile );		// Load teh climatology with teh 2-D array
		ok = ok && SetParticleSizeClimatology  ( lognormalprofile );				// Give the climatology to AerosolProfile parent class, it will delete it.
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfileH2SO4::SetModeRadiusAndWidth, Error setting Sulphate log normal parameter profile");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::SetGammaProfileClimatology		2013-1-25*/
/** Sets the particle size climatology to a single gamma profile of
 *	effective radius and rate
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::SetGammaProfileClimatology( const double* altmeters, const double* effectiveradius_microns, const double* rate, size_t numalt)
{
	nx2dArray<double>					profile;
	bool								ok;
	size_t								idx;
	skClimatology_UserDefinedTable*		gammaprofile;							//	The height profile of log normal mode radius and mode width
	CLIMATOLOGY_HANDLE					species[2] = {SKCLIMATOLOGY_GAMMA_EFFECTIVERADIUS_MICRONS, SKCLIMATOLOGY_GAMMA_EFFECTIVEVARIANCE_PERMICRON};


	gammaprofile = new skClimatology_UserDefinedTable;							// Create a new climatology
	ok = profile.SetSize( numalt, 3 );												// Copy the height profiles n to an 2D array
	if (ok)	
	{
		for (idx = 0; idx < numalt; idx++ )
		{
			profile.At(idx,0) = altmeters[idx];
			profile.At(idx,1) = effectiveradius_microns[idx];
			profile.At(idx,2) = rate[idx];
		}
		ok =       gammaprofile->LoadProfileFrom2dArray( species, 2, profile );		// Load teh climatology with teh 2-D array
		ok = ok && SetParticleSizeClimatology  ( gammaprofile );					// Give the climatology to AerosolProfile parent class, it will delete it.
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfile::SetGammaProfileClimatology, Error setting the gamma climatology profile");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::GetDistributionParameter		2013-1-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::GetDistributionParameter( const CLIMATOLOGY_HANDLE& species, const GEODETIC_INSTANT& pt, double* value)
{
	bool	ok;

	CheckDirtyAndUpdate();
	ok = (m_moderadiusclimatology != NULL);
	ok = ok && m_moderadiusclimatology->GetParameter( species,pt, value, false);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfile::GetDistributionParameter, Error fetching aerosol particle dsitribution paramater. Perhaps the climatology is not defined");
		*value = 0.0;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfile::SetModeRadiusAndWidthProfileFromFile		2008-8-12*/
/** Sets the particle size climatology to a single log normal profile of
 *	mode radius in microns and mode width read from a text file.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_AerosolProfile::SetLogNormalProfileClimatologyFromFile ( const char* filename )
{
	nx2dArray<double>  aerosolparams;
	nx1dArray<double>	buffer;
	nx1dArray<double>	ah;
	nx1dArray<double>   moderad;
	nx1dArray<double>	modewidth;
	bool				ok;

	ok = aerosolparams.InputColumnMajorText( filename, 3, 0 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_AerosolProfileH2SO4::SetModeRadiusAndWidthProfileFromFile, Error reading aeraolsol log normal parameters from file %s", (const char*)filename);
	}
	else
	{
		ah        = aerosolparams.XSlice(0, &buffer); // * 1000.0;
		moderad   = aerosolparams.XSlice(1, &buffer);
		modewidth = aerosolparams.XSlice(2, &buffer);
		if (nxarray::Max(ah) < 999.0)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfileH2SO4::SetModeRadiusAndWidthProfileFromFile, It looks like altitude in file %s is in kilometers. I am changing it to meters", (const char*)filename);
			ah *= 1000.0;
		}
		ok = SetLogNormalProfileClimatology( ah.UnsafeArrayBasePtr(), moderad.UnsafeArrayBasePtr(), modewidth.UnsafeArrayBasePtr(), ah.size() );

	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *																	2008-2-29*/
/** Default Height profile of Log Normal Parameters for sulphate aerosol.
**/
/*---------------------------------------------------------------------------*/

//static double g_defaultlognormal_h[]          = {	0,       500,     1500,    2500,    3500,    4500,    5500,    6500,    7500,    8500,    9500,   10500,   11500,   12500,   13500,   14500,   15500,   16500,   17500,   18500,   19500,   20500,   21500,   22500,   23500,   24500,   25500,   26500,   27500,   28500,   29500,   30500,   31500,   32500,   33500,   34500,   35500,   36500,   37500,   38500,   39500,    40500};
//static double g_defaultlognormal_moderadius[] = {	0.18000, 0.18000, 0.17675, 0.17350, 0.17025, 0.16700, 0.16375, 0.16050, 0.15725, 0.15400, 0.15075, 0.14750, 0.14425, 0.14100, 0.13775, 0.13450, 0.13125, 0.12800, 0.12475, 0.12150, 0.11825, 0.11500, 0.11175, 0.10850, 0.10525, 0.10200, 0.09875, 0.09550, 0.09225, 0.08900, 0.08575, 0.08250, 0.07925, 0.07600, 0.07275, 0.06950, 0.06625, 0.06300, 0.05975, 0.05650, 0.05325, 0.05325};
//static double g_defaultlognormal_modewidth[]  = {   1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6,     1.6};



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileH2SO4::skOpticalProperties_AerosolProfileH2SO4		2008-2-29*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_AerosolProfileH2SO4::skOpticalProperties_AerosolProfileH2SO4()
{
	skOpticalProperties_MieAerosolCached*	mieaerosol;
	skRTRefractiveIndex_H2SO4*				ri;
	skRTParticleDist_LogNormal*				distribution;
	bool									ok;
	static const double g_sizeparamalts_data[2]       =  { 0.00, 100000.0};
	static const double	g_moderadius_sulphate_data[2] =  { 0.08, 0.08 };
	static const double g_modewidth_sulphate_data [2] =  { 1.6,  1.6};

	mieaerosol       = new skOpticalProperties_MieAerosolCached;
	distribution     = new skRTParticleDist_LogNormal;
	ri               = new skRTRefractiveIndex_H2SO4;

	distribution->SetDistributionParameters( 0.08, 1.6, 0.0);
	ok =       SetOpticalProperties				( mieaerosol );
	ok = ok && SetRefractiveIndex				( ri );
	ok = ok && SetParticleDistribution			( distribution );
	ok = ok && SetLogNormalProfileClimatology  ( g_sizeparamalts_data, g_moderadius_sulphate_data, g_modewidth_sulphate_data, N_ELEMENTS(g_sizeparamalts_data) );
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileH2SO4::~skOpticalProperties_AerosolProfileH2SO4		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_AerosolProfileH2SO4::~skOpticalProperties_AerosolProfileH2SO4()
{
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileH2SO4::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_AerosolProfileH2SO4::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_AerosolProfileH2SO4*	clone;
	bool								ok;

	clone = new skOpticalProperties_AerosolProfileH2SO4;
	ok    = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfileH2SO4::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->DeepCopy(*this);
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "skOpticalProperties_AerosolProfileH2SO4::CreateClone, error copying this object to clone");
		}
	}
	*userclone = clone;
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileDust::skOpticalProperties_AerosolProfileDust		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_AerosolProfileDust::skOpticalProperties_AerosolProfileDust()
{
	skOpticalProperties_MieAerosolCached*		dustaerosol;
	skRTRefractiveIndex_Dust*					ri;
	skRTParticleDist_LogNormal*					distribution;
	bool										ok;
	static const double g_sizeparamalts_data[2]       =  { 0.00, 100000.0};
	static const double	g_moderadius_sulphate_data[2] =  { 0.08, 0.08 };
	static const double g_modewidth_sulphate_data [2] =  { 1.6,  1.6};


	dustaerosol      = new skOpticalProperties_MieAerosolCached;
	distribution     = new skRTParticleDist_LogNormal;
	ri               = new skRTRefractiveIndex_Dust;

	NXTRACE_ONCEONLY(firsta,("skOpticalProperties_AerosolProfileDust::skOpticalProperties_AerosolProfileDust. Need to put in sensible dust distribution parameters in the constructor\n"));
	distribution->SetDistributionParameters( 0.08, 1.6, 0.0);
	ok =       SetOpticalProperties				( dustaerosol );
	ok = ok && SetRefractiveIndex				( ri );
	ok = ok && SetParticleDistribution			( distribution );
	ok = ok && SetLogNormalProfileClimatology	( g_sizeparamalts_data, g_moderadius_sulphate_data, g_modewidth_sulphate_data, N_ELEMENTS(g_sizeparamalts_data) );
}

/*-----------------------------------------------------------------------------
*					skOpticalProperties_AerosolProfileWater::skOpticalProperties_AerosolProfileWater		2017-03-01*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_AerosolProfileWater::skOpticalProperties_AerosolProfileWater()
{
	skOpticalProperties_MieAerosolCached*		wateraerosol;
	skRTRefractiveIndex_Water*					ri;
	skRTParticleDist_LogNormal*					distribution;
	bool										ok;
	static const double g_sizeparamalts_data[2] = { 0.00, 100000.0 };
	static const double	g_moderadius_sulphate_data[2] = { 0.08, 0.08 };
	static const double g_modewidth_sulphate_data[2] = { 1.6,  1.6 };


	wateraerosol = new skOpticalProperties_MieAerosolCached;
	distribution = new skRTParticleDist_LogNormal;
	ri = new skRTRefractiveIndex_Water;

	NXTRACE_ONCEONLY(firsta, ("skOpticalProperties_AerosolProfileWater::skOpticalProperties_AerosolProfileWater. Need to put in sensible dust distribution parameters in the constructor\n"));
	distribution->SetDistributionParameters(0.08, 1.6, 0.0);
	ok = SetOpticalProperties(wateraerosol);
	ok = ok && SetRefractiveIndex(ri);
	ok = ok && SetParticleDistribution(distribution);
	ok = ok && SetLogNormalProfileClimatology(g_sizeparamalts_data, g_moderadius_sulphate_data, g_modewidth_sulphate_data, N_ELEMENTS(g_sizeparamalts_data));
}

/*-----------------------------------------------------------------------------
*					skOpticalProperties_AerosolProfileWater::skOpticalProperties_AerosolProfileWater		2017-03-01*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_AerosolProfileIce_Mie::skOpticalProperties_AerosolProfileIce_Mie()
{
	skOpticalProperties_MieAerosolCached*		iceaerosol;
	skRTRefractiveIndex_ICE*					ri;
	skRTParticleDist_LogNormal*					distribution;
	bool										ok;
	static const double g_sizeparamalts_data[2] = { 0.00, 100000.0 };
	static const double	g_moderadius_sulphate_data[2] = { 0.08, 0.08 };
	static const double g_modewidth_sulphate_data[2] = { 1.6,  1.6 };


	iceaerosol = new skOpticalProperties_MieAerosolCached;
	distribution = new skRTParticleDist_LogNormal;
	ri = new skRTRefractiveIndex_ICE;

	NXTRACE_ONCEONLY(firsta, ("skOpticalProperties_AerosolProfileIce_Mie::skOpticalProperties_AerosolProfileIce_Mie. Need to put in sensible dust distribution parameters in the constructor\n"));
	distribution->SetDistributionParameters(0.08, 1.6, 0.0);
	ok = SetOpticalProperties(iceaerosol);
	ok = ok && SetRefractiveIndex(ri);
	ok = ok && SetParticleDistribution(distribution);
	ok = ok && SetLogNormalProfileClimatology(g_sizeparamalts_data, g_moderadius_sulphate_data, g_modewidth_sulphate_data, N_ELEMENTS(g_sizeparamalts_data));
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileDust::CreateClone		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_AerosolProfileDust::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_AerosolProfileDust*	clone;
	bool									ok;

	clone = new skOpticalProperties_AerosolProfileDust;
	ok    = (clone != NULL ) && clone->DeepCopy( *this );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfileDust::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->DeepCopy(*this);
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "skOpticalProperties_AerosolProfileDust::CreateClone, error copying this object to clone");
		}
	}
	*userclone = clone;
	return ok;
}
*/


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileIce::skOpticalProperties_AerosolProfileIce		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_AerosolProfileIce::skOpticalProperties_AerosolProfileIce()
{
	skOpticalProperties_IceCrystalCached*	iceaerosol;
	skRTRefractiveIndex_ICE*				ri;
	skRTParticleDist_2Gamma*				distribution;
	bool									ok;
	static const double g_sizeparamalts_data     [2]     =  { 0.00, 100000.0};
	static const double	g_effectiveradius_microns[2]     =  { 0.5, 0.5 };
	static const double g_rate                   [2]     =  { 0.113,  0.113};

	NXTRACE_ONCEONLY(firsta,("skOpticalProperties_AerosolProfileIce::skOpticalProperties_AerosolProfileIce. Need to put in sensible cirrus cloud ice distribution parameters in the constructor\n"));
	iceaerosol       = new skOpticalProperties_IceCrystalCached;
	distribution     = new skRTParticleDist_2Gamma;
	ri               = new skRTRefractiveIndex_ICE;

	distribution->SetDistributionParameters(0.5, 0.113, 0.0 );
	ok =       SetOpticalProperties				( iceaerosol );
	ok = ok && SetRefractiveIndex				( ri );
	ok = ok && SetParticleDistribution			( distribution );
	ok = ok && SetGammaProfileClimatology		( g_sizeparamalts_data, g_effectiveradius_microns, g_rate, N_ELEMENTS(g_sizeparamalts_data));
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileIce::CreateClone		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_AerosolProfileIce::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_AerosolProfileIce*	clone;
	bool								ok;

	clone = new skOpticalProperties_AerosolProfileIce;
	ok    = (clone != NULL ) && clone->DeepCopy( *this );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_AerosolProfileIce::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->DeepCopy(*this);
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "skOpticalProperties_AerosolProfileIce::CreateClone, error copying this object to clone");
		}
	}
	*userclone = clone;
	return ok;
}
*/
