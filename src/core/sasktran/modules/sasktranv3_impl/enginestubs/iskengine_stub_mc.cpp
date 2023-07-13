#include "../dllimplementation/stdafx.h"
#include "skengine_stubs.h"
//#include <regex>

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::ISKEngine_Stub_MC		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_MC::ISKEngine_Stub_MC()
{
	m_radiance.SetReuseMemory( true );
	m_numordersofscatter    = 50;
	m_updateclimatology     = true;
	m_geometryisconfigured  = false;
	m_storeStokes           = false;
	m_storeAMF = false;
	m_storeSecondary = false;
	m_numthreads = 0;
	MakeScalarSetFunctions  ( );
	MakeVectorSetFunctions  ( );
	MakeStringSetFunctions  ( );
	MakeVectorGetFunctions  ( );
	
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::~ISKEngine_Stub_MC		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_MC::~ISKEngine_Stub_MC()
{
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::AddLineOfSight		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::AddLineOfSight( double mjd,  const nxVector& obs, const nxVector& look, int* losindex  )
{
	bool		ok;

	ok = m_linesofsight.AddLineOfSight( obs, look, mjd );
	if (ok)
	{
		*losindex = (int)m_linesofsight.NumRays() - 1;
	}
	else
	{
		*losindex = -999999;
	}
	m_geometryisconfigured = false;
	m_radiance.erase();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::AddSpecies		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::AddSpecies( const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty)
{
	skClimatology*				climptr;
	nxUnknown*					optbaseptr;
	skOpticalProperties*		optptr;
	bool						ok;

	optbaseptr = (opticalproperty != nullptr)? opticalproperty->RawObjectPointer() : nullptr;
	climptr   = dynamic_cast<skClimatology*>( climatology->RawObjectPointer() );
	optptr    = (optbaseptr != nullptr) ? dynamic_cast<skOpticalProperties*>( optbaseptr) : nullptr;
	ok        = m_opticalstate.AddSpecies( species, climptr, optptr);

	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::AddEmission		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::AddEmission( const EMISSION_HANDLE& species, ISKEmission_Stub* emissionobject)
{
	skEmission*				emission;
	bool						ok;

	emission = dynamic_cast<skEmission*>( emissionobject->RawObjectPointer() );
	ok      = m_opticalstate.AddEmission( species, emission );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::SetAlbedo		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::SetAlbedo( double albedo )
{
	return m_opticalstate.SetAlbedo( albedo );
}

/*-----------------------------------------------------------------------------
*					ISKEngine_Stub_MC::SetBRDF		2017-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::SetBRDF(ISKBrdf_Stub* brdf)
{
	skBRDF* brdf_cast = dynamic_cast<skBRDF*>(brdf->RawObjectPointer());

	if (!brdf_cast)
	{
		nxLog::Record(NXLOG_WARNING, "Error in SetBRDF, input BRDF is not a valid BRDF object");
		return false;
	}
	else
	{
		return m_opticalstate.SetAlbedoObject(brdf_cast);
	}
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::SetWavelengths		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::SetWavelengths( const double* wavelen, int numwavelen )
{
	m_wavelen.assign( wavelen, wavelen + numwavelen );
	m_radiance.erase();
	m_specs.SetRadianceWavelengths(m_wavelen);
	return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::InitializeModel		2014-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::InitializeModel()
{
	bool ok = true;
	
	if( !m_geometryisconfigured )
	{
		if(!m_specs.HasBeenFinalized()) ok = ok && m_specs.FinalizeSpecs();
		ok = ok && m_engine.ConfigureModel(m_specs, m_linesofsight, m_numthreads);
		m_geometryisconfigured = ok;
	}
	if( !m_geometryisconfigured ) nxLog::Record(NXLOG_ERROR, "ISKEngine_Stub_MC::CalculateRadiance, Could not configure geometry.");

	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::CalculateRadiance		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::CalculateRadiance(const double** radiance, int* numwavelens, int* numlinesofsight )
{
	std::vector<double>		r;
	bool					ok = true;
	bool					ok1;
	size_t					numrays = m_linesofsight.NumRays();
	size_t					numwave = m_wavelen.size();
	size_t					numcells = m_specs.GetNumAMFCells();
	std::vector<double>		amf, amfvar;

	std::vector< skRTStokesVector >* vecptr = nullptr;
    ok = ok && InitializeModel();
	ok = ok && m_engine.ConfigureMinorResolutionUpdates( m_specs );

	if( m_storeStokes ){
		vecptr = new std::vector< skRTStokesVector >;
		vecptr->resize( numrays );
		m_radiancePol.SetSize( numrays, numwave );
	}

	r.resize( numrays );
	m_radiance.SetSize( numrays, numwave );
	m_variance.SetSize( numrays, numwave );
	if (m_storeAMF) {
		m_airMassFactor.SetSize(numcells, numrays, numwave);
		m_airMassFactorVariance.SetSize(numcells, numrays, numwave);
	}
	if (m_storeSecondary)
	{
		m_secondary.SetSize(numrays, numwave);
		m_secondaryVariance.SetSize(numrays, numwave);
	}

	for (size_t wavidx = 0; wavidx <numwave; ++wavidx)
	{
		ok1 = m_engine.CalculateRadiance( &r, m_wavelen.at(wavidx), m_numordersofscatter, &m_opticalstate, vecptr, m_updateclimatology );
		for (size_t losidx = 0; losidx < numrays; ++losidx)
		{
			m_radiance.At(losidx,wavidx ) = r.at(losidx);
			m_engine.GetMeasurementVariance(losidx, m_variance.At(losidx,wavidx));

			if( m_storeStokes ) m_radiancePol.At(losidx,wavidx) = vecptr->at(losidx);

			if (m_storeAMF) {
				m_engine.GetAirMassFactors(losidx, amf);
				m_engine.GetAirMassFactorVariance(losidx, amfvar);
				for (size_t amfidx = 0; amfidx < numcells; ++amfidx) {
					m_airMassFactor.At(amfidx, losidx, wavidx) = amf[amfidx];
					m_airMassFactorVariance.At(amfidx, losidx, wavidx) = amfvar[amfidx];
				}
			}

			if (m_storeSecondary)
			{
				m_engine.GetSecondaryMeasurement(losidx, m_secondary.At(losidx, wavidx));
				m_engine.GetSecondaryVariance(losidx, m_secondaryVariance.At(losidx, wavidx));
			}
		}
		m_updateclimatology = false;
		ok = ok && ok1;
	}
	*radiance = m_radiance.ArrayBasePtr();
	*numwavelens = (int)numwave;
	*numlinesofsight = (int)numrays;

	delete vecptr;

	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::CalculateRadiancePolarized		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::CalculateStokesVector(const ISKStokesVector** radiance, int* numwavelens , int* numlinesofsight)
{
	std::vector<double>					r;
	bool								ok = true;
	bool								ok1;
	size_t								numrays = m_linesofsight.NumRays();
	size_t								numwave = m_wavelen.size();
	IQUV								iquv;
	GEOGRAPHIC_BASIS					basis_helio, basis_geo;
	ISKBasisDirection					basdir;
	std::vector< skRTStokesVector >		vecptr;
	const SKTRAN_CoordinateTransform_V2*	coords = m_engine.Coordinates();


	ok = ok && m_engine.ConfigureMinorResolutionUpdates( m_specs );

	vecptr.resize( numrays );
	m_radiancePol.SetSize( numrays, numwave );
	m_radiancepolarized.SetSize( numrays, numwave);
	
	r.resize( numrays );
	m_radiance.SetSize( numrays, numwave );
	m_variance.SetSize( numrays, numwave );
	for (size_t wavidx = 0; wavidx <numwave; ++wavidx)
	{
		ok1 = m_engine.CalculateRadiance( &r, m_wavelen.at(wavidx), m_numordersofscatter, &m_opticalstate, &vecptr, m_updateclimatology );
		for (size_t losidx = 0; losidx < numrays; ++losidx)
		{
			m_radiance.At(losidx,wavidx ) = r.at(losidx);
			m_engine.GetMeasurementVariance(losidx, m_variance.At(losidx,wavidx));
			m_radiancePol.At(losidx,wavidx) = vecptr.at(losidx);
			iquv.I = vecptr.at(losidx).At(1);
			iquv.Q = vecptr.at(losidx).At(2);
			iquv.U = vecptr.at(losidx).At(3);
			iquv.V = vecptr.at(losidx).At(4);
			ok = ok && GetBasisHelio( &basis_helio, losidx);
			ok = ok && GetBasisGeo(&basis_geo, losidx);
			basdir.Assign( basis_helio.X(), basis_helio.Y(), basis_helio.Z() );
			m_radiancepolarized.At(losidx,wavidx).Assign( iquv, basdir);
			basdir.Assign( basis_geo.X(), basis_geo.Y(), basis_geo.Z() );
			m_radiancepolarized.At(losidx,wavidx).to_new_basis(basdir);
		}
		m_updateclimatology = false;
		ok = ok && ok1;
	}
	*radiance = m_radiancepolarized.ArrayBasePtr();
	*numwavelens = (int)numwave;
	*numlinesofsight = (int)numrays;


	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::Radiance		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/
//
//double ISKEngine_Stub_MC::Radiance( int userlosindex, int userwavelenindex)
//{
//	bool	ok;
//	size_t	losindex = userlosindex;
//	size_t	windex   = userwavelenindex;
//	double  value;
//
//	ok = (windex < m_radiance.YSize()) && ( losindex < m_radiance.XSize() );
//	value = ok ? m_radiance.At(losindex,windex) : std::numeric_limits<double>::quiet_NaN();
//	return value;
//};

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::SetAtmosphericState		2014-02-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::SetAtmosphericState( ISKClimatology_Stub* climatology )
{
	bool ok = true;

	ok = ok && m_opticalstate.SetAtmosphericStateModel( dynamic_cast<skClimatology*>(climatology->RawObjectPointer()) );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::MakeDefaultOpticalState		2014-3-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::MakeDefaultOpticalState()
{
	skOpticalProperties_RayleighDryAir*		rayleigh;					// Optical properties of one air molecule
	skClimatology_MSIS90*					msis90;	
	skClimatology_LabowOzoneVMR*			o3numberdensity;
	skOpticalProperties_O3_OSIRISRes*		o3_opticalprops;			// optical properties of one O3 molecule
	bool									ok;

	rayleigh         = new skOpticalProperties_RayleighDryAir;
	msis90           = new skClimatology_MSIS90;
	o3numberdensity  = new skClimatology_LabowOzoneVMR;
	o3_opticalprops  = new skOpticalProperties_O3_OSIRISRes;

	m_opticalstate.erase();
	ok =       (rayleigh != NULL) && (msis90 != NULL) && (o3_opticalprops != NULL) && (o3numberdensity != NULL);

	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis90,          rayleigh );
	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_PRESSURE_PA, msis90, rayleigh );
	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_TEMPERATURE_K, msis90, rayleigh );
	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_O3_CM3, o3numberdensity, o3_opticalprops);
	ok = ok && m_opticalstate.SetAlbedo( 1.0 );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_MC::MakeDefaultOpticalState, there was an error making the species list");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::MakeScalarSetFunctions( )
{
	bool ok = true;

	AddSetScalarFunction( "numordersofscatter",
		[&, this](double d)
		{ 
			bool ok = true;
			int scatorder = (int)ceil(d-0.5);
			ok = ok && 0 < scatorder;
			if(ok)
				m_numordersofscatter = scatorder;
			else
				nxLog::Record(NXLOG_WARNING, "ISKEngine_MC, Needs positive order of scatter");
		
			return ok;
		}
	);

	AddSetScalarFunction("numthreads",
		[&, this](double d)
		{ 
			bool ok = true;
			size_t numthreads = (size_t)ceil(d-0.5);
			ok = ok && 0 <= numthreads;
			if(ok)
				m_numthreads = numthreads;
			else
				nxLog::Record(NXLOG_WARNING, "ISKEngine_MC, Needs non-negative number of threads");
		
			return ok;
		}
	);

	AddSetScalarFunction( "setnumthreads",
		[&, this](double d)
		{ 
			bool ok = true;
			size_t numthreads = (size_t)ceil(d-0.5);
			ok = ok && 0 < numthreads;
			if(ok)
				m_numthreads = numthreads;
			else
				nxLog::Record(NXLOG_WARNING, "ISKEngine_MC, Needs non-negative number of threads");
		
			return ok;
		}
	);


	AddSetScalarFunction( "setprecision", 
		[&, this](double d)
		{
			bool ok = true;
			//ok = ok && m_specs.SetNumPhotonsPerLOS      ( (size_t)(1.0/(d*d)) );
			nxLog::Record(NXLOG_WARNING, "ISKEngine_MC, ``SetPrecision'' is deprecated, use ``SetTargetStd''");
			return ok;
		}
	);

	AddSetScalarFunction( "setnumoptpropalts",
		[&, this](double d)
		{
			m_specs.SetNumOptPropAlts( (size_t) d );
			return true;
		}
	);

	AddSetScalarFunction( "settargetstd", 
		[&, this](double d)
		{
			bool ok = true;
			ok = ok && m_specs.SetPrecisionMC( d );
			return ok;
		}
	);

	AddSetScalarFunction( "setminfractionhigherorder", 
		[&, this](double d)
		{
			bool ok = true;
			ok = ok && m_specs.SetMinFractionHigherOrder( d );
			return ok;
		}
	);

	AddSetScalarFunction( "setnumphotonsperlos",
		[&, this](double d)
		{
			bool ok = true;
			ok = ok && m_specs.SetNumPhotonsPerLOS( (size_t)d );
			return ok;
		}
	);

	AddSetScalarFunction( "setnumraytracingshells",
		[&, this](double d)
		{
			bool ok = true;
			ok = ok && m_specs.SetNumRayTracingAlts( (size_t)d );
			return ok;
		}
	);
	
	AddSetScalarFunction( "setsolartablealtitudedelta",
		[&, this](double d)
		{
			bool ok = true;
			//ok = ok && m_specs.SetNumSolTableAlts( (size_t)d );
			ok = ok && d>0.0;
			if(ok){
				ok = ok && m_specs.SetSolarTableAltDelta ( d );
			} else {
				ok = false;
				nxLog::Record(NXLOG_WARNING, "setsolartablealtitudedelta, Solar table altitude spacing must be positive. ");
			}

			return            ok;
		}
	);

	AddSetScalarFunction( "setminrelpathweight", 
		[&, this](double d)
		{
			bool ok = true;
			ok = ok && m_specs.SetMinimumRelPathWeight( d );
			return ok;
		}
	);

	AddSetScalarFunction( "setscatterpositionres",
		[&, this](double d)
		{
			bool ok = true;
			ok = ok && m_specs.SetScatterPositionRes( d );
			return ok;
		}
	);

	AddSetScalarFunction( "setscattangleresolution", 
		[&, this](double d)
		{
			bool ok = true;
			ok = ok = ok && m_specs.SetScattAngleResolution( d );
			return ok;
		}
	);

	AddSetScalarFunction( "setrngseed",
		[&, this](double d)
		{
			bool ok = true;
			std::uint32_t	n = (std::uint32_t)d;
			ok = ok && m_specs.SetRngSeed( n );
			return ok;
		}
	);

	AddSetScalarFunction( "minextinctionratioofcell",
		[&, this](double d)
		{
			bool ok = true;
			ok = ok && m_specs.SetAdaptOptDepthMinRatio( d );
			return ok;
		}
	);

	AddSetScalarFunction( "maxopticaldepthofcell",
		[&, this](double d)
		{
			return m_specs.SetAdaptOptDepthMax( d );
		}
	);

	AddSetScalarFunction( "setmcdebugmode",
		[&, this](double d)
		{
			bool ok = true;
			if( 0.0 != d){
				ok = ok && m_specs.SetAllowDynamicThreads   ( false );
				//ok = ok && m_specs.SetChunkSize             ( 0 );
				this->m_numthreads = 1;
				ok = ok && m_specs.SetRngSeed               ( (std::uint32_t)std::floor(d+0.6) );
			} else{
				ok = ok && m_specs.SetAllowDynamicThreads   ( true );
				//ok = ok && m_specs.SetChunkSize             ( 1 );
				ok = ok && m_specs.SetRngSeed               ( 0 );
			}
			return ok;
		}
	);

	AddSetScalarFunction( "setoptpropinttype",
		[&, this](double d)
		{
			bool ok = true;
			int	quadraturetype = (int)d;

			switch (quadraturetype)
			{
			case 0  :	ok = ok && m_specs.SetOptPropIntType( SKTRAN_Specifications_MC::OptPropIntType::straight );
						break;
			case 1  :	ok = ok && m_specs.SetOptPropIntType( SKTRAN_Specifications_MC::OptPropIntType::adaptive ); 
						break;
			case 2  :   ok = ok && m_specs.SetOptPropIntType( SKTRAN_Specifications_MC::OptPropIntType::constant );
						break;
			default :	ok = false;
						nxLog::Record( NXLOG_WARNING,"ISKEngineMC, invalid integer value [%d] for property SetOptPropintType", (int) quadraturetype);
						break;
			};
			return ok;
		}
	);

	AddSetScalarFunction( "setphotonlogtype",
		[&, this](double d)
		{
			bool ok = true;
			int photonlogtype = (int) ceil(d-0.5);

			switch (photonlogtype)
			{
			case 0  :	ok = ok && m_specs.SetKernelType( SKTRAN_Specifications_MC::LogType::none );
						break;
			case 1  :	ok = ok && m_specs.SetKernelType( SKTRAN_Specifications_MC::LogType::obsPlane );
						break;
			case 2  :	ok = ok && m_specs.SetKernelType( SKTRAN_Specifications_MC::LogType::stDev );
						break;
			case 3  :	ok = ok && m_specs.SetKernelType( SKTRAN_Specifications_MC::LogType::radAlongLOS );
						break;
			case 4  :	ok = ok && m_specs.SetKernelType( SKTRAN_Specifications_MC::LogType::photAlongLOS );
						break;
			case 5  :	ok = ok && m_specs.SetKernelType( SKTRAN_Specifications_MC::LogType::scatPtOnLOS );
						break;
			default :	ok = false;
						nxLog::Record( NXLOG_WARNING,"ISKEngineMC, invalid integer value [%d] for property SetPhotonLogType", (int) photonlogtype);
						break;
			};
			return ok;
		}
	);

	AddSetScalarFunction( "setsolarraytracertype",
		[&, this](double d)
		{
			bool ok = true;
			int raytracertype = (int)d;

			switch (raytracertype)
			{
			case 0 :	ok = ok && m_specs.SetSolarRayTracerType( SKTRAN_Specifications_MC::RayTracerType::shell );
						break;
			case 1 :	ok = ok && m_specs.SetSolarRayTracerType( SKTRAN_Specifications_MC::RayTracerType::curved );
						break;
			case 2 :	ok = ok && m_specs.SetSolarRayTracerType( SKTRAN_Specifications_MC::RayTracerType::generic );
						break;
			default :	ok = false;
						nxLog::Record( NXLOG_WARNING,"ISKEngineMC, invalid integer value [%d] for property SetSolarRayTracerType", (int) raytracertype);
						break;
			};
			return ok;
		}
	);

	AddSetScalarFunction( "setlosraytracertype",
		[&, this](double d)
		{
			bool ok = true;
			int raytracertype = (int)d;

			switch (raytracertype)
			{
			case 0 :	ok = ok && m_specs.SetLOSRayTracerType( SKTRAN_Specifications_MC::RayTracerType::shell );
						break;
			case 1 :	ok = ok && m_specs.SetLOSRayTracerType( SKTRAN_Specifications_MC::RayTracerType::curved );
						break;
			case 2 :	ok = ok && m_specs.SetLOSRayTracerType( SKTRAN_Specifications_MC::RayTracerType::generic );
						break;
			default :	ok = false;
						nxLog::Record( NXLOG_WARNING,"ISKEngineMC, invalid integer value [%d] for property SetLOSRayTracerType", (int) raytracertype);
						break;
			};
			return ok;
		}
	);
	
	AddSetScalarFunction( "setmsraytracertype", 
		[&, this](double d)
		{
			bool ok = true;
			int raytracertype = (int)d;

			switch (raytracertype)
			{
			case 0 :	ok = ok && m_specs.SetMSRayTracerType( SKTRAN_Specifications_MC::RayTracerType::shell );
						break;
			case 1 :	ok = ok && m_specs.SetMSRayTracerType( SKTRAN_Specifications_MC::RayTracerType::curved );
						break;
			case 2 :	ok = ok && m_specs.SetLOSRayTracerType( SKTRAN_Specifications_MC::RayTracerType::generic );
						break;
			default :	ok = false;
						nxLog::Record( NXLOG_WARNING,"ISKEngineMC, invalid integer value [%d] for property SetMSRayTracerType", (int) raytracertype);
						break;
			};
			return ok;
		}
	);

	AddSetScalarFunction( "setpolarizationtype",
		[&, this](double d)
		{
			bool ok = true;
			int opticaltabletype = (int) ceil(d-0.5);

			ok = SetPolarizationMode( opticaltabletype);
			return ok;
		}
	);

	AddSetScalarFunction( "setsolartabletype",
		[&, this](double d)
		{
			bool ok = true;
			int specifier = (int)ceil(d-0.5);
		
			switch(specifier)
			{
			case 0 :	ok = ok && m_specs.SetSolarTableType( SKTRAN_Specifications_MC::SolarTableType::noTable );
				break;
			case 1 :    ok = ok && m_specs.SetSolarTableType( SKTRAN_Specifications_MC::SolarTableType::dim2 );
				break;
			case 2 :    ok = ok && m_specs.SetSolarTableType( SKTRAN_Specifications_MC::SolarTableType::dim3 );
						break;
			case 3 :	ok = ok && m_specs.SetSolarTableType( SKTRAN_Specifications_MC::SolarTableType::doNothing );
						break;
			default :   ok = false;
						break;
			};
			return ok;
		}
	);

	AddSetScalarFunction("setemissiontabletype",
		[&, this](double d)
		{
			bool ok = true;
			int specifier = (int)ceil(d - 0.5);

			switch (specifier)
			{
			case 0:		ok = ok && m_specs.SetEmissionTableType(SKTRAN_Specifications_MC::EmissionTableType::doNothing);
						break;
			case 1:		ok = ok && m_specs.SetEmissionTableType(SKTRAN_Specifications_MC::EmissionTableType::dim1);
						break;
			default:	ok = false;
						break;
			};
			return ok;
		}
	);

	AddSetScalarFunction( "setopticaltabletype",
		[&, this](double d)
		{
			bool ok = true;
			int specifier =(int)ceil(d-0.5);
			switch(specifier)
			{
			case 0:	ok = ok && m_specs.SetOptTableType( SKTRAN_Specifications_MC::OptTableType::dim1 ); break;
			case 1: ok = ok && m_specs.SetOptTableType( SKTRAN_Specifications_MC::OptTableType::dim3_delaunay ); break;
			case 2: ok = ok && m_specs.SetOptTableType( SKTRAN_Specifications_MC::OptTableType::dim1_constant); break;
			default: ok = false; break;
			};
			return ok;
		}
	);

	AddSetScalarFunction( "secondaryoutput",
		[&, this](double d)
		{
			bool ok = true;
			int specifier = (int)ceil(d - 0.5);
			switch (specifier)
			{
				case 0: 
					ok = ok && m_specs.SetSecondaryOutput( SKTRAN_Specifications_MC::SecondaryOutput::none );
					m_storeAMF = false;
					m_storeSecondary = false;
					break;
				case 1: 
					ok = ok && m_specs.SetSecondaryOutput( SKTRAN_Specifications_MC::SecondaryOutput::lengthAMF );
					m_storeAMF = true;
					m_storeSecondary = false;
					break;
				case 2: 
					ok = ok && m_specs.SetSecondaryOutput( SKTRAN_Specifications_MC::SecondaryOutput::opticalDepthAMF );
					m_storeAMF = true;
					m_storeSecondary = false;
					break;
				case 3:
					ok = ok && m_specs.SetSecondaryOutput( SKTRAN_Specifications_MC::SecondaryOutput::elasticRaman );
					m_storeAMF = false;
					m_storeSecondary = true;
					break;
				case 4: 
					ok = ok && m_specs.SetSecondaryOutput( SKTRAN_Specifications_MC::SecondaryOutput::ringSpectrum );
					m_storeAMF = false;
					m_storeSecondary = true;
					break;
				case 5: 
					ok = ok && m_specs.SetSecondaryOutput( SKTRAN_Specifications_MC::SecondaryOutput::fillingInParameter );
					m_storeAMF = false;
					m_storeSecondary = true;
					break;
				default: ok = false; break;
			};
			return ok;
		}
	);

	AddSetScalarFunction( "surfaceheight",
		[&, this](double d)
		{
			return m_specs.SetGroundShiftAlt(d);
		}
	);

	AddSetScalarFunction( "toaheight", 
		[&, this](double d)
		{
			return m_specs.SetTOAHeight(d);
		}
	);

	AddSetScalarFunction("scattertype",
		[&, this](double d)
		{
			bool ok = true;
			int scattertype = (int)ceil(d - 0.5);

			switch (scattertype)
			{
			case 0:	ok = ok && m_specs.SetScatterType(SKTRAN_Specifications_MC::ScatterType::elastic);
				break;
			case 1:	ok = ok && m_specs.SetScatterType(SKTRAN_Specifications_MC::ScatterType::inelastic);
				break;
			case 2: ok = ok && m_specs.SetScatterType(SKTRAN_Specifications_MC::ScatterType::manualInelastic);
				break;
			case 3: ok = ok && m_specs.SetScatterType(SKTRAN_Specifications_MC::ScatterType::both);
				break;
			case 4: ok = ok && m_specs.SetScatterType(SKTRAN_Specifications_MC::ScatterType::manualBoth);
				break;
			default:	ok = false;
				nxLog::Record(NXLOG_WARNING, "ISKEngineMC, invalid integer value [%d] for property SetScatterType", (int)scattertype);
				break;
			};
			return ok;
		}
	);

	AddSetScalarFunction("wavelengthtype",
		[&, this](double d)
		{
			bool ok = true;
			int wltype = (int)ceil(d - 0.5);

			switch (wltype)
			{
			case 0:	ok = ok && m_specs.SetWavelengthType(SKTRAN_Specifications_MC::WavelengthType::single);
				break;
			case 1:	ok = ok && m_specs.SetWavelengthType(SKTRAN_Specifications_MC::WavelengthType::simultaneous);
				break;
			default:	ok = false;
				nxLog::Record(NXLOG_WARNING, "ISKEngineMC, invalid integer value [%d] for property SetWavelengthType", (int)wltype);
				break;
			};
			return ok;
		}
	);

	AddSetScalarFunction("setprimarywavelength",
		[&, this](double d)
		{
			return m_specs.SetPrimaryWavelength(d);
		}
	);

	AddSetScalarFunction("setnadirreferencepointonground",
		[&, this](double d)
		{
			bool ok = true;
			int specifier = (int)ceil(d - 0.5);

			switch (specifier)
			{
			case 0:		ok = ok && m_specs.SetNadirReferencePointOnGround(false);
				break;
			case 1:		ok = ok && m_specs.SetNadirReferencePointOnGround(true);
				break;
			default:	ok = false;
				break;
			};
			return ok;
		}
	);

	return ok;
}

//bool ISKEngine_Stub_MC::SetPropertyScalar( const char* propertyname, double value )
//{
//	bool ok = true;
//
//	nxString str(propertyname);
//	str.MakeLower();
//	auto funciterator = m_scalarsetfunctions.find( str );
//	if( funciterator == std::end(m_scalarsetfunctions) )
//	{
//		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_MC::SetPropertyScalar, this object does not support the scalar property [%s]\n", (const char*)propertyname);
//		ok = false;
//	} else{
//		ok = ok && funciterator->second(value);
//	}
//
//	return ok;
//}

bool ISKEngine_Stub_MC::MakeVectorSetFunctions( )
{
	bool ok = true;

	AddSetVectorFunction( "setsun",
		[&, this](const double* sun, int n) 
		{
			bool ok = true;
			ok = ok && 3==n;
			ok = ok && m_specs.SetSunGeographicPosition( nxVector(sun[0],sun[1],sun[2]) );
			return ok;
		}
	);

	AddSetVectorFunction( "mcprecision_vector",
		[&, this](const double* precisions, int n)
		{
			bool ok = true;

			vector<double> pvec (precisions, precisions + n);
			vector<size_t> nphotons; nphotons.resize(n);
			std::transform(pvec.begin(), pvec.end(), nphotons.begin(), [](double d) { return (size_t)(1.0/(d*d));});
			m_specs.SetNumPhotonsPerLOS( nphotons );
			return ok;
		}
	);

	AddSetVectorFunction( "manualopticalheights", 
		[&, this](const double* h, int n)
		{
			m_specs.SetManualOpticalHeights( std::vector<double>(h,h+n) );
			return true;
		}
	);
	
	AddSetVectorFunction( "manualamfshells", 
		[&, this](const double* h, int n)
		{
			return m_specs.SetManualAMFHeights(std::vector<double>(h, h + n));
		}
	);

	AddSetVectorFunction( "manualraytracingshells",
		[&, this](const double* h, int n)
		{
			return m_specs.SetManualRayTracingHeights(std::vector<double>(h, h + n));
		}
	);

	AddSetVectorFunction( "manualsolartableheights", 
		[&, this](const double* h, int n)
		{
			return m_specs.SetManualSolTableHeights(std::vector<double>(h, h + n));
		}
	);

	AddSetVectorFunction( "opticalpropertieswavelengths",
		[&, this](const double* wl, int n)
		{
			return m_specs.SetOpticalPropertiesWavelengths(std::vector<double>(wl, wl + n));
		}
	);

	AddSetVectorFunction( "radiancewavelengths",
		[&, this](const double* wl, int n)
		{
			return m_specs.SetRadianceWavelengths(std::vector<double>(wl, wl + n));
		}
	);

	AddSetVectorFunction( "solarspectrum",
		[&, this](const double* irr, int n)
		{
			return m_specs.SetSolarSpectrum(std::vector<double>(irr, irr + n));
		}
	);

	AddSetVectorFunction("maxramanorders",
		[&, this](const double* maxorders, int n)
		{
			std::vector<size_t> v(n);
			for (size_t i = 0; i < n; i++) v[i] = (size_t)*(maxorders + i);
			return m_specs.SetMaxRamanOrders(v);
		}
	);

	AddSetVectorFunction("minfractionhigherramanorders",
		[&, this](const double* minfrac, int n)
		{
			return m_specs.SetMinFractionHigherOrder(std::vector<double>(minfrac, minfrac + n));
		}
	);

	AddSetVectorFunction("setreferencepoint",
		[&, this](const double* value, int n)
		{
			bool ok = true;

			ok = ok && (n == 4);
			if (ok)
			{
				ok = m_specs.SetReferencePoint(value[0], value[1], value[2], value[3]);
			}
			return ok;
		}
	);
	return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool ISKEngine_Stub_MC::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
//{
//	bool ok = true;
//
//	nxString str(propertyname);
//	str.MakeLower();
//	auto funciterator = m_vectorsetfunctions.find( str );
//	if( funciterator == std::end(m_vectorsetfunctions) )
//	{
//		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_MC::SetPropertyArray, this object does not support the vector property [%s]\n", (const char*)propertyname);
//		ok = false;
//	} else{
//		ok = ok && funciterator->second( value, numpoints );
//	}
//
//	return ok;
//}



bool ISKEngine_Stub_MC::MakeStringSetFunctions()
{

	AddSetStringFunction( "amfspecies", 
		[&, this](const char* cstr)
		{
			nxString s(cstr);
			if (!s)
			{
				m_specs.SetAMFSpeciesHandle(SKCLIMATOLOGY_UNDEFINED);
			}
			else
			{
				auto* clim = FindGlobalClimatologyHandle((const char*)s);
				if (*clim == SKCLIMATOLOGY_UNDEFINED) {
					nxLog::Record(NXLOG_WARNING, "ISKEngine_MC, Unknown skclimatology handle %s, configuring for geometric AMF", (const char*)s);
					m_specs.SetAMFSpeciesHandle(SKCLIMATOLOGY_UNDEFINED);
				}
				else
				{
					m_specs.SetAMFSpeciesHandle(*clim);
				}
			}
			return true;
		}
	);
	return true;
}

//
//bool ISKEngine_Stub_MC::SetPropertyString(const char* propertyname, const char* value)
//{
//	bool ok = true;
//
//	nxString str(propertyname);
//	str.MakeLower();
//	auto funciterator = m_stringsetfunctions.find(str);
//	if (funciterator == std::end(m_stringsetfunctions))
//	{
//		nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_MC::SetPropertyString, this object does not support the string property [%s]\n", (const char*)propertyname);
//		ok = false;
//	}
//	else {
//		ok = ok && funciterator->second(value);
//	}
//
//	return ok;
//}
//

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool ISKEngine_Stub_MC::SetPropertyObject( const char* propertyname, nxUnknown* object )
//{
//	nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_MC::Set, this object does not support any object properties including [%s]\n", (const char*)propertyname);
//	return false;
//}
//

bool ISKEngine_Stub_MC::GetWeightingFunctions(const double** wf, int* numwavel, int* numlinesofsight, int* numwf)
{	
	bool ok = false;

	nxLog::Record(NXLOG_WARNING, "MC Engine does not support calculation of weighting functions");

	*wf = nullptr;
	*numwavel = 0;
	*numlinesofsight = 0;
	*numwf = 0;

	return ok;
}



//bool ISKEngine_Stub_MC::ParseCommandAndIndex( const nxString& input, nxString& cmd, int& index )
//{
//	nxStringArray	tokens;
//	int				numtoken;
//	bool			ok;
//
//	numtoken = nxStrtok( input, &tokens, "([]) ,:;");
//	ok = (numtoken == 2);
//	if (ok)
//	{
//		cmd = tokens.GetAt(0);
//		index = atoi( tokens.GetAt(1) );
//		cmd.MakeLower();
//	}
//

	//// first check if we have an input with an index
	//std::regex e("^([^\\[\\(]*)[\\[\\(](\\d*)[\\]\\)]");			// i feel bad
	//std::smatch m;
	//std::string str = (const char*)input;

	//bool ok = true;
	//ok = ok && std::regex_search( str, m, e );
	//if( ok )
	//{
	//	// we have an input with an index
	//	ok = (m.size() == 3);
	//	if ( ok )
	//	{
	//		cmd = m[1].str().c_str();
	//		index = atoi(m[2].str().c_str());
	//		cmd.MakeLower();
	//	}
	//}
//	if( !ok )
//	{
//		cmd = input;
//		cmd.MakeLower();
//		index = -1;
//	}
//	return true;
//}


bool ISKEngine_Stub_MC::MakeVectorGetFunctions()
{
	AddGetVectorFunction( "referencepoint",
		[&, this]( int index )
		{
			GEODETIC_INSTANT pt = m_engine.ReferencePoint();
			m_getpropertybuffer.resize(4);
			m_getpropertybuffer[0] = pt.latitude;
			m_getpropertybuffer[1] = pt.longitude;
			m_getpropertybuffer[2] = pt.heightm;
			m_getpropertybuffer[3] = pt.mjd;
			return true;
		}
	);
	//m_vectorgetfunctions[nxString("sun")] = [&, this]( int index )
	//{
	//	nxVector sun = m_engine.GetSun();
	//	m_getpropertybuffer.resize(3);
	//	m_getpropertybuffer[0] = sun.X();
	//	m_getpropertybuffer[1] = sun.Y();
	//	m_getpropertybuffer[2] = sun.Z();
	//	return true;
	//};
	//m_vectorgetfunctions[nxString("observer")] = [&, this]( int index )
	//{
	//	bool ok = true;

	//	ok = ok && m_linesofsight.NumRays() > index;
	//	ok = ok && index >= 0;
	//	if( !ok )
	//	{
	//		nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_HR::GetPropertyArray, The observer array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)index, (int)m_linesofsight.NumRays());
	//	}
	//	/* else */
	//	const SKTRAN_LineOfSightEntry_V2* entry;
	//	ok = ok && m_linesofsight.GetRay( index, &entry );
	//	nxVector pt = entry->Observer();
	//	m_getpropertybuffer.resize(3);
	//	m_getpropertybuffer[0] = pt.X();
	//	m_getpropertybuffer[1] = pt.Y();
	//	m_getpropertybuffer[2] = pt.Z();
	//	return ok;
	//};

	//m_vectorgetfunctions[nxString("look")] = [&, this]( int index )
	//{
	//	bool ok = true;

	//	ok = ok && m_linesofsight.NumRays() > index;
	//	ok = ok && index >= 0;
	//	if( !ok )
	//	{
	//		nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_HR::GetPropertyArray, The observer array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)index, (int)m_linesofsight.NumRays());
	//	}
	//	/* else */
	//	const SKTRAN_LineOfSightEntry_V2* entry;
	//	ok = ok && m_linesofsight.GetRay( index, &entry );
	//	nxVector pt = entry->Look();
	//	m_getpropertybuffer.resize(3);
	//	m_getpropertybuffer[0] = pt.X();
	//	m_getpropertybuffer[1] = pt.Y();
	//	m_getpropertybuffer[2] = pt.Z();
	//	return ok;
	//};

	AddGetVectorFunction( "stokesvec",
		[&, this]( int index )
		{
			bool ok = true;

			ok = ok && 0<=index && index< (int)(m_radiancePol.XSize()*m_radiancePol.YSize());
			if(ok){
				int numrays = (int) m_radiance.XSize();
				int numwav  = (int) m_radiance.YSize();
				int xindex = index % numrays;
				int yindex = index / numrays;
				m_getpropertybuffer.resize(4);
			
				if( m_storeStokes ){
					m_getpropertybuffer[0] = m_radiancePol.At(xindex, yindex).At(1);
					m_getpropertybuffer[1] = m_radiancePol.At(xindex, yindex).At(2);
					m_getpropertybuffer[2] = m_radiancePol.At(xindex, yindex).At(3);
					m_getpropertybuffer[3] = m_radiancePol.At(xindex, yindex).At(4);
				} else{
					m_getpropertybuffer[0] = m_radiance.At(xindex, yindex);
					m_getpropertybuffer[1] = 0.0;
					m_getpropertybuffer[2] = 0.0;
					m_getpropertybuffer[3] = 0.0;
				}
			} else{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_MC::GetPropertyArray, No stokes vector stored for index %i.", (int)index );
			}
			return ok;
		}
	);

	AddGetVectorFunction( "basis",
		[&, this]( int index )
		{
			bool ok = true;

			ok = ok && 0 <= index && index < (int)m_linesofsight.NumRays();
			if( ok )
			{
				GEOGRAPHIC_BASIS basis;
				ok = ok && GetBasisGeo( &basis, index);
				m_getpropertybuffer.resize(9);
				m_getpropertybuffer[0] = basis.X().X();
				m_getpropertybuffer[1] = basis.X().Y();
				m_getpropertybuffer[2] = basis.X().Z();
				m_getpropertybuffer[3] = basis.Y().X();
				m_getpropertybuffer[4] = basis.Y().Y();
				m_getpropertybuffer[5] = basis.Y().Z();
				m_getpropertybuffer[6] = basis.Z().X();
				m_getpropertybuffer[7] = basis.Z().Y();
				m_getpropertybuffer[8] = basis.Z().Z();
			} else{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_MC::GetPropertyArray, The observer array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)index, (int)m_linesofsight.NumRays());
			}
			return ok;
		}
	);

	AddGetVectorFunction( "variance",
		[&, this]( int  )
		{
			bool ok = true;

			int numrays = (int) m_variance.XSize();
			int numwav  = (int) m_variance.YSize();

			m_getpropertybuffer.resize( numrays*numwav );
			for( int index=0; index<numrays*numwav; ++index){
				int yindex = index / numrays;
				int xindex = index % numrays;
				m_getpropertybuffer[index] = m_variance.At(xindex, yindex);  
			}
			return ok;
		}
	);

	AddGetVectorFunction( "airmassfactor", 
		[&, this](int)
		{
			bool ok = true;
			ok = ok && m_storeAMF;

			if (ok)
			{
				int numlayers = (int)m_airMassFactor.XSize();
				int numrays = (int)m_airMassFactor.YSize();
				int numwav = (int)m_airMassFactor.ZSize();

				m_getpropertybuffer.resize(numrays*numwav*numlayers);
				for (int index = 0; index < numrays*numwav*numlayers; ++index) {
					int zindex = index / (numlayers * numrays);
					int yindex = (index % (numlayers * numrays)) / numlayers;
					int xindex = (index % (numlayers * numrays)) % numlayers;

					m_getpropertybuffer[index] = m_airMassFactor.At(xindex, yindex, zindex);
				}
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_MC::GetPropertyArray, no AMFs were stored.");
			}

			return ok;
		}
	);

	AddGetVectorFunction( "airmassfactorvariance", 
		[&, this](int)
		{
			bool ok = true;
			ok = ok && m_storeAMF;

			if (ok)
			{
				int numlayers = (int)m_airMassFactorVariance.XSize();
				int numrays = (int)m_airMassFactorVariance.YSize();
				int numwav = (int)m_airMassFactorVariance.ZSize();

				m_getpropertybuffer.resize(numrays*numwav*numlayers);
				for (int index = 0; index < numrays*numwav*numlayers; ++index) {
					int zindex = index / (numlayers * numrays);
					int yindex = (index % (numlayers * numrays)) / numlayers;
					int xindex = (index % (numlayers * numrays)) % numlayers;

					m_getpropertybuffer[index] = m_airMassFactorVariance.At(xindex, yindex, zindex);
				}
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_MC::GetPropertyArray, no AMFs were stored.");
			}
			
			return ok;
		}
	);

	AddGetVectorFunction("secondarymeasurement",
		[&, this](int)
		{
			bool ok = true;
			ok = ok && m_storeSecondary;

			if (ok)
			{
				int numrays = (int)m_secondary.XSize();
				int numwav = (int)m_secondary.YSize();

				m_getpropertybuffer.resize(numrays*numwav);
				for (int index = 0; index < numrays*numwav; ++index) {
					int yindex = index / numrays;
					int xindex = index % numrays;;

					m_getpropertybuffer[index] = m_secondary.At(xindex, yindex);
				}
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_MC::GetPropertyArray, no secondary measurements were stored.");
			}

			return ok;
		}
	);

	AddGetVectorFunction("secondaryvariance",
		[&, this](int)
		{
			bool ok = true;
			ok = ok && m_storeSecondary;

			if (ok)
			{
				int numrays = (int)m_secondaryVariance.XSize();
				int numwav = (int)m_secondaryVariance.YSize();

				m_getpropertybuffer.resize(numrays*numwav);
				for (int index = 0; index < numrays*numwav; ++index) {
					int yindex = index / numrays;
					int xindex = index % numrays;;

					m_getpropertybuffer[index] = m_secondaryVariance.At(xindex, yindex);
				}
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_MC::GetPropertyArray, no secondary measurements were stored.");
			}

			return ok;
		}
	);

	return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::GetBasis		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::GetBasisHelio( GEOGRAPHIC_BASIS* basis, size_t losindex)
{
	bool								ok = true;
	const SKTRAN_LineOfSightEntry_V2*	entry;
			
	ok = ok && nullptr!=basis;
	ok = ok && m_linesofsight.GetRay( losindex, &entry );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_HR::GetBasis, error fetching ray (%d). Cannot fetch BASIS for this ray", (int)losindex);
	}
	else
	{
		const HELIODETIC_POINT observerh;

		HELIODETIC_VECTOR observervector = m_engine.Coordinates()->GeographicToHelio(m_engine.Coordinates()->TranslateGeoidToOsculatingSphere(entry->Observer()));
		HELIODETIC_POINT observerpoint;

		observerpoint.FromVector(observervector, m_engine.Coordinates());

		HELIODETIC_UNITVECTOR lookh;
		lookh.FromVector(m_engine.Coordinates()->GeographicToHelio(entry->Look()));

		HELIODETIC_BASIS heliobasis(observerpoint, lookh);

		const nxVector propagation = m_engine.Coordinates()->HelioUnitVectorToGeographic(heliobasis.X());
		const nxVector solar_theta = m_engine.Coordinates()->HelioUnitVectorToGeographic(heliobasis.Y());
		const nxVector solar_phi = m_engine.Coordinates()->HelioUnitVectorToGeographic(heliobasis.Z());

		const nxVector geo_phi = entry->Observer().UnitVector().Cross(propagation).UnitVector();
		const nxVector geo_theta = entry->Observer();

		*basis = GEOGRAPHIC_BASIS(propagation, solar_theta, solar_phi);
	}
	return ok;
}

bool ISKEngine_Stub_MC::GetBasisGeo( GEOGRAPHIC_BASIS* basis, size_t losindex)
{
	bool								ok = true;
	const SKTRAN_LineOfSightEntry_V2*	entry;
			
	ok = ok && nullptr!=basis;
	ok = ok && m_linesofsight.GetRay( losindex, &entry );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_HR::GetBasis, error fetching ray (%d). Cannot fetch BASIS for this ray", (int)losindex);
	}
	else
	{
		const nxVector obs_translated = m_engine.Coordinates()->TranslateGeoidToOsculatingSphere(entry->Observer());
		const nxVector propagation = -1 * entry->Look();
		const nxVector geo_theta = obs_translated.UnitVector().Cross(propagation).UnitVector();
		const nxVector geo_phi = propagation.Cross(geo_theta);

		*basis = GEOGRAPHIC_BASIS(propagation, geo_theta, geo_phi);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::GetPropertyArray		 2014- 5- 16*/
/** **/
/*---------------------------------------------------------------------------*/

//bool ISKEngine_Stub_MC::GetProperty( const char* propertyname, const double** value, int* numpoints )
//{
//	nxString name( propertyname);
//	nxString cmd;
//	int index;
//	bool ok = ParseCommandAndIndex( name, cmd, index );
//	auto funcinterator = m_vectorgetfunctions.find( cmd );
//	if( funcinterator == std::end( m_vectorgetfunctions ) )
//	{
//		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_MC::GetPropertyArray, The MC Array Property does not support parameter %s", (const char*) cmd );
//		*numpoints = 0;
//		*value     = NULL;
//		return false;
//	}
//	else
//	{
//		ok = funcinterator->second(index);
//		*numpoints = (int)m_getpropertybuffer.size();
//		*value = &m_getpropertybuffer[0];
//		return ok;
//	}
//}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::SetPolarizationMode		 2015- 10- 30*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_MC::SetPolarizationMode( int opticaltabletype)
{
	bool ok = true;

	if (opticaltabletype > 2) opticaltabletype = 3;
	switch (opticaltabletype)
	{
		case 0 :	
			m_storeStokes = false;
			ok = ok && m_specs.SetPolType( SKTRAN_Specifications_MC::PolType::none );
			break;
		case 1 :	
			m_storeStokes = true;
			ok = ok && m_specs.SetPolType( SKTRAN_Specifications_MC::PolType::pv1 );
			break;
		case 2:
			m_storeStokes = true;
			ok = false;
			nxLog::Record( NXLOG_WARNING, "ISKEngine MC, pseudovector2 is not implemented.");
			break;
		case 3 :	
			m_storeStokes = true;
			ok = ok && m_specs.SetPolType( SKTRAN_Specifications_MC::PolType::pol );
			break;
		default :
			ok = false;
			nxLog::Record( NXLOG_WARNING,"ISKEngine MC, invalid integer value [%d] for property SetPolarizationType", (int) opticaltabletype);
			break;
	};
	return ok;
};


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC::GetPropertyObject		 2014- 5- 16*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool ISKEngine_Stub_MC::GetPropertyObject( const char* propertyname, nxUnknown** object )
//{
//	nxString name( propertyname);
//
//	*object = nullptr;
//	nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_MC::GetPropertyObject, The MC Object Property does not support parameter %s", (const char*) name );
//	return false;
//}
//

/*
%%%%%%%%%%%%%%
% MC-specific:
%%%%%%%%%%%%%%

SetNumPhotonsPerLOS( size_t N )
  Set how many ray paths will be simulated. Precision (stdev/result) varies as N^{-1/2}. 
Values:
  N>=0 -- The number of photons to simulate. 
Default:
  N = 1e4, corresponds to precision 1%. 


SetMinRelPathWeight( double w )
  MC simulation of a single photon path stops at some order if each higher-order scatter could add to the measured radiance by less than w*{amountAlreadyMeasured}. This can save a lot of time, as only the odd "important" high-order path has to be simulated, but it does introduce a small bias. The importance of each order of scattering decreases approximately geometrically (except if there are bright clouds, etc), so only the first truncated photon contributes significantly to this bias. 
Values:
  w>=0.0 -- Specifies a criterion for truncating the simulation of a path.
Default:
  w = 0.0 -- No truncation. 
Recommended values:
  Choosing  w = 1 / ( 3 * numPhotonsPerLOS^2 )  makes the negative bias approximately one third of the magnitude of random noise in the experiment. Very often a photon could make at most a contribution much less than w, so the bias is actually much smaller than might be expected. 


SetScatterPositionRes( double d )
  Simulated scatter positions are determined to within d meters of the "physical" scatter positions. 
Values:
  d>=0.0 -- Maximum error between simulated scattering positions and physical solution, in meters. 
Default:
  d = 50.0;
Recommended values: 
  d = 50.0. Decreasing d has little impact on performance. Scatter positions are also guaranteed to be between the correct set of quadrature points (e.g. if ray tracing altitudes spaced 10.0 meters apart the scatter position is determined to at least 10.0 meters, regardless of d). 


SetScattAngleResolution  ( double d )
  All optical scattering properties are precalculated at some distribution uniform in cos(scatteringAngle). The number of points in this distribution is the same as the number of points on a circle separated by d degrees. 
Values:
  d > 0.0
Default:
  d = 0.5 corresponds to 361 points distributed evenly in cos(scatterAngle).


SetRngSeed( std::uint32_t n )
    Control how random number generator values are seeded. 
Values:
  n==0 -- Seed is generated from system clock at runtime. Seeds are guaranteed to be different for each RNG. 
  n >0 -- Seed is n for all RNGs. This implies all threads will perform an identical sequence of trials, and precision worsens by sqrt(n). 
Default:
  n=0;
Recommended values: 
  n>0 should only be used for debugging. To guarantee exactly the same results are produced between runs the user has to use options:
	ok = ok && mc_specs.SetAllowDynamicThreads   ( false );
	ok = ok && mc_specs.SetChunkSize             ( 0 );
	ok = ok && mc_specs.SetRngSeed               ( n );
Maybe there should be just one interface option "SetMCDebugMode" that sets all those options at once so the user doesn't have to see the RNG seed.
The non-debug-mode options are:
	ok = ok && mc_specs.SetAllowDynamicThreads   ( true );
	ok = ok && mc_specs.SetChunkSize             ( 1 );
	ok = ok && mc_specs.SetRngSeed               ( 0 );


%%%%%%%%%%%%%%%%%%%%
%%% Not MC-specific: 
%%%%%%%%%%%%%%%%%%%%

SetOptPropIntType (SKTRAN_Specifications_MC_OptPropIntType t)
  Set the optical property table integrator (quadrature) type.
Values:
  SKTRAN_Specifications_MC::opiType_straight -- Use the v21-style integrator
  SKTRAN_Specifications_MC::opiType_adaptive -- Use Dan's adaptive integrator
Default:
  opiType_straight;
Recommended values:
  d = opiType_straight for thin atmospheres, opiType_adaptive for higher precision in optically thick atmospheres. 


SetSolarRayTracerType    (SKTRAN_Specifications_MC_RayTracerType  t)
SetLOSRayTracerType      (SKTRAN_Specifications_MC_RayTracerType  t)
SetMSRayTracerType       (SKTRAN_Specifications_MC_RayTracerType  t)
  Set the ray tracer type for rays from the sun into the atmosphere, rays for the observer line of sight, and rays traced during multiple scattering. 
Values:
  SKTRAN_Specifications_MC::rtType_shell  -- V21-style straight rays. 
  SKTRAN_Specifications_MC::rtType_curved -- Lorne's curved ray algorithm.
Default:
  rtType_shell for all.
Recommended values:
  Solar rays should be straight -- we don't have a good algorithm for tracing curved rays from the sun to a point. 
  MS rays should probably be straight -- using curved rays has almost no effect and is a big performance hit. Using curved MS rays might be desireable for some people. 
  LOS rays can be either straight or curved -- the performance hit for curved rays is pretty small here.
  */
