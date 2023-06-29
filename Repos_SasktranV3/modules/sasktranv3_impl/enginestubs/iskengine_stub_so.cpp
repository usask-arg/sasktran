#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::ISKEngine_Stub_SO		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_SO::ISKEngine_Stub_SO()
{
	m_radiance.SetReuseMemory(true);
	m_numordersofscatter = 50;
	m_updateclimatology  = true;
	m_modelisinitialized  = false;
	MakeScalarSetFunctions();
	MakeVectorSetFunctions();
	MakeGetFunctions();
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::~ISKEngine_Stub_SO		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_SO::~ISKEngine_Stub_SO()
{
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::AddLineOfSight		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::AddLineOfSight( double mjd,  const nxVector& obs, const nxVector& look, int* losindex  )
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
	m_radiance.erase();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::AddSpecies		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::AddSpecies( const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty)
{
	skClimatology*				climptr;
	nxUnknown*					optbaseptr;
	skOpticalProperties*		optptr;
	bool						ok;

	optbaseptr = (opticalproperty != nullptr)? opticalproperty->RawObjectPointer() : nullptr;
	climptr   = dynamic_cast<skClimatology*>( climatology->RawObjectPointer() );
	optptr    = (optbaseptr != nullptr) ? dynamic_cast<skOpticalProperties*>( optbaseptr) : nullptr;

	ok      = (climptr != nullptr); // && (optptr != nullptr);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_SO::AddSpecies, Cannot convert the climatology and or optical property objects to corresponding C++ pointers");
	}
	else
	{
		ok      = m_opticalstate.AddSpecies( species, climptr, optptr);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::AddSpecies		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::AddEmission( const EMISSION_HANDLE& species, ISKEmission_Stub* emissionobject)
{
	skEmission*				emission;
	bool						ok;

	emission = dynamic_cast<skEmission*>( emissionobject->RawObjectPointer() );
	ok      = m_opticalstate.AddEmission( species, emission );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::SetAlbedo		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::SetAlbedo( double albedo )
{
	return m_opticalstate.SetAlbedo( albedo );
}

/*-----------------------------------------------------------------------------
*					ISKEngine_Stub_SO::SetBRDF		2017-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::SetBRDF(ISKBrdf_Stub* brdf)
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
 *					ISKEngine_Stub_SO::SetPolarizationMode		 2015- 10- 30*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::SetPolarizationMode( int mode )
{
	bool ok;

	ok = (mode == 0);
	if (!ok) nxLog::Record(NXLOG_WARNING," ISKEngine SO, The SO engine only supports the scalar mode for polarization settings");
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::SetWavelengths		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::SetWavelengths( const double* wavelen, int numwavelen )
{
	m_wavelen.assign( wavelen, wavelen + numwavelen );
	m_radiance.erase();
	return true;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::InitializeModel		2014-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::InitializeModel()
{
	m_modelisinitialized = m_engine.ConfigureModel( m_specs, m_linesofsight, 0 );
	return m_modelisinitialized;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::CalculateRadiance		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::CalculateRadiance( const double** radiance, int* numwavelens , int* numlinesofsight)
{
	std::vector<double>		r;
	bool					ok = true;
	bool					ok1;
	size_t					numrays = m_linesofsight.NumRays();
	size_t					numwave = m_wavelen.size();

	if (!m_modelisinitialized) ok = InitializeModel();
	r.resize( numrays );
	ok = ok && m_radiance.SetSize(numrays, numwave);
	ok = ok && (numrays > 0) && (numwave > 0);

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_SO::CalculateRadiance, There is nothing to do as either the number of rays (%d) is zero or the number of wavelengths (%d) is zero", (int)numrays, (int)numwave);
	}
	else
	{
		for (size_t i = 0; i <numwave; i++)
		{
			ok1 = m_engine.CalculateRadiance( &r, m_wavelen.at(i), m_numordersofscatter, &m_opticalstate, nullptr, m_updateclimatology );
			for (size_t l = 0; l < numrays; l++)
			{
				m_radiance.At( l, i ) = ok1 ? r.at(l) : std::numeric_limits<double>::quiet_NaN();
			}
			m_updateclimatology = false;
			ok = ok && ok1;
		}
	}
	*radiance = m_radiance.ArrayBasePtr();
	*numwavelens = (int)numwave;
	*numlinesofsight = (int)numrays;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::CalculateRadiancePolarized		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::CalculateStokesVector( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight)
{
	const double*						r;
	bool								ok;
	size_t								numwave;
	size_t								numrays;
	IQUV								iquv;
	ISKBasisDirection					basdir;
	double								nan = std::numeric_limits<double>::quiet_NaN();
	nxVector							nullvector(nan, nan, nan);
		
	basdir.Assign( nullvector, nullvector, nullvector);
	ok = CalculateRadiance( &r, numwavelens, numlinesofsight );
	if (ok)
	{
		numwave = *numwavelens;
		numrays = *numlinesofsight;
		ok = (numrays > 0) && (numwave > 0);
		ok = ok && m_radiancepolarized.SetSize(numrays,numwave); 
		if (ok)
		{
			for (size_t i = 0; i <numwave; i++)
			{
				for (size_t l = 0; l < numrays; l++)
				{
					iquv.I = m_radiance.At( l, i );
					iquv.Q = 0.0;
					iquv.U = 0.0;
					iquv.V = 0.0;
					m_radiancepolarized.At(l,i).Assign( iquv, basdir);
				}
			}
			*radiancep = m_radiancepolarized.ArrayBasePtr();
		}
	}
	if (!ok) *radiancep = nullptr;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::GetWeightingFunctions		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::GetWeightingFunctions(const double** wf, int* numwavel, int* numlinesofsight, int* numwf)
{
	bool ok = false;

	nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_SO::GetWeightingFunctions, The SO Engine does not support calculation of weighting functions");

	*wf = nullptr;
	*numwavel = 0;
	*numlinesofsight = 0;
	*numwf = 0;

	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::Radiance		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

//double ISKEngine_Stub_SO::Radiance( int userlosindex, int userwavelenindex)
//{
//	bool	ok;
//	size_t	losindex = userlosindex;
//	size_t	windex   = userwavelenindex;
//	double  value;
//
//	ok = (windex < m_radiance.YSize()) && ( losindex < m_radiance.XSize());
//	value = ok ? m_radiance.At(losindex,windex) : std::numeric_limits<double>::quiet_NaN();
//	return value;
//};


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::SetAtmosphericState		2014-02-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_SO::SetAtmosphericState( ISKClimatology_Stub* climatology )
{
	bool	ok;
	ok = m_opticalstate.SetAtmosphericStateModel( dynamic_cast<skClimatology*>(climatology->RawObjectPointer()) );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::CheckModelNotInitialized		 2015- 11- 12*/
/** **/
/*---------------------------------------------------------------------------*/

bool  ISKEngine_Stub_SO::CheckModelNotInitialized()
{
	bool	ok;
	ok = !m_modelisinitialized;
	if (!ok) nxLog::Record(NXLOG_WARNING,"ISKEngine SO, You cannot set this property after the model after calling InitializeModel or CalculateRadiance");
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/


void ISKEngine_Stub_SO::MakeScalarSetFunctions( )
{

	AddSetScalarFunction( "numordersofscatter",
		[&, this](double value)
		{ 
			bool ok;
			m_numordersofscatter = (size_t)(value + 0.5);
			if (m_numordersofscatter < 1)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(SO)::SetPropertyScalar(NumOrdersOfScatter), The number of orders of scatter must be greater than 0. I am setting it to 1");
				m_numordersofscatter = 1;
				ok = false;
			}
			if 	( m_numordersofscatter > 1000)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(SO)::SetPropertyScalar(NumOrdersOfScatter), The number of orders of scatter (%d) should not be greater than 1000, Truncating to 1000", (int)m_numordersofscatter);
				m_numordersofscatter = 1000;
				ok = false;
			}
			return ok;
		}
	);
	
	AddSetScalarFunction( "numdiffuseprofiles",
		[&, this](double value)
		{ 
			bool ok;
			ok = CheckModelNotInitialized() && m_specs.DiffuseSpecificationsVar()->ConfigureNumberDiffuseProfiles( (size_t)(value+0.5) );
			return ok;
		}
	);

	AddSetScalarFunction( "setreferencepoint_targetaltitude",
		[&, this](double value)
		{ 
			bool ok;
			ok = CheckModelNotInitialized() && m_specs.RayTracingRegionManagerVar()->SetReferencePointTargetAltitude( value);
			return ok;
		}
	);

	AddSetScalarFunction( "setreferencepoint_targetrange",
		[&, this](double value)
		{ 
			bool ok;
			ok = CheckModelNotInitialized() && m_specs.RayTracingRegionManagerVar()->SetReferencePointTargetRange(value);
			return ok;
		}
	);

	AddSetScalarFunction( "setgroundaltitude",
		[&, this](double value)
		{ 
			bool ok;
			ok = CheckModelNotInitialized() && m_specs.RayTracingRegionManagerVar()->SetGroundAltitude( value);
			return ok;
		}
	);

	AddSetScalarFunction( "setupperboundaltitude",
		[&, this](double value)
		{ 
			bool ok;
			ok = CheckModelNotInitialized() && m_specs.RayTracingRegionManagerVar()->SetUpperBoundAltitude(value);
			return ok;
		}
	);

	AddSetScalarFunction( "setlowerboundaltitude",
		[&, this](double value)
		{ 
			bool ok;
			ok = CheckModelNotInitialized() && m_specs.RayTracingRegionManagerVar()->SetLowerBoundAltitude( value);
			return ok;
		}
	);

}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::SetPropertyArray		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKEngine_Stub_SO::MakeVectorSetFunctions()
{

	AddSetVectorFunction( "setsun",
		[&, this](const double* value, int numpoints) 
		{
			bool ok;
			ok = (numpoints == 3);
			if (ok)
			{
				nxVector	sun( value[0], value[1], value[2] );
				ok = CheckModelNotInitialized() && m_specs.RayTracingRegionManagerVar()->SetSun( sun );
			}
			else
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(SO)::SetPropertyArray(SetSun), The sun must be specied as a 3 element double array");
			}
			return ok;
		}
	);

	AddSetVectorFunction( "setreferencepoint",
		[&, this](const double* value, int numpoints) 
		{
			bool ok;
			ok = (numpoints == 4);
			if (ok)
			{
				ok = CheckModelNotInitialized() && m_specs.RayTracingRegionManagerVar()->SetReferencePoint( value[0], value[1], value[2], value[3] );
			}
			else
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(SO)::SetPropertyArray(SetReferencePoint), The reference must be specied as a 2 element double array [latitude, longitude]");
			}
			return ok;
		}
	);

	AddSetVectorFunction( "configureuserdefinedshells",
		[&, this](const double* value, int numpoints) 
		{
			bool	ok;
			std::vector<double>	shellalts;

			shellalts.assign( value, value+numpoints);
			ok = CheckModelNotInitialized() && m_specs.ConfigureUserDefinedShells(shellalts);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(SO)::SetPropertyArray(ConfigureUserDefinedShells), Error configuring the user defined shell altitudes");
			}
			return ok;
		}
	);

	AddSetVectorFunction( "manualdiffuseheights",
		[&, this](const double* value, int numpoints) 
		{
			bool ok;
			ok = (numpoints > 0);
			ok = ok && CheckModelNotInitialized() && m_specs.DiffuseSpecificationsVar()->ConfigureDiffuseAltitudeResolution( value, numpoints );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(SO)::SetPropertyArray(ManualDiffuseHeights), Error configuring the user defined shell altitudes");
			}
			return ok;
		}
	);

	AddSetVectorFunction( "diffuseincomingresolution",
		[&, this](const double* value, int numpoints) 
		{
			bool ok;
			bool ok1;
			ok = (numpoints == 4);
			ok = ok && CheckModelNotInitialized();
			if (ok)
			{
				size_t	groundres  = (size_t) ceil(value[0]-0.5);
				size_t	horizonres = (size_t) ceil(value[1]-0.5);
				size_t	atmosres   = (size_t) ceil(value[2]-0.5);
				size_t	numazi     = (size_t) ceil(value[3]-0.5);
				nx1dArray<double>	diffuseazi;

				if (( groundres > 0) && (horizonres > 0) & (atmosres > 0))
				{
					ok =  m_specs.DiffuseSpecificationsVar()->ConfigureIncomingZenithResolutions	( groundres, horizonres, atmosres, false );
				}

				if (numazi > 0)
				{
					diffuseazi.Indgen(numazi)   *=  (360.0/numazi);	// Get the azimuth angles same as Sasktran V1
					ok1 = m_specs.DiffuseSpecificationsVar()->ConfigureIncomingAzimuthResolution( diffuseazi.ArrayBasePtr(), numazi );
					ok = ok && ok1;
				}
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(SO)::SetPropertyArray(DiffuseIncomingResolution), Error configuring the user defined resolutions");
			}
			return ok;
		}
	);
	AddSetVectorFunction( "configureincomingzenithangles",
		[&, this](const double* zenithdegrees, int numpoints) 
		{
			bool ok;
			double * zenptr = (double *)(intptr_t)zenithdegrees;

			ok = (numpoints >  0);
			ok = ok && CheckModelNotInitialized();
			if (ok)
			{
				ok =  m_specs.DiffuseSpecificationsVar()->ConfigureUserDefinedIncomingZenithAngle( zenptr, numpoints);
			}

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(SO)::SetProperty(diffuseincomingzenithangles), Error configuring the user defined incoming zenith angleresolutions");
			}
			return ok;
		}
	);
	AddSetVectorFunction( "configureincomingazimuthangles",
		[&, this](const double* azidegrees, int numpoints) 
		{
			bool ok;
			double * aziptr = (double *)(intptr_t)azidegrees;

			ok = (numpoints >  0);
			ok = ok && CheckModelNotInitialized();
			if (ok)
			{
				ok =  m_specs.DiffuseSpecificationsVar()->ConfigureIncomingAzimuthResolution( aziptr, numpoints);
			}

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(SO)::SetProperty(diffuseincomingzenithangles), Error configuring the user defined incoming zenith angleresolutions");
			}
			return ok;
		}
	);

}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::GetPropertyArray		 2014- 5- 16*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKEngine_Stub_SO::MakeGetFunctions()
{
	AddGetScalarFunction( "numlinesofsight",
		[&, this]( double* val )
		{
			*val = (double)m_linesofsight.NumRays();
			return true;
		}
	);

	AddGetVectorFunction( "referencepoint",
		[&, this]( int index )
		{
			GEODETIC_INSTANT pt;
			pt = m_engine.ReferencePoint();
			m_getpropertybuffer.resize(4);
			m_getpropertybuffer.at(0) = pt.latitude;
			m_getpropertybuffer.at(1) = pt.longitude;
			m_getpropertybuffer.at(2) = pt.heightm;
			m_getpropertybuffer.at(3) = pt.mjd;
			return true;
		}
	);

	AddGetVectorFunction( "wavel",
		[&, this]( int index )
		{
			m_getpropertybuffer = m_wavelen;
			return true;
		}
	);

	AddGetVectorFunction( "sun",
		[&, this]( int index )
		{
			nxVector sun;
			sun = m_engine.InternalSpecs()->CoordinateSystemPtr()->GetSunUnitVector();
			m_getpropertybuffer.resize(3);
			m_getpropertybuffer.at(0) = sun.X();
			m_getpropertybuffer.at(1) = sun.Y();
			m_getpropertybuffer.at(2) = sun.Z();
			return true;
		}
	);

	AddGetVectorFunction( "observer",
		[&, this]( int losindex )
		{
			bool ok;
			ok = (losindex >= 0) && (losindex < (int)m_linesofsight.NumRays());
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_SO::GetPropertyArray, The Observer array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)losindex, (int)m_engine.LinesOfSight()->NumRays());
			}
			else
			{
				const SKTRAN_LineOfSightEntry_V2*	entry;
				ok = m_linesofsight.GetRay( losindex, &entry);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_SO::GetPropertyArray, The Observer array getproperty error retrieving entry for line of sight index (%d)", (int)losindex);
				}
				else
				{
					nxVector pt;

					pt = entry->Observer();
					m_getpropertybuffer.resize(3);
					m_getpropertybuffer.at(0) = pt.X();
					m_getpropertybuffer.at(1) = pt.Y();
					m_getpropertybuffer.at(2) = pt.Z();
				}
			}
			return ok;
		}
	);

	AddGetVectorFunction( "look",
		[&, this]( int losindex )
		{
			bool ok;
			ok = (losindex >= 0) && (losindex < (int)m_linesofsight.NumRays());
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_SO::GetPropertyArray, The Look array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)losindex, (int)m_engine.LinesOfSight()->NumRays());
			}
			else
			{
				const SKTRAN_LineOfSightEntry_V2*	entry;
				ok = m_linesofsight.GetRay( losindex, &entry);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_SO::GetPropertyArray, The Look array getproperty error retrieving entry for line of sight index (%d)", (int)losindex);
				}
				else
				{
					nxVector pt;

					pt = entry->Look();
					m_getpropertybuffer.resize(3);
					m_getpropertybuffer.at(0) = pt.X();
					m_getpropertybuffer.at(1) = pt.Y();
					m_getpropertybuffer.at(2) = pt.Z();
				}
			}
			return ok;
		}
	);
}
