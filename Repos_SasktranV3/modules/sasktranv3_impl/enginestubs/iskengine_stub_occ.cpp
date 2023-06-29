#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::ISKEngine_Stub_OCC		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_OCC::ISKEngine_Stub_OCC()
{
	m_updateclimatology  = true;
	m_isconfigured       = false;
	m_extinction.SetReuseMemory(true);
	MakeScalarSetFunctions();
	MakeVectorSetFunctions();
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::~ISKEngine_Stub_OCC		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_OCC::~ISKEngine_Stub_OCC()
{
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::AddLineOfSight		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::AddLineOfSight( double mjd,  const nxVector& obs, const nxVector& look, int* losindex  )
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
	m_extinction.erase();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::AddSpecies		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::AddSpecies( const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty)
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
 *					ISKEngine_Stub_OCC::AddEmission		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::AddEmission( const EMISSION_HANDLE& species, ISKEmission_Stub* emissionobject)
{
	skEmission*				emission;
	bool						ok;

	emission = dynamic_cast<skEmission*>( emissionobject->RawObjectPointer() );
	ok      = m_opticalstate.AddEmission( species, emission );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::SetAlbedo		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::SetAlbedo( double albedo )
{
	return m_opticalstate.SetAlbedo( albedo );
}

/*-----------------------------------------------------------------------------
*					ISKEngine_Stub_HR::SetBRDF		2017-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::SetBRDF(ISKBrdf_Stub* brdf)
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
 *					ISKEngine_Stub_OCC::SetPolarizationMode		 2015- 10- 30*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::SetPolarizationMode( int mode )
{
	bool ok;

	ok = (mode == 0);
	if (!ok) nxLog::Record(NXLOG_WARNING," ISKEngine OCC, The OCC engine only supports the scalar mode for polarization settings");
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::SetWavelengths		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::SetWavelengths( const double* wavelen, int numwavelen )
{
	bool ok;

	m_wavenumber.assign( wavelen, wavelen + numwavelen );
	std::for_each( m_wavenumber.begin(), m_wavenumber.end(), [](double& value) { value = 1.0E07/value;});

	m_extinction.erase();
	ok = (numwavelen ==0) || ContainerIsAscendingOrder( m_wavenumber.begin(), m_wavenumber.end() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKEngine(OCC)::SetWavelengths. The array of wavelengths should be in descending order to make the code work efficiently");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::InitializeModel		 2015- 11- 25*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::InitializeModel()
{
	if (!m_isconfigured)
	{
		m_isconfigured = m_engine.ConfigureModel( m_specs, m_linesofsight, 0 );
		
	}
	return m_isconfigured;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::CalculateRadiance		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::CalculateRadiance( const double** extinction, int* numwavelens, int* numlinesofsight )
{
	bool					ok;
	size_t					numrays = m_linesofsight.NumRays() + m_specs.TangentPointLinesOfSight().size();
	size_t					numwave = m_wavenumber.size();

	ok = InitializeModel();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_OCC::CalculateRadiance, Error initializing model");
	}
	else
	{
		m_extinction.SetSize(numrays, numwave);
		ok = m_engine.CalculateMultiWavelengthRadiance( &m_extinctionbuffer, m_wavenumber, 0, &m_opticalstate, nullptr, m_updateclimatology );
	}
	if (ok)
	{
		for (size_t iw = 0; iw <numwave; iw++)
		{
			for (size_t ih = 0; ih < numrays; ih++)
			{
				m_extinction.At( ih, iw ) = m_extinctionbuffer.at(iw).at(ih);
			}
		}
	}
	else
	{
		m_extinction.erase();
	}
	m_updateclimatology = false;
	*extinction      = m_extinction.ArrayBasePtr();
	*numwavelens     = (int)numwave;
	*numlinesofsight = (int)numrays;
	return ok;
}



/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::CalculateRadiancePolarized		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::CalculateStokesVector( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight)
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
		ok = ok && m_extinctionpolarized.SetSize(numrays,numwave); 
		if (ok)
		{
			for (size_t i = 0; i <numwave; i++)
			{
				for (size_t l = 0; l < numrays; l++)
				{
					iquv.I = m_extinction.At( l, i );
					iquv.Q = 0.0;
					iquv.U = 0.0;
					iquv.V = 0.0;
					m_extinctionpolarized.At(l,i).Assign( iquv, basdir);
				}
			}
			*radiancep = m_extinctionpolarized.ArrayBasePtr();
		}
	}
	if (!ok) *radiancep = nullptr;
	return ok;
}




/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::GetWeightingFunctions		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::GetWeightingFunctions(const double** wf, int* numwavel, int* numlinesofsight, int* numwf)
{
	bool ok = false;

	nxLog::Record(NXLOG_WARNING, "ISKEngine_Stub_OCC::GetWeightingFunctions, The OCC engine does not support calculation of weighting functions");

	*wf = nullptr;
	*numwavel = 0;
	*numlinesofsight = 0;
	*numwf = 0;

	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::Radiance		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

//double ISKEngine_Stub_OCC::Radiance( int userlosindex, int userwavelenindex)
//{
//	bool	ok;
//	size_t	losindex = userlosindex;
//	size_t	windex   = userwavelenindex;
//	double  value;
//
//	ok = (windex < m_extinction.YSize()) && ( losindex < m_extinction.XSize());
//	value = ok ? m_extinction.At(losindex,windex) : std::numeric_limits<double>::quiet_NaN();
//	return value;
//};


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::SetAtmosphericState		2014-02-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_OCC::SetAtmosphericState( ISKClimatology_Stub* climatology )
{
	bool	ok;
	ok = m_opticalstate.SetAtmosphericStateModel( dynamic_cast<skClimatology*>(climatology->RawObjectPointer()) );
	return ok;
}


/*---------------------------------------------------------------------------
 *           ISKEngine_Stub_OCC::MakeScalarSetFunctions           2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/
void ISKEngine_Stub_OCC::MakeScalarSetFunctions() 
{
	AddSetScalarFunction( "SetReferencePoint_TargetAltitude",
		[&, this](double value)
		{ 
			return m_specs.RayTracingRegionManagerVar()->SetReferencePointTargetAltitude( value);
		}
	);

	AddSetScalarFunction( "SetReferencePoint_TargetRange",
		[&, this](double value)
		{ 
			return m_specs.RayTracingRegionManagerVar()->SetReferencePointTargetRange(value);
		}
	);

	AddSetScalarFunction( "SetGroundAltitude",
		[&, this](double value)
		{ 
			return m_specs.RayTracingRegionManagerVar()->SetGroundAltitude( value);
		}
	);

	AddSetScalarFunction( "SetUpperBoundAltitude",
		[&, this](double value)
		{ 
			return m_specs.RayTracingRegionManagerVar()->SetUpperBoundAltitude(value);
		}
	);

	AddSetScalarFunction( "SetLowerBoundAltitude",
		[&, this](double value)
		{ 
			return m_specs.RayTracingRegionManagerVar()->SetLowerBoundAltitude( value);
		}
	);

	AddSetScalarFunction( "SetRayTracingWaveNumber",
		[&, this](double value)
		{ 
			return m_engine.SetRayTracingWavenumber(value);
		}
	);
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_OCC::SetPropertyArray		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKEngine_Stub_OCC::MakeVectorSetFunctions()
{
	AddSetVectorFunction( "SetReferencePoint",
		[&, this](const double* value, int n) 
		{
			bool ok;
			ok = (n == 4);
			if (ok)
			{
				ok = m_specs.SetReferencePointManually( value[0], value[1],  value[3] );		// We dont need the height variable
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(OCC)::SetPropertyArray(SetReferencePoint), The reference must be specied as a 4 element double array [latitude, longitude, height, mjd]");
			}
			return ok;
		}
	);

	AddSetVectorFunction( "SetSun",
		[&, this](const double* value, int numpoints) 
		{
			bool ok;
			ok = (numpoints == 3);
			if (ok)
			{
				nxVector sunpos( value[0], value[1], value[2]);
				ok = m_specs.SetSunManually( sunpos );
			}
			if (!ok) 
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(OCC)::SetPropertyArray(SetReferencePoint), The sun unit vector must be specied as a 3 element double array [X, Y, Z]");
			}
			return ok;
		}
	);

	AddSetVectorFunction( "AddLinesOfSightFromTangentAlt",
		[&, this](const double* value, int numpoints) 
		{
			int n = numpoints/2;
			bool ok  = true;
			bool ok1;

			for (int i = 0; i < n; i++)
			{
				ok1 = m_specs.AddLineOfSightFromTangentAlt( value[2*i], value[2*i+1] );	
				ok = ok && ok1;
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(OCC)::SetProperty(AddLinesOfSightFromTangentAlt), there were errors adding the (%d) line of sight entries ", (int)n );
			}
			return ok;
		}
	);

	AddSetVectorFunction( "SetRayTracingShells",
		[&, this](const double* value, int numpoints) 
		{
			bool ok;
			std::vector<double>	shellalts;

			shellalts.assign( value, value+numpoints);
			ok = m_specs.ConfigureUserDefinedShells(shellalts);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEngine(OCC):::SetPropertyArray(SetRayTracingShells), Error configuring the user defined shell altitudes");
			}
			return ok;
		}
	);

}

