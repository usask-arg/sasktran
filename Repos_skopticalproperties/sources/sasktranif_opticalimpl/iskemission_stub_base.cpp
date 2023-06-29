#include <skopticalproperties21.h>
#include "sources/sasktranif_opticalimpl/skoptprop_stubs.h"
#include "sources/sasktranif_opticalimpl/skemission_stubs.h"


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Base::ISKEmission_Stub_Base		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEmission_Stub_Base::ISKEmission_Stub_Base	( skEmission* optprop)
{
	m_emission = optprop;
	if (m_emission != NULL) m_emission->AddRef();
}


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Base::~ISKEmission_Stub_Base		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEmission_Stub_Base::~ISKEmission_Stub_Base	()
{
	if (m_emission != NULL) m_emission->Release();
}

/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Base::UpdateLocation		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission_Stub_Base::UpdateLocation( const GEODETIC_INSTANT& pt, bool isground )
{
	return m_emission->UpdateLocation( pt, isground );
}


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Base::InternalClimatology_UpdateCache		 2014- 10- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission_Stub_Base::UpdateCache( const GEODETIC_INSTANT& pt )
{
	return m_emission->UpdateCache(pt);
}


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Base::GetExtinctionPerCm		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission_Stub_Base::IsotropicEmission( double wavenumber, double* isotropicradiance)
{
	bool	ok;

	ok = m_emission->IsotropicEmission( wavenumber, isotropicradiance);
	return ok;
}



/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Base::CalculateCrossSectionsArray		2014-3-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission_Stub_Base::IsotropicEmissionArray( const double * wavenumber, int numwavenumber, double *isorad, int numisorad)
{
	bool	ok;
	std::vector<double>	w;
	std::vector<double>	rad;

	w.assign( wavenumber, wavenumber+numwavenumber);
	ok =	m_emission->IsotropicEmissionArray( w, &rad);
	std::copy( rad.begin(), rad.end(), isorad);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Tabulated_HeightWavelength		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEmission_Stub_Tabulated_HeightWavelength::ISKEmission_Stub_Tabulated_HeightWavelength( skEmission_Tabulated_HeightWavelength* emission)
	                            :ISKEmission_Stub_Base(emission)
{
	m_userdefinedemission = emission;
	MakeSetPropertyFunctions();

}


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Tabulated_HeightWavelength::~ISKEmission_Stub_Tabulated_HeightWavelength		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEmission_Stub_Tabulated_HeightWavelength::~ISKEmission_Stub_Tabulated_HeightWavelength	()
{
}


/*---------------------------------------------------------------------------
 * ISKEmission_Stub_Tabulated_HeightWavelength::MakeSetPropertyFunctions2020-08-24 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKEmission_Stub_Tabulated_HeightWavelength::MakeSetPropertyFunctions()
{
	AddSetVectorFunction( "heights",
		[&, this](const double* value, int numpoints )
		{
			nx1dArray<double>  values;
			values.Attach( numpoints, (double *)value);
			m_currentheightarray = values;
			return true;
		}
	);
	AddSetVectorFunction( "wavelengths",
		[&, this](const double* value, int numpoints )
		{
			nx1dArray<double>  values;
			values.Attach( numpoints, (double *)value);
			m_currentwavelengtharray = values;
			return true;
		}
	);
	AddSetVectorFunction( "emissiontable",
		[&, this](const double* value, int numpoints )
		{
			nx2dArray<double>				emit;
			size_t							numy;
			size_t							numx;
			bool							ok;

			numy = m_currentheightarray.size();
			numx = m_currentwavelengtharray.size();

			ok = (numy*numx) == numpoints;
			ok = ok && emit.Attach(numx,numy, (double *)value);
			ok = ok && m_userdefinedemission->SetEmissionTable( emit, m_currentwavelengtharray,m_currentheightarray);

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEmission_Stub_Tabulated_HeightWavelength::SetPropertyArray, Error setting the emission table. make sure you have previously set Array Property Heights and Wavelengths");
			}
			return ok;
		}
	);
}


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Thermal						2017-8-30*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEmission_Stub_Thermal::ISKEmission_Stub_Thermal(skEmission_Thermal* emission)
							:ISKEmission_Stub_Base(emission)
{
	m_thermalemission = emission;
	MakeSetPropertyFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Thermal::~ISKEmission_Stub_Thermal		2017-8-30*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEmission_Stub_Thermal::~ISKEmission_Stub_Thermal()
{
}


/*---------------------------------------------------------------------------
 *       ISKEmission_Stub_Thermal::MakeSetPropertyFunctions       2020-08-21 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKEmission_Stub_Thermal::MakeSetPropertyFunctions()
{

	AddSetScalarFunction( "emissivity",
		[&, this](double value )
		{
			m_thermalemission->SetGroundEmissivity(value);
			return true;
		}
	);

	AddSetObjectFunction("atmospheric_state",
		[&, this](nxUnknown* object )
		{
			skClimatology * climate = ( skClimatology*) object;
			return m_thermalemission->SetAtmosphericState( climate );
		}
	);
}


/*-----------------------------------------------------------------------------
 *					:ISKEmission_Stub_HitranChemical		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEmission_Stub_HitranChemical::ISKEmission_Stub_HitranChemical	( skEmission_HitranChemical* emission)
	                           :ISKEmission_Stub_Base (emission)
{
	m_thermalemission = emission;
	MakeSetPropertyFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_HitranChemical::~ISKEmission_Stub_HitranChemical		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEmission_Stub_HitranChemical::~ISKEmission_Stub_HitranChemical	()
{
}

/*---------------------------------------------------------------------------
 *   ISKEmission_Stub_HitranChemical::MakeSetPropertyFunctions    2020-08-21 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKEmission_Stub_HitranChemical::MakeSetPropertyFunctions()
{

	AddSetStringFunction( "set_chemical_name",
		[&, this](const char *value )
		{
			bool				ok;
				
			ok  = m_thermalemission->SetChemicalName( value );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEmission_Stub_HitranChemical::SetProperty(set_chemical_name) failed for name <%s>.", (const char*)value);
			}
			return ok;
		}
	);


	AddSetVectorFunction( "set_wavenumber_range",
		[&, this] ( const double* value, int numpoints)
		{
			bool ok;
			ok =      ( numpoints == 2); 
			ok = ok && m_thermalemission->SetWavenumberRange( value[0], value[1] );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKEmission_Stub_HitranChemical::SetProperty(setwavenumberrange). There were errors setting the array properties.");
			}
			return ok;
		}
	);

	AddSetVectorFunction( "enable_cached_emissions",
		[&, this] ( const double* value, int numpoints)
		{
			bool ok;
			double *	wavenumptr = (double *)(intptr_t)value;
			ok = m_thermalemission->EnableCachedEmissions( wavenumptr, numpoints );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKEmission_Stub_HitranChemical::SetProperty(enablecachedcrosssections). There were errors setting the array properties.");
			}
			return ok;
		}
	);

	AddSetScalarFunction( "setisotopefilter",
		[&, this](double value )
		{
			bool	ok;
			int	filterid;

			filterid = (int)value;
			ok = m_thermalemission->SetIsotopeIdFilter( filterid );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKEmission_Stub_HitranChemical::SetProperty(setisotopefilter). There were errors setting the property of SetIsotopeFilter to %15.8e",(double)value);
			}
			return ok;
		}
	);
	
	AddSetStringFunction( "set_lower_state_global_quanta_filter",
		[&, this](const char *value )
		{
			bool				ok;
				
			ok  = m_thermalemission->SetLowerStateGlobalQuantaFilter( std::string(value) );

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEmission_Stub_HitranChemical::SetProperty(set_lower_state_global_quanta_filter) failed.");
			}
			return ok;
		}
	);

	AddSetStringFunction( "set_upper_state_global_quanta_filter",
		[&, this](const char *value )
		{
			bool				ok;
				
			ok  = m_thermalemission->SetUpperStateGlobalQuantaFilter( std::string(value) );

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEmission_Stub_HitranChemical::SetProperty(set_upper_state_global_quanta_filter) failed.");
			}
			return ok;
		}
	);

	AddSetObjectFunction( "set_self_broadening_climatology", 
		[&, this]( nxUnknown* object )
		{
			bool			ok = false;
			skClimatology*	climate;

			climate = dynamic_cast<skClimatology*>(object);
			ok      = m_thermalemission->SetSelfBroadeningClimatology( climate );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEmission_Stub_HitranChemical::SetProperty( set_self_broadening_climatology ) failed");
			}
			return ok;
		}
	);

	AddSetStringFunction( "set_self_broadening_climatology_handle",
		[&, this](const char *value )
		{
			bool				ok;
			CLIMATOLOGY_HANDLE*	handleptr;
				
			handleptr = FindGlobalClimatologyHandle ( value );					// MAke sure this handle already exists. If it does not then ** IN THE FUTURE *** we should perhaps create the entry
			ok  = (handleptr != nullptr);
			ok  = ok && m_thermalemission->SetSelfBroadeningClimatologyHandle( *handleptr );

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEmission_Stub_HitranChemical::SetProperty(set_self_broadening_climatology_handle) failed as <%s> was not a recognised existing climatology handle.", (const char*)value);
			}
			return ok;
		}
	);	

	AddSetObjectFunction( "set_excited_upper_state_climatology", 
		[&, this]( nxUnknown* object )
		{
			bool			ok = false;
			skClimatology*	climate;

			climate = dynamic_cast<skClimatology*>(object);
			ok      = m_thermalemission->SetUpperStateNumberDensity( climate );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEmission_Stub_HitranChemical::SetProperty( set_excited_upper_state_climatology ) failed");
			}
			return ok;
		}
	);

	AddSetStringFunction( "set_excited_upper_state_climatology_handle",
		[&, this](const char *value )
		{
			bool				ok;
			CLIMATOLOGY_HANDLE*	handleptr;
				
			handleptr = FindGlobalClimatologyHandle ( value );					// MAke sure this handle already exists. If it does not then ** IN THE FUTURE *** we should perhaps create the entry
			ok  = (handleptr != nullptr);
			ok  = ok && m_thermalemission->SetUpperStateNumberDensityHandle( *handleptr );

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEmission_Stub_HitranChemical::SetProperty(set_excited_upper_state_climatology_handle) failed as <%s> was not a recognised existing climatology handle.", (const char*)value);
			}
			return ok;
		}
	);

	AddSetObjectFunction( "set_atmospheric_state_climatology", 
		[&, this]( nxUnknown* object )
		{
			bool			ok = false;
			skClimatology*	climate;

			climate = dynamic_cast<skClimatology*>(object);
			ok      = m_thermalemission->SetAtmosphericState( climate );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKEmission_Stub_HitranChemical::SetProperty( set_atmospheric_state_climatology ) failed");
			}
			return ok;
		}
	);
}

