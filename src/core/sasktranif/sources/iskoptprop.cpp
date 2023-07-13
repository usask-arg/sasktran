#include "sasktranif_internals.h"


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::ISKOpticalProperty		 2016- 9- 27*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty::ISKOpticalProperty( )
{
	m_optprop = nullptr;
}

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::ISKOpticalProperty		2014-2-8*/
/**  **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty::ISKOpticalProperty( const char* name )
{
	SasktranIF_ClassFactoryLocator	classfactory;

	classfactory.CreateISKOpticalProperty( name, &m_optprop, DllNamePtr() );
}

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::~ISKOpticalProperty		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

ISKOpticalProperty::~ISKOpticalProperty					()
{
	if (m_optprop != NULL) m_optprop->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::SkOpticalPropertyPointer		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

nxUnknown* ISKOpticalProperty::RawObjectUnknown() 
{
	return m_optprop->RawObjectPointer();
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::SetAtmosphericState		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty::SetAtmosphericState ( ISKClimatology& atmosphere)
{
	return m_optprop->SetAtmosphericState( atmosphere.Stub() );
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::SetLocation		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty::SetLocation ( const GEODETIC_INSTANT& pt )
{
	return m_optprop->SetLocation( pt );
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::InternalClimatology_UpdateCache		 2014- 10- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt)
{
	return m_optprop->InternalClimatology_UpdateCache( pt );
}

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::CalculateCrossSections		2014-3-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty::CalculateCrossSections( const double * wavenumber,  double *absxs, double* extxs, double* scattxs, int numortype)
{
	bool ok = false;

	if ( numortype >= 0)
	{
		ok = m_optprop->CalculateCrossSectionsArray( wavenumber, numortype, absxs, extxs, scattxs );
	}
	else if (numortype == -1)
	{
		ok = m_optprop->CalculateCrossSections( *wavenumber, absxs, extxs, scattxs);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKOpticalProperty::CalculateCrossSections, error invoking ISKOpticalProperty::CalculateCrossSections with numortype = %d", (int)numortype);
	}
	return ok;
}

 bool ISKOpticalProperty::CalculatePhaseMatrix(double wavenumber, double cosscatterangle, double phasematrix[16])
 {
	 return m_optprop->CalculatePhaseMatrix(&wavenumber, &cosscatterangle, phasematrix);
 }

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::SetPropertyScalar		2014-2-9*/
/** Set engine specific scalar properties. Each engine provides
 *	additional properties that can be set by the user. Documentation for
 *	each engine describes what properties are available. The desired property 
 *	is identified by a string and the value of the property is passed as floating point
 *	value.
 *
 *	\param propertyname
 *		The name of the property to be set
 *
 *	\param value
 *		The scalar value to assign to the property.
 *
 *	\returns
 *		True if successful otherwise false.
**/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty::SetPropertyScalar( const char* propertyname, double value )
{
	return m_optprop->SetPropertyScalar( propertyname, value);
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyArray		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
{
	return m_optprop->SetPropertyArray( propertyname, value, numpoints);
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyObject		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty::SetPropertyObject( const char* propertyname, ISKModuleBase* object )
{
	return m_optprop->SetPropertyObject( propertyname, object->RawObjectUnknown());
}



/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::SetPropertyString		 2016- 9- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty::SetPropertyString( const char* propertyname, const char* str ) 
{
	bool	ok;

	ok = (m_optprop != NULL) && m_optprop->SetPropertyString( propertyname, str);
	return ok;
}



/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::AddUserDefined		2014-4-2*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKOpticalProperty::AddUserDefined( double temperature,  double* wavelen_nm, int numwave, double* crosssection, int numcross)
{
	return m_optprop->AddUserDefined( temperature, wavelen_nm, numwave, crosssection, numcross);
}

bool ISKOpticalProperty::AddUserDefinedPressure(double* pressure, int numpressure, double* temperature, int numtemperature, double* wavelen_nm, int numwavel, double* crosssection, int numcross, double broadnervmr)
{
	return m_optprop->AddUserDefinedPressure(pressure, numpressure, temperature, numtemperature, wavelen_nm, numwavel, crosssection, numcross, broadnervmr);
}



/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::ISKSolarSpectrum		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKSolarSpectrum::ISKSolarSpectrum( const char* name )
{
	SasktranIF_ClassFactoryLocator	classfactory;

	classfactory.CreateISKSolarSpectrum( name, &m_solarspectrum, DllNamePtr() );
}

/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::~ISKSolarSpectrum		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKSolarSpectrum::~ISKSolarSpectrum					()
{
	if (m_solarspectrum != NULL) m_solarspectrum->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty::SkOpticalPropertyPointer		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

nxUnknown* ISKSolarSpectrum::RawObjectUnknown() 
{
	return m_solarspectrum->RawObjectPointer();
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::Irradiance		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/


bool ISKSolarSpectrum::Irradiance( const double* wavelen_nm_vacuum, double* irradiance, int numpoints )
{
	bool ok;
	if      (numpoints == -1) ok = m_solarspectrum->Irradiance     ( *wavelen_nm_vacuum, irradiance);
	else if (numpoints >=  0) ok = m_solarspectrum->IrradianceArray( wavelen_nm_vacuum, irradiance, numpoints);
	else
	{
		nxLog::Record(NXLOG_WARNING,"ISKSolarSpectrum::Irradiance, number of points (%d) is invalid", (int)numpoints);
		ok = false;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::IrradianceAt1AU		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::IrradianceAt1AU( const double* wavelen_nm_vacuum, double* irradiance, int numpoints )
{
	bool ok;
	if      (numpoints == -1) ok = m_solarspectrum->IrradianceAt1AU     ( *wavelen_nm_vacuum, irradiance);
	else if (numpoints >=  0) ok = m_solarspectrum->IrradianceAt1AUArray( wavelen_nm_vacuum, irradiance, numpoints);
	else
	{
		nxLog::Record(NXLOG_WARNING,"ISKSolarSpectrum::IrradianceAt1AU, number of points (%d) is invalid", (int)numpoints);
		ok = false;
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::SetSolarDistanceFromMjd		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::SetSolarDistanceFromMjd( double mjd )
{
	return m_solarspectrum->SetSolarDistanceFromMjd( mjd );
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::SetSolarDistanceFromAU		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::SetSolarDistanceFromAU( double au )
{
	return m_solarspectrum->SetSolarDistanceFromAU(au);
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::MinValidWavelength		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::MinValidWavelength( double* minwavelength_nm)
{
	return m_solarspectrum->MinValidWavelength( minwavelength_nm );
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::MaxValidWavelength		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::MaxValidWavelength( double* minwavelength_nm)
{
	return m_solarspectrum->MaxValidWavelength( minwavelength_nm );
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::NanometerResolutionFWHM		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::NanometerResolutionFWHM( const double* wavelen_nm_vacuum, double* resolution_nm_fwhm, int numpoints)
{
		bool ok;
	if      (numpoints == -1) ok = m_solarspectrum->NanometerResolutionFWHM     ( *wavelen_nm_vacuum, resolution_nm_fwhm);
	else if (numpoints >=  0) ok = m_solarspectrum->NanometerResolutionFWHMArray( wavelen_nm_vacuum, resolution_nm_fwhm, numpoints);
	else
	{
		nxLog::Record(NXLOG_WARNING,"ISKSolarSpectrum::NanometerResolutionFWHM, number of points (%d) is invalid", (int)numpoints);
		ok = false;
	}
	return ok;

}

/*---------------------------------------------------------------------------
 *                ISKSolarSpectrum::SampleSpacing                 2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/
bool ISKSolarSpectrum::SampleSpacing ( const double* wavelen_nm_vacuum, double* sample_spacing, int numpoints    )
{
	bool ok;
	if      (numpoints == -1) ok = m_solarspectrum->SampleSpacing     ( *wavelen_nm_vacuum, sample_spacing);
	else if (numpoints >=  0) ok = m_solarspectrum->SampleSpacingArray( wavelen_nm_vacuum,  sample_spacing, numpoints);
	else
	{
		nxLog::Record(NXLOG_WARNING,"ISKSolarSpectrum::NanometerResolutionFWHM, number of points (%d) is invalid", (int)numpoints);
		ok = false;
	}
	return ok;

}

/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::SetPropertyScalar		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::SetPropertyScalar( const char* propertyname, double value )
{
	return m_solarspectrum->SetPropertyScalar( propertyname, value );
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::SetPropertyArray		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
{
	return m_solarspectrum->SetPropertyArray( propertyname, value, numpoints );
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::SetPropertyObject		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::SetPropertyObject( const char* propertyname, ISKModuleBase* object )
{
	return m_solarspectrum->SetPropertyObject( propertyname, object->RawObjectUnknown() );
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum::SetPropertyString		 2016- 9- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum::SetPropertyString( const char* propertyname, const char* str ) 
{
	bool	ok;

	ok = (m_solarspectrum != NULL) && m_solarspectrum->SetPropertyString( propertyname, str);
	return ok;
}




