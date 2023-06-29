#include "sasktranif_internals.h"


/*-----------------------------------------------------------------------------
 *					ISKEmission::ISKEmission		2014-2-8*/
/**  **/
/*---------------------------------------------------------------------------*/

ISKEmission::ISKEmission( const char* name )
{
	SasktranIF_ClassFactoryLocator	classfactory;

	classfactory.CreateISKEmission( name, &m_emissionprop,DllNamePtr() );
}

/*-----------------------------------------------------------------------------
 *					ISKEmission::~ISKEmission		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEmission::~ISKEmission					()
{
	if (m_emissionprop != NULL) m_emissionprop->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKEmission::SkOpticalPropertyPointer		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

nxUnknown* ISKEmission::RawObjectUnknown() 
{
	return m_emissionprop->RawObjectPointer();
}


/*-----------------------------------------------------------------------------
 *					ISKEmission::SetAtmosphericState		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission::UpdateLocation ( const GEODETIC_INSTANT& pt, bool isground )
{
	return m_emissionprop->UpdateLocation( pt, isground );
}


/*-----------------------------------------------------------------------------
 *					ISKEmission::InternalClimatology_UpdateCache		 2014- 10- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission::UpdateCache( const GEODETIC_INSTANT& pt)
{
	return m_emissionprop->UpdateCache( pt );
}

/*-----------------------------------------------------------------------------
 *					ISKEmission::IsotropicEmission		 2015- 3- 11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission::IsotropicEmission( const double* wavenumber, double* isotropicradiance, int numortype)
{
	bool	ok = false;

	if (numortype >= 0)
	{
		ok = m_emissionprop->IsotropicEmissionArray( wavenumber, numortype, isotropicradiance, numortype );
	}
	else if (numortype == -1)
	{
		ok = m_emissionprop->IsotropicEmission( *wavenumber, isotropicradiance );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKEmission::IsotropicEmission, Error callling ISKEmission::IsotropicEmission with numortype = %d", (int)numortype);
	}
	return ok;
	
}

/*-----------------------------------------------------------------------------
 *					ISKEmission::SetPropertyScalar		2014-2-9*/
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

bool ISKEmission::SetPropertyScalar( const char* propertyname, double value )
{
	return m_emissionprop->SetPropertyScalar( propertyname, value);
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyArray		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
{
	return m_emissionprop->SetPropertyArray( propertyname, value, numpoints);
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyObject		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission::SetPropertyObject( const char* propertyname, ISKModuleBase* object )
{
	return m_emissionprop->SetPropertyObject( propertyname, object->RawObjectUnknown());
}


/*-----------------------------------------------------------------------------
 *					ISKEmission::SetPropertyString		 2016- 9- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEmission::SetPropertyString( const char* propertyname, const char* str ) 
{
	bool	ok;

	ok = (m_emissionprop != NULL) && m_emissionprop->SetPropertyString( propertyname, str);
	return ok;
}




