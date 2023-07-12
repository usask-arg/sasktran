#include "sasktranif_internals.h"


/*-----------------------------------------------------------------------------
 *					ISKBrdf::ISKBrdf		2014-2-8*/
/**  **/
/*---------------------------------------------------------------------------*/

ISKBrdf::ISKBrdf( const char* name )
{
	SasktranIF_ClassFactoryLocator	classfactory;

	classfactory.CreateISKBrdf( name, &m_brdfprop,DllNamePtr() );
}

/*-----------------------------------------------------------------------------
 *					ISKBrdf::~ISKBrdf		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf::~ISKBrdf()
{
	if (m_brdfprop != NULL) m_brdfprop->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKBrdf::SkOpticalPropertyPointer		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

nxUnknown* ISKBrdf::RawObjectUnknown() 
{
	return (m_brdfprop != nullptr) ? m_brdfprop->RawObjectPointer() : nullptr;
}



/*-----------------------------------------------------------------------------
 *					ISKBrdf::IsotropicEmission		 2015- 3- 11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf::BRDF( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* return_brdf)
{
	bool	ok = m_brdfprop != nullptr;

	ok = ok && m_brdfprop->BRDF( wavelennm, pt, MU_in, MU_out, COSDPHI, return_brdf);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKBrdf::BRDF, Error calling ISKBrdf::BRDF");
	}
	return ok;
	
}

/*-----------------------------------------------------------------------------
 *					ISKBrdf::SetPropertyScalar		2014-2-9*/
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

bool ISKBrdf::SetPropertyScalar( const char* propertyname, double value )
{
	bool	ok = m_brdfprop != nullptr;

	ok = ok && m_brdfprop->SetPropertyScalar( propertyname, value);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyArray		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
{
	bool	ok = m_brdfprop != nullptr;

	ok = ok && m_brdfprop->SetPropertyArray( propertyname, value, numpoints);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyObject		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf::SetPropertyObject( const char* propertyname, ISKModuleBase* object )
{
	bool	ok = m_brdfprop != nullptr;

	ok = ok && m_brdfprop->SetPropertyObject( propertyname, object->RawObjectUnknown());
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKBrdf::SetPropertyString		 2016- 9- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf::SetPropertyString( const char* propertyname, const char* str ) 
{
	bool	ok = m_brdfprop != nullptr;

	ok = ok && m_brdfprop->SetPropertyString( propertyname, str);
	return ok;
}




