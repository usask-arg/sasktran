#include "sasktranif_internals.h"



static CLIMATOLOGY_HANDLE g_userdefinedhandlebase = { 0x00000000, 0x0000, 0x0002, 0x02, 0x02, 0x02, 0x04, 0x02, 0x02, 0x02, 0x02};

/*-----------------------------------------------------------------------------
 *					SASKTRANIF_Create_NewUserDefined_ClimatologyName		 2016- 7- 8*/
/** Creates a new Climatology Name that is recognized by the climatology interfaces.**/
/*---------------------------------------------------------------------------*/

bool ISKClimatology::Create_New_ClimatologyName( const char* name )
{
	bool	ok; 

	ok = HasKey_InGlobalClimatologyHandle( name);
	if (!ok)
	{
		g_userdefinedhandlebase.Data1++;
		ok = AddGlobalClimatologyHandle ( name, g_userdefinedhandlebase);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::ISKEngine		2014-2-8*/
/** Create an instance of the specified engine. Current values for name
 *	are:
 *	- SO, successive orders engine
 *  - HR, High resolution successive orders engine.
 *  - MC, Monte-Carlo engine
 *  - OCC, Occultation engine.
 *	The constructor immediately locates and loads the engine stub from its
 *	DLL.
 **/
/*---------------------------------------------------------------------------*/

ISKClimatology::ISKClimatology( const char* name )
{
	SasktranIF_ClassFactoryLocator	classfactory;

	classfactory.CreateISKClimatology( name, &m_climatology, DllNamePtr() );
}



/*-----------------------------------------------------------------------------
 *					ISKClimatology::~ISKClimatology		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology::~ISKClimatology()
{
	if (m_climatology != NULL) m_climatology->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology::SkClimatologyPointer		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

nxUnknown* ISKClimatology::RawObjectUnknown	()
{
	return (m_climatology != NULL) ? m_climatology->RawObjectPointer() : NULL;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology::UpdateCache		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology::UpdateCache( const GEODETIC_INSTANT& location )
{
	return (m_climatology != NULL) ?  m_climatology->UpdateCache( location ) : false;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology::GetParameter		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology::GetParameter( const char* speciesname,  const GEODETIC_INSTANT& location, double* value )
{
	nxString	name(speciesname);

	name.RemoveWhiteSpace();
	const CLIMATOLOGY_HANDLE& species = *FindGlobalClimatologyHandle(name);

	return (m_climatology != NULL) ?  m_climatology->GetParameter( species, location, value ) : false;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology::GetHeightProfile		 2016- 9- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology::GetHeightProfile( const char* speciesname,  const GEODETIC_INSTANT& userlocation, const double* altitude, double *profile, int numalt )
{
	bool	ok = true;
	bool	ok1;
	size_t	i;
	size_t	n = numalt;
	GEODETIC_INSTANT	location(userlocation);

	for (i = 0; i < n; i++)
	{
		location.heightm = altitude[i];
		ok1 = GetParameter( speciesname,  location, &profile[i] );
		ok = ok && ok1;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology::SetPropertyScalar		2014-2-9*/
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

bool ISKClimatology::SetPropertyScalar( const char* propertyname, double value )
{
	bool	ok;

	ok = (m_climatology != NULL) && m_climatology->SetPropertyScalar( propertyname, value);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology::SetPropertyArray		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
{
	bool	ok;

	ok = (m_climatology != NULL) && m_climatology->SetPropertyArray( propertyname, value, numpoints);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology::SetPropertyObject		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology::SetPropertyObject( const char* propertyname, ISKModuleBase* object )
{
	bool	ok;

	ok = (m_climatology != NULL) && m_climatology->SetPropertyObject( propertyname, object->RawObjectUnknown());
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology::SetPropertyString		 2016- 9- 26*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology::SetPropertyString( const char* propertyname, const char* str ) 
{
	bool	ok;

	ok = (m_climatology != NULL) && m_climatology->SetPropertyString( propertyname, str);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology::SetPropertyUserDefined		2014-3-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology::SetPropertyUserDefined( const char* speciesname,  double* profilevalues, int numpoints)
{
	bool	ok;
	nxString	name(speciesname);

	name.RemoveWhiteSpace();
	const CLIMATOLOGY_HANDLE& species = *FindGlobalClimatologyHandle(name);

	ok = (m_climatology != NULL) && m_climatology->SetPropertyUserDefined( species, profilevalues, numpoints);
	return ok;
}

