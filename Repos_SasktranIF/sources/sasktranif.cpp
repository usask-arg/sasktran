#include "sasktranif_internals.h"


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

ISKEngine::ISKEngine( const char* name )
{
	SasktranIF_ClassFactoryLocator	classfactory;

	classfactory.CreateISKEngine( name, &m_engine,DllNamePtr() );
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::~ISKEngine		2014-2-8*/
/** Destroys the ISKEngine object and all of the associated resources.
**/
/*---------------------------------------------------------------------------*/

ISKEngine::~ISKEngine()
{
	if (m_engine != NULL) m_engine->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::AddLineOfSight		2014-2-8*/
/** Adds the given line of sight and mjd to the engine. The engine will calculate
 *	the radiance along each line of sight in the array for each requested
 *	wavelength at the next call to #CalculateRadiance. The lines of sight
 *	are added to an internal array. 
 *
 *  \par Geographic Extent of Lines of Sight.
 *	Most engines we provide implicitly assume that you are passing in lines
 *	of sight and mjd's that are sensible for the problem that the engine is
 *	trying to solve. In particular, most engines assume 
 *	that the array of lines of sight are collected over a relatively narrow 
 *	region of the Earth, e.g. one scan of an atmospheric limb sounding instrument.
 *	As a guideline for time duration we normally say that angular size of the sun is
 *	0.5 degrees and the Earth will rotate through this angle in 2 minutes.  Thus it
 *	makes sense to choose lines of sight collected within a two minute window for
 *	the successive order engines. Similarly a (north-south) ray at 50 km tangent
 *	altitude spans about 15 degrees of latitude or about 1/24th of an orbit or
 *	~4 minutes of time for low earth orbits. Again it makes sense to pick lines of sight
 *	collected within 2-4 minutes of each other from a low earth orbit so we have most of the
 *	lines of sight mostly within the region spanned by 1 ray tangent at 50 km altitude.
 *
  *	\par Representative Instant
 *	The engines will normally choose an average mjd and an "average" reference
 *	point on the Earth from the lines of sight you provide.
 *	The average mjd is used as the representative instant which used by the
 *	climatologies and optical properties that define the atmospheric state. 
 *	The representative instant is also used to define the position of the sun.
 *
 *	\par Reference Point
 *	The reference point, which is calculated from the lines of sight, is considered
 *	to be the point on the surface of the Earth that best represents all of the
 *	lines of sight. An osculating sphere is fitted to this point with the intent of
 *	finding the sphere that has the curvature that best matches the scenario being calculated.
 *
 *	\param mjd
 *		The modified julian date of this line of sight
 *
 *	\param observer
 *		The geographic location of the observer expressed in meters as a cartesian position vector
 *		from the center of the geoid earth. The coordinates are provided in the array in the order [X,Y,Z].
 *		Class nxGeodetic in Repos_BaseCode can be used to convert latitudes, longitudes and heights to
 *		cartesian position vectors.
 *
 *	\param lookvector
 *		The look direction of the ray away from the observer expressed as a cartesian unit vector.
 *		Note that the ray actually travels in the opposite direction towards the observer. The
 *		coordinates are provided in the array in the order [X,Y,Z].
 *
 *	\param	losindex
 *		Returns the index of this line of sight in the internal array. This number is only
 *		informational as the line of sight is always added to the end of the internal array.
 *
 *	\returns
 *		True if success otherwise false.
 **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::AddLineOfSight( double mjd,  const nxVector& observer, const nxVector& lookvector, int* losindex  )
{
	bool	ok;

	*losindex = 0;
	ok = (m_engine != NULL) && m_engine->AddLineOfSight( mjd, observer, lookvector, losindex);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::AddSpecies		2014-2-8*/
/** Adds a species to the atmospheric optical state or replaces an existing one.
 *	The atmospheric optical state is used to calculate absorption, extinction and
 *	scattering within the atmosphere. Users will typically add all of the species
 *	that have important contributions to absorption, extinction and scattering in the
 *	wavelength region of interest to them. Note that the engines internally define another
 *	atmospheric state to provide information for ray-tracing and spectral line-shape,
 *	e.g. pressure, temperature, number density, water vapour and CO2 vmr.
 *	
 *	\param species
 *		The species being added to the atmospheric optical state. The climatology
 *		must support this species.
 *
 *	\param climatology
 *		The climatology that provides support for "species". Note that the engine
 *		uses the same internal climatology structure instance as the one that is stored within the 
 *		ISKClimatology "climatology" object. Thus any changes made to the "climatology"
 *		variable will also apply in the subsequent  calls to CalculateRadiance.
 *
 *	\param opticalproperty
 *		The optical properties, (i.e. absorption, scattering and extinction), of the
 *		species being passed in. Note that the engine uses the same internal optical
 *		properties structure instance as the one that is stored within the 
 *		ISKOpticalProperty "opticalproperty" object. Thus any changes made to the "opticalproperty"
 *		variable will also apply in subsequent calls to CalculateRadiance.
 *
 *	\returns
 *		True if successful otherwise false.
**/
/*---------------------------------------------------------------------------*/

bool ISKEngine::AddSpecies( const char* speciesname, ISKClimatology& climatology, ISKOpticalProperty& opticalproperty)
{
	bool ok;
	nxString	name(speciesname);

	name.RemoveWhiteSpace();
	const CLIMATOLOGY_HANDLE& species = *FindGlobalClimatologyHandle(name);

	ok =  (m_engine != NULL) && m_engine->AddSpecies( species, climatology.Stub(), opticalproperty.Stub());
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::AddEmission		2015-3-12*/
/** Adds an emission species to the atmospheric optical state or replaces an existing one.
 *	The atmospheric optical state can be used to calculate the effects of volume emission
 *	within the atmosphere.
 *	
 *	\param speciesid
 *		The id code of the emission object being added to the atmospheric optical state.
 *
 *	\param emission
 *		The emission object being added to the engine. The emission object be added to
 *		the engines internal table
 *
 *	\returns
 *		True if successful otherwise false.
**/
/*---------------------------------------------------------------------------*/

bool ISKEngine::AddEmission( const char* speciesname, ISKEmission& emission )
{
	bool ok;
	nxString	name(speciesname);

	name.RemoveWhiteSpace();
	const EMISSION_HANDLE& speciesid = *FindGlobalClimatologyHandle(name);

	ok =  (m_engine != NULL) && m_engine->AddEmission( speciesid,  emission.Stub() );
	return ok;
}
/*-----------------------------------------------------------------------------
 *					ISKEngine::SetAtmosphericState		2014-2-9*/
/** The atmospheric state used by the engine to provide atmospheric conditions
 *	for line-shape and ray tracing calculations. The atmospheric state climatology
 *	object should provide pressure (Pascals) and temperature (Kelvins). The object
 *	is used by the optical property objects to calculate the dependency of 
 *	cross-sections on atmospheric conditions (eg. line shape and temperature dependence).
 *	It is also used by any possible curved ray tracing code to calculate refractive index parameters,
 *	which will need pressure, temperature, number density and possibly water vapour and CO2. All engines
 *	provide a default atmospheric state that is reasonable, e.g. MSIS/CIRA and this method allows
 *	users to improve upon the default value. It must be called before  
 *	the first call by this instance to CalculateRadiance.
 *
 *	\par Who provides exoctic information
 *	The atmospheric state object is only intended to provide general purpose atmopsheric 
 *	information such as pressure and temperature. Any optical properties within the engine
 *	that need more esoteric information will implement their own internal climatological
 *	property. For example, aerosol classes require climatologies of
 *	mode radius and mode width: these are provided  by an internal class specific and not by the atmospheric state.
 *
 *	\param climatology
 *		The ISKClimatology object to be used by the engine for setting atmospheric state
 *		for line-shape and ray tracing calculations 
 **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::SetAtmosphericState( ISKClimatology& climatology )
{
	bool ok;

	ok =  (m_engine != NULL) && m_engine->SetAtmosphericState( climatology.Stub() );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine::SetAlbedo		2014-2-8*/
/** Set the albedo that will be used in the next call
 *	to CalculateRadiance. The albedo is assumed to be Lambertian and constant for
 *	all wavelengths.
 *
 *	\param albedo
 *		The albedo to be used. Usually a value betweem 0 and 1. Some engines may allow
 *		you to have albedos greater than 1.0. This can be useful to offset excessive
 *		extinction in the troposphere.
**/
/*---------------------------------------------------------------------------*/

bool ISKEngine::SetAlbedo( double albedo )
{
	bool	ok;

	ok = (m_engine != NULL) && m_engine->SetAlbedo( albedo );
	return ok;
}

/*-----------------------------------------------------------------------------
*					ISKEngine::SetBRDF		2017-3-10*/
/** Set the BRDF that will be used in the next call
*	to CalculateRadiance. 
*
*	\param brdf
*		The brdf to be used. An instance of ISKBrdf
**/
/*---------------------------------------------------------------------------*/

bool ISKEngine::SetBRDF(ISKBrdf* brdf)
{
	bool	ok;

	ok = (m_engine != NULL) && m_engine->SetBRDF(brdf->Stub());
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyScalar		2014-2-9*/
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

bool ISKEngine::SetPropertyScalar( const char* propertyname, double value )
{
	bool	ok;

	ok = (m_engine != NULL) && m_engine->SetPropertyScalar( propertyname, value);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyArray		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
{
	bool	ok;

	ok = (m_engine != NULL) && m_engine->SetPropertyArray( propertyname, value, numpoints);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyObject		2014-2-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::SetPropertyObject( const char* propertyname, ISKModuleBase* object )
{
	bool	ok;

	ok = (m_engine != NULL) && m_engine->SetPropertyObject( propertyname, object->RawObjectUnknown());
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPropertyString		 2016- 9- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::SetPropertyString( const char* propertyname, const char* str ) 
{
	bool	ok;

	ok = (m_engine != NULL) && m_engine->SetPropertyString( propertyname, str);
	return ok;
}




/*-----------------------------------------------------------------------------
 *					ISKEngine::GetPropertyArray		 2014- 5- 16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::GetProperty( const char* propertyname, const double** value, int* numpoints )
{
	bool	ok;

	*value = NULL;
	*numpoints = 0;
	ok = (m_engine != NULL) && m_engine->GetProperty( propertyname, value, numpoints);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine::SetPolarizationMode		 2015- 10- 30*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::SetPolarizationMode( int polarizationmode)
{
	bool ok;

	ok =  (m_engine != NULL) && m_engine->SetPolarizationMode( polarizationmode);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine::SetWavelengths		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::SetWavelengths( const double* wavelen, int numwavelen )
{
	bool ok;

	ok =  (m_engine != NULL) && m_engine->SetWavelengths( wavelen, numwavelen );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::InitializeModel		2014-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::InitializeModel()
{
	bool ok;

	ok =  (m_engine != nullptr) && m_engine->InitializeModel();

	return ok;

}
/*-----------------------------------------------------------------------------
 *					ISKEngine::CalculateRadiance		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::CalculateRadiance(const double** radiance, int* numwavelens, int* numlinesofsight)
{
	bool ok;

	ok =  (m_engine != nullptr) && m_engine->CalculateRadiance( radiance, numwavelens, numlinesofsight );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::CalculateRadiancePolarized		 2015- 10- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::CalculateStokesVector( const ISKStokesVector** radiancep, int* numwavelens , int* numlinesofsight)
{
	bool ok;

	ok =  (m_engine != nullptr) && m_engine->CalculateStokesVector( radiancep, numwavelens, numlinesofsight );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::GetWeightingFunctions		 2015- 10- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine::GetWeightingFunctions( const double** wf, int* numwavel, int* numlinesofsight, int* numwf)
{
	bool ok;

	ok = (m_engine != nullptr) && m_engine->GetWeightingFunctions(wf, numwavel, numlinesofsight, numwf);

	return ok;
}


