#include <skclimatology21.h>

/*-----------------------------------------------------------------------------
 *					skClimatology_Constant::skClimatology_Constant					2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_Constant::skClimatology_Constant( )
{
	m_constantvalue = 1.0;
}

skClimatology_Constant::skClimatology_Constant(double value )
{
	m_constantvalue = value;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Constant::~skClimatology_Constant					2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_Constant::~skClimatology_Constant( )
{
}

/*-----------------------------------------------------------------------------
 *					skClimatology_Constant::UpdateCache		2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Constant::UpdateCache( const GEODETIC_INSTANT& placeandtime )
{
	SetCacheIsLoaded(true );
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Constant::GetParameter			2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Constant::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	*value = m_constantvalue;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Constant::CreateClone		2008-12-16*/
/** **/
/*---------------------------------------------------------------------------*/

//bool skClimatology_Constant::CreateClone( skClimatology** userclone) const
//{
//	bool	ok;
//	skClimatology_Constant*	clone;
//
//	clone = new skClimatology_Constant;
//
//	ok = (clone != NULL );
//	ok = ok && (clone)->AddRef();
//	if (ok)
//	{
//		(clone)->m_constantvalue = m_constantvalue;
//	}
//	ok = ok && (clone)->DeepCopy(*this);
//	*userclone = clone;
//	return ok;
//}

/*-----------------------------------------------------------------------------
 *					skClimatology_Constant::IsSupportedSpecies		2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Constant::IsSupportedSpecies(const CLIMATOLOGY_HANDLE& species )
{
	return true;
}

/*-----------------------------------------------------------------------------
 *					skClimatology_OneTemperatureAndPressure::skClimatology_OneTemperatureAndPressure		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/


skClimatology_OneTemperatureAndPressure::skClimatology_OneTemperatureAndPressure()
{
	m_T = 273.15;
	m_P = 101000.0;
	SetCacheIsLoaded(true );
}


/*-----------------------------------------------------------------------------
 *					skClimatology_OneTemperatureAndPressure::~skClimatology_OneTemperatureAndPressure		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_OneTemperatureAndPressure::~skClimatology_OneTemperatureAndPressure()
{
}



/*-----------------------------------------------------------------------------
 *					skClimatology_OneTemperatureAndPressure::UpdateCache		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OneTemperatureAndPressure::UpdateCache( const GEODETIC_INSTANT& placeandtime )
{
	SetCacheIsLoaded(true );
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_OneTemperatureAndPressure::GetParameter		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OneTemperatureAndPressure::GetParameter( const CLIMATOLOGY_HANDLE& species, const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{

	bool	ok = true;

	if (species == SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3)
	{
		*value = m_P/( nxSI::KBOLTZMAN * m_T )*1.0E-06;	
	}
	else if (species == SKCLIMATOLOGY_PRESSURE_PA)
	{
		*value = m_P;
	}
	else if (species == SKCLIMATOLOGY_TEMPERATURE_K)
	{
		*value = m_T;
	}
	else
	{
		*value = MissingValue();
		ok     = false;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_OneTemperatureAndPressure::IsSupportedSpecies		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OneTemperatureAndPressure::IsSupportedSpecies  ( const CLIMATOLOGY_HANDLE& species )
{
	bool	ok;

	ok  =    (species == SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3)
		  || (species == SKCLIMATOLOGY_PRESSURE_PA)
		  || (species == SKCLIMATOLOGY_TEMPERATURE_K);
	return ok;
}

