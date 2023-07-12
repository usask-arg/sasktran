#include <skclimatology21.h>

/*-----------------------------------------------------------------------------
 *					skClimatology_Zero::skClimatology_Zero					2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_Zero::skClimatology_Zero( )
{
}

/*-----------------------------------------------------------------------------
 *					skClimatology_Zero::~skClimatology_Zero					2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_Zero::~skClimatology_Zero( )
{
}

/*-----------------------------------------------------------------------------
 *					skClimatology_Zero::UpdateCache		2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Zero::UpdateCache( const GEODETIC_INSTANT& placeandtime )
{
	SetCacheIsLoaded(true );
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Zero::GetParameter			2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Zero::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	*value = 0.0;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Zero::CreateClone		2008-12-16*/
/** **/
/*---------------------------------------------------------------------------*/

//bool skClimatology_Zero::CreateClone( skClimatology** clone) const
//{
//	bool	ok;
//
//	*clone = new skClimatology_Zero;
//
//	ok = (*clone != NULL );
//	ok = ok && (*clone)->AddRef();
//	ok = ok && (*clone)->DeepCopy(*this);
//	return ok;
//}
//
/*-----------------------------------------------------------------------------
 *					skClimatology_Zero::IsSupportedSpecies		2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Zero::IsSupportedSpecies(const CLIMATOLOGY_HANDLE& species )
{
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Undefined::skClimatology_Undefined					2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_Undefined::skClimatology_Undefined( )
{
}

/*-----------------------------------------------------------------------------
 *					skClimatology_Undefined::~skClimatology_Undefined					2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_Undefined::~skClimatology_Undefined( )
{
}

/*-----------------------------------------------------------------------------
 *					skClimatology_Undefined::UpdateCache		2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Undefined::UpdateCache( const GEODETIC_INSTANT& placeandtime )
{
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Undefined::GetParameter			2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Undefined::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	nxLog::Record(NXLOG_WARNING,"skClimatology_Undefined::GetParameter, You are using the Undefined climatology. This class is usually used to remined you that you have forgot to set  climatology or atmopsheric state");
	*value = 0.0;
	return true;
}

//
/*-----------------------------------------------------------------------------
 *					skClimatology_Undefined::IsSupportedSpecies		2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Undefined::IsSupportedSpecies(const CLIMATOLOGY_HANDLE& species )
{
	return true;
}

