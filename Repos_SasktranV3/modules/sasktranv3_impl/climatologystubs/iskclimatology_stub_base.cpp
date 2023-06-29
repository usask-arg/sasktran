#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"

/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_Base::ISKClimatology_Stub_Base		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_Base::ISKClimatology_Stub_Base	( skClimatology* climate)
{
	m_climatology = climate;
	if (m_climatology != NULL) m_climatology->AddRef();
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_Base::~ISKClimatology_Stub_Base		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_Base::~ISKClimatology_Stub_Base	()
{
	if (m_climatology != NULL) m_climatology->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_Base::SkOpticalPropertyPointer		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/
//
//nxUnknown*	ISKClimatology_Stub_Base::RawObjectPointer()
//{ 
//	return m_climatology;
//}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_Base::UpdateCache		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology_Stub_Base::UpdateCache( const GEODETIC_INSTANT& location )
{
	return m_climatology->UpdateCache( location );
}

/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_Base::GetParameter		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology_Stub_Base::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& location, double* value )
{
	return m_climatology->GetParameter( species, location, value, false);
}

/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_Base::Set		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology_Stub_Base::SetPropertyUserDefined( const CLIMATOLOGY_HANDLE& species,  double* profile, int numpoints )
{
	nxLog::Record(NXLOG_WARNING,"ISKClimatology_Stub_Base::Set, this object does not support user defined profiles \n");
	return false;
}

