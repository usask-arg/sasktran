#include <skopticalproperties21.h>

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::skOpticalProperties_HenyeyGreenstein		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_HenyeyGreenstein::skOpticalProperties_HenyeyGreenstein()
{
	m_extinction = NULL;
	m_g          = 0.7;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::~skOpticalProperties_HenyeyGreenstein		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_HenyeyGreenstein::~skOpticalProperties_HenyeyGreenstein()
{
	if (m_extinction != NULL) m_extinction->Release();
	m_extinction = NULL;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::SetCrossSectionObject		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HenyeyGreenstein::SetCrossSectionObject( skOpticalProperties* object )
{
	if (object != NULL) object->AddRef();
	if (m_extinction != NULL) m_extinction->Release();
	m_extinction  = object;
	if (m_extinction != NULL)
	{
		if (!m_extinction->IsScatterer())
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HenyeyGreenstein::SetCrossSectionObject, You are setting the cross-sections of the Henyey Greenstein object with an object that does not scatter!!. This will have no effect on the calculation");
		}
	}
	return true;
}




/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::Set_HG_AsymmetryFactor		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HenyeyGreenstein::Set_HG_AsymmetryFactor( double g )
{
	m_g = g;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::DeepCopy		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_HenyeyGreenstein::DeepCopy( const skOpticalProperties_HenyeyGreenstein& other )
{
	bool	ok;

	if (m_extinction != NULL ) m_extinction->Release();						// release the current extinction object
	m_extinction = NULL;													// and set it to NULL

	ok = (other.m_extinction == NULL);										// everything is good if other objects extinction object is NULL
	if (!ok) ok = other.m_extinction->CreateClone( &m_extinction );			// otherwise clone the object

	m_g = other.m_g;														// copy over the asymmetry factor
	return ok;																// and that is that
}
*/
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::SetAtmosphericState		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HenyeyGreenstein::SetAtmosphericState( skClimatology* neutralatmosphere)
{
	bool	ok;

	ok = (m_extinction != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_HenyeyGreenstein::SetAtmosphericState, You must define the extinction object, before calling this function");
	}
	else
	{
		ok = m_extinction->SetAtmosphericState( neutralatmosphere);
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "skOpticalProperties_HenyeyGreenstein::SetAtmosphericState, Error setting atmospheric state of associated extinction object");
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::SetAtmosphericState		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HenyeyGreenstein::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged)
{
	bool	ok;

	ok = (m_extinction != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_HenyeyGreenstein::SetAtmosphericState, You must define the extinction object, before calling this function");
	}
	else
	{
		ok = m_extinction->SetLocation( pt, crosssectionschanged );
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "skOpticalProperties_HenyeyGreenstein::SetAtmosphericState, Error setting atmospheric state of associated extinction object");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::InternalClimatology_UpdateCache		2011-8-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HenyeyGreenstein::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt)
{
	bool	ok;

	ok = (m_extinction != NULL);
	if (ok)
	{
		ok = m_extinction->InternalClimatology_UpdateCache( pt );
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "skOpticalProperties_HenyeyGreenstein::InternalClimatology_UpdateCache, Error setting atmospheric state of associated extinction object");
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::IsScatterer		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HenyeyGreenstein::IsScatterer() const
{
	bool	ok;

	ok = (m_extinction != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_HenyeyGreenstein::IsScatterer, You must define the extinction object, before calling this function");
	}
	else
	{
		ok = m_extinction->IsScatterer();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::IsAbsorber		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HenyeyGreenstein::IsAbsorber() const
{
	bool	ok;

	ok = (m_extinction != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_HenyeyGreenstein::IsAbsorber, You must define the extinction object, before calling this function");
	}
	else
	{
		ok = m_extinction->IsAbsorber();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::CalculateCrossSections		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HenyeyGreenstein::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs)
{
	bool	ok;

	ok = (m_extinction != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_HenyeyGreenstein::CalculateCrossSections, You must define the extinction object, before calling this function");
	}
	else
	{
		ok = m_extinction->CalculateCrossSections( wavenumber, absxs, extxs, scattxs);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::CreateClone		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_HenyeyGreenstein::CreateClone ( skOpticalProperties** userclone) const
{
	skOpticalProperties_HenyeyGreenstein*	clone;
	bool								ok;

	clone = new skOpticalProperties_HenyeyGreenstein;
	ok = (clone != NULL);
	if (ok)
	{
		clone->AddRef();
		ok = ok && clone->DeepCopy(*this);
	}
	*userclone = clone;
	return ok;
}

*/


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HenyeyGreenstein::CalculatePhaseMatrix		2008-11-17*/
/** Claculates the Henyey Greenstein phase matrix. Normailized so it equals 4pi
*	when integrated over the unit sphere **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HenyeyGreenstein::CalculatePhaseMatrix( double /*wavenumber*/, double cosscatterangle, skRTPhaseMatrix* P)
{
	double	s_ph;
	double	mg2;

	mg2  = m_g*m_g;
	s_ph = (1.0 - mg2)/ pow( (1.0 + mg2 - 2.0 * m_g * cosscatterangle), 1.5);

	P->At(1,1) = (SKRTFLOAT)( s_ph);
	P->At(1,2) = (SKRTFLOAT)( 0 );
	P->At(1,3) = (SKRTFLOAT)( 0 );
	P->At(1,4) = (SKRTFLOAT)( 0 );
	
	P->At(2,1) = (SKRTFLOAT)( 0 );
	P->At(2,2) = (SKRTFLOAT)( 0 );
	P->At(2,3) = (SKRTFLOAT)( 0 );
	P->At(2,4) = (SKRTFLOAT)( 0 );

	P->At(3,1) = (SKRTFLOAT)( 0 );
	P->At(3,2) = (SKRTFLOAT)( 0 );
	P->At(3,3) = (SKRTFLOAT)( 0 );
	P->At(3,4) = (SKRTFLOAT)( 0 );

	P->At(4,1) = (SKRTFLOAT)( 0 );
	P->At(4,2) = (SKRTFLOAT)( 0 );
	P->At(4,3) = (SKRTFLOAT)( 0 );
	P->At(4,4) = (SKRTFLOAT)( 0 );
	return true;
} 




