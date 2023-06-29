#include "../sasktranv21_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileOptical_V21::SKTRAN_DiffuseProfileOptical_V21		2008-1-4*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_DiffuseProfileOptical_V21::SKTRAN_DiffuseProfileOptical_V21	( SKTRAN_DiffuseProfileGeometry_V21*	geometry )
{
	m_geometry  = geometry;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileOptical_V21::~SKTRAN_DiffuseProfileOptical_V21		2008-1-4*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_DiffuseProfileOptical_V21::~SKTRAN_DiffuseProfileOptical_V21	()
{
	ReleaseResources();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileOptical_V21::ReleaseResources		2008-1-4*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_DiffuseProfileOptical_V21::ReleaseResources( )
{
	size_t	idx;

	for (idx = 0; idx < m_points.size(); idx++)
	{
		if (m_points[idx] != NULL )
		{
			m_points[idx]->Release();
		}
	}
	m_points.clear();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileOptical_V21::AttachToGeometry		2008-1-4*/
/** Attaches this optical, wavelength dependent profile to the geometry profile.
 *	This simply allocates the profile and initializes the points
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffuseProfileOptical_V21::AttachToGeometry( )
{
	bool							ok = true;
	bool							ok1;
	size_t							idx;
	size_t							numpoints;
	SKTRAN_DiffusePointOptical_V21*	point;

	numpoints = m_geometry->NumPoints();									// Getthe number of points in this profile
	ReleaseResources();														// release any existing points
	m_points.reserve( numpoints); 											// do the allocation
	for (idx=0; idx < numpoints; idx++)										// for each point
	{																		// in the profile
		point =  m_geometry->AtVar(idx)->OpticalPoint();
		point->AddRef();
		ok1   =  point->AttachToGeometry();									// Attach the point to the geometry
		if (ok1) m_points.push_back( point );								// and save the point in our array of points
		else     point->Release();
		ok = ok && ok1;														
	}
	if (!ok)																// if it did not work
	{																		// then tell the user
		nxLog::Record( NXLOG_WARNING, "SKTRAN_DiffuseProfileOptical_V21::AttachToGeometry, Error attaching points to geometry");
		ReleaseResources();													// if it did not work the clear this object
	}
	return ok;																// return the status
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileOptical_V21::ConfigureOptical		2008-1-4
 * I left this in to emphasize that we dont actually need it.
 *---------------------------------------------------------------------------*/

/*
bool SKTRAN_DiffuseProfileOptical_V21::ConfigureOptical()
{
	return true;
}
*/


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileOptical_V21::ScatterToOutgoingRadiation		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_DiffuseProfileOptical_V21::ScatterIncomingRadiance( bool ignorehighaltdiffuse )
{
	bool	ok1;
	bool	ok;
	size_t	idx;

	ok = true;
	for( idx = 0; idx < m_points.size(); idx++)
	{
		ok1 = m_points[idx]->ScatterIncomingRadiance( ignorehighaltdiffuse );
		ok = ok && ok1;
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileOptical_V21::CalculateIncomingRadiation		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_DiffuseProfileOptical_V21::CalculateIncomingRadiance( bool ignorehighaltdiffuse )
{
	bool	ok1;
	bool	ok;
	size_t	idx;

	ok = true;
	for( idx = 0; idx < m_points.size(); idx++)					// idx == 0 is the ground point
	{
		ok1 = m_points[idx]->CalculateIncomingRadiances( ignorehighaltdiffuse );
		ok = ok && ok1;
	}
	return ok;
}
*/


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileOptical_V21::InitializeFirstOrderIncomingRadiances		2008-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_DiffuseProfileOptical_V21::InitializeFirstOrderIncomingRadiances()
{
	bool	ok1;
	bool	ok;
	size_t	idx;

	ok = true;
	for( idx = 0; idx < m_points.size(); idx++)
	{
		ok1 = m_points[idx]->InitializeFirstOrderIncomingRadiances();
		ok = ok && ok1;
	}
	return ok;
}
*/
