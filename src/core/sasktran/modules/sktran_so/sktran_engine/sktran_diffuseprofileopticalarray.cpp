#include "../sasktranv21_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::SKTRAN_TableDiffusePointsOptical_V21		2007-12-18*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableDiffusePointsOptical_V21::SKTRAN_TableDiffusePointsOptical_V21(SKTRANSO_TableDiffusePoints*	geometry  )
{
	m_geometry  = geometry;
}	



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::~SKTRAN_TableDiffusePointsOptical_V21		2007-12-18*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableDiffusePointsOptical_V21::~SKTRAN_TableDiffusePointsOptical_V21( )
{
	ReleaseDiffuseArray();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::ReleaseDiffuseArray		2007-12-18*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableDiffusePointsOptical_V21::ReleaseDiffuseArray()
{
	size_t	idx;

	for (idx = 0; idx < m_profiles.size(); idx++)
	{
		delete m_profiles[idx];
	}
	m_profiles.clear();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::AllocateDiffuseArrays		2007-12-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableDiffusePointsOptical_V21::AllocateDiffuseArrays ()
{
	bool								ok = true;
	size_t								numdiffuse;
	size_t								idx;
	SKTRAN_DiffuseProfileOptical_V21*	entry;
	bool								ok1;

	ReleaseDiffuseArray();															// Release the current array of diffuse profiles
	numdiffuse = m_geometry->NumProfiles();											// Get the number of diffuse profiles in teh geometry section
	m_profiles.reserve(numdiffuse);													// reserve this number of profiles
	for (idx = 0; idx < numdiffuse; idx++)											// then for each profile
	{																				// in the diffuse geometry
		ok1 = m_geometry->ProfileAtVar(idx)->CreateOpticalProfile( &entry );			// create an optical counterpart
		if (ok1) m_profiles.push_back(entry );										// and if ok save the entry
		ok = ok && ok1;																// and check the status
	}
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "SKTRAN_TableDiffusePointsOptical_V21::AllocateDiffuseArrays, error allocating the diffuse optical profiles");
		ReleaseDiffuseArray();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::AttachToGeometry		2007-12-18*/
/** Attach this object to the diffuse geometry table.  This allows it to
 *	allocate storage and otehr aspect that are not wavelength dependent
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableDiffusePointsOptical_V21::AttachToGeometry( )
{
	bool									ok;
	bool									ok1;
	size_t									idx;

	ok  = AllocateDiffuseArrays();
	if (ok)
	{
		for (idx = 0; idx <m_profiles.size(); idx++)
		{
			ok1 = m_profiles[idx]->AttachToGeometry( );
			ok = ok && ok1;
		}
	}
	if (!ok)
	{
			nxLog::Record( NXLOG_WARNING, "SKTRAN_TableDiffusePointsOptical_V21::AttachToGeometry, Error configuring the diffuse profile array");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::ConfigureOptical		2008-1-4*/
/** Configure the optical properties of the Diffuse profile table. This
 *	should only be called after a call to AttachToGeometry.  The basic philosophy
 *	is to make a one time call to AttachToGeometry for each active thread in the
 *	process and make multiple calls to ConfigureOptical for each wavelength
 *	processed by each thread.  This strategy limits the amount of memory
 *	allocation and deallocation made by the system.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableDiffusePointsOptical_V21::ConfigureOptical(bool singlescatter, const SKTRAN_TableOpticalProperties_V21*	optProp, SKTRAN_ThreadManager* threadmanager )
{
	bool	ok = true;

	ok = ok && threadmanager->DiffusePointsTable_ConfigureOpticalStage2( this, optProp );	// Now configure the individual points
	if (!ok)																				// and if there is an error
	{																						// tell the user and quit.
		nxLog::Record( NXLOG_WARNING, "SKTRAN_TableDiffusePointsOptical_V21::ConfigureOptical, Error configuring the diffuse profile array");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::ScatterToOutgoingRadiation		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_TableDiffusePointsOptical_V21::ScatterIncomingRadiance( bool ignorehighaltdiffuse)
{
	bool	ok = true;
	bool	ok1;
	size_t	idx;
	
	NXTRACE_ONCEONLY(firsttime,("SKTRAN_TableDiffusePointsOptical_V21::ScatterIncomingRadiance, This must be moved into the multi-thread code\n"));

	for (idx = 0; idx < m_profiles.size(); idx++)					// for each of 
	{																	// the diffuse profiles				
		ok1 = m_profiles[idx]->ScatterIncomingRadiance(ignorehighaltdiffuse);
		ok  = ok && ok1;												// and keep track of the status
	}																	// do all of the profiles in the diffuse table
	if (!ok)															// and if there is an error
	{																	// tell the user and quit.
		nxLog::Record( NXLOG_WARNING, "SKTRAN_TableDiffusePointsOptical_V21::ScatterToOutgoingRadiation, Error scattering incoming to outgoing profiles");
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::InitializeFirstOrderIncomingRadiances		2008-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_TableDiffusePointsOptical_V21::InitializeFirstOrderIncomingRadiances()
{
	bool	ok = true;
	bool	ok1;
	size_t	idx;
																		// and
	NXTRACE_ONCEONLY(firsttime,("SKTRAN_TableDiffusePointsOptical_V21::InitializeFirstOrderIncomingRadiances, This must be moved into the multi-thread code\n"));

	for (idx = 0; idx < m_profiles.size(); idx++)					// for each of 
	{																	// the diffuse profiles				
		ok1 = m_profiles[idx]->InitializeFirstOrderIncomingRadiances();
		ok  = ok && ok1;												// and keep track of the status
	}																	// do all of the profiles in the diffuse table
	if (!ok)															// and if there is an error
	{																	// tell the user and quit.
		nxLog::Record( NXLOG_WARNING, "SKTRAN_TableDiffusePointsOptical_V21::InitializeFirstOrderIncomingRadiances, Error initializing incoming first order radiance");
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::CalculateIncomingRadiation		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_TableDiffusePointsOptical_V21::CalculateIncomingRadiance(bool ignorehighaltdiffuse)
{
	bool	ok = true;
	bool	ok1;
	size_t	idx;
																				// and
	NXTRACE_ONCEONLY(firsttime,("SKTRAN_TableDiffusePointsOptical_V21::CalculateIncomingRadiance, This must be moved into the multi-thread code\n"));
	for (idx = 0; idx < m_profiles.size(); idx++)							// for each of 
	{																			// the diffuse profiles				
		ok1 = m_profiles[idx]->CalculateIncomingRadiance(ignorehighaltdiffuse);
		ok  = ok && ok1;														// and keep track of the status
	}																			// do all of the profiles in the diffuse table
	if (!ok)																	// and if there is an error
	{																			// tell the user and quit.
		nxLog::Record( NXLOG_WARNING, "SKTRAN_TableDiffusePointsOptical_V21::ScatterToOutgoingRadiation, Error scattering incoming to outgoing profiles");
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::LookUpJindex		2008-2-9*/
/** **/
/*---------------------------------------------------------------------------*/

/* const SKTRAN_StokesScalar* SKTRAN_TableDiffusePointsOptical_V21::LookUpJindex( const SKTRANSO_JIndex* index, ENUM_SKTRAN_JSOURCE tablesource ) const
{
	const SKTRAN_DiffusePointOptical_V21*	point;
	const SKTRAN_StokesScalar*				radiance = NULL;

	point = PointAt( index->PositionIndex(), index->HeightIndex() );
	NXASSERT(( (point != NULL) && ( (tablesource == SKTRAN_JDIFFUSE) || ( tablesource == SKTRAN_JGROUND_UPFLUX)) ));
	if (tablesource == SKTRAN_JDIFFUSE ) radiance =  point->OutgoingJPtr( index->VertexIndex() );
	else                                 radiance =  point->IncomingRadianceArray() + index->VertexIndex();
	return radiance;
}

*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::NumPoints		2010-3-31*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_TableDiffusePointsOptical_V21::NumDiffusePoints() const
{

	return m_geometry->NumDiffusePoints();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::LookupDiffusePoint		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableDiffusePointsOptical_V21::LookupDiffusePoint( size_t pointindex, SKTRAN_DiffusePointOptical_V21** point )
{
	size_t							profileindex;
	size_t							heightindex;
	bool							ok;

	ok = m_geometry->PointIndexToProfileAndHeight( pointindex,	&profileindex, &heightindex );
	if (ok)  *point  = m_profiles[profileindex]->AtVar( heightindex );
	else     *point  = NULL;
	return ok;
}
