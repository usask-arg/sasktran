#include "../sasktranv21_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::SKTRANSO_TableDiffusePoints		2007-12-13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_TableDiffusePoints::SKTRANSO_TableDiffusePoints()
{
	m_slongrid   = NULL;
	m_cosszagrid = NULL;
	m_opticaltable = NULL; //new SKTRAN_TableDiffusePointsOptical_V21(this);

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePoints_2D_Height_SZA::~SKTRAN_TableDiffusePoints_2D_Height_SZA		2007-12-13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_TableDiffusePoints::~SKTRANSO_TableDiffusePoints()
{
	ReleaseResources();
	if (m_opticaltable != NULL) delete m_opticaltable;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::ReleaseResources		2007-12-13*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_TableDiffusePoints::ReleaseResources()
{
	DeleteProfiles();
	if (m_slongrid   != NULL) m_slongrid->Release();
	if (m_cosszagrid != NULL ) m_cosszagrid->Release();
	m_slongrid   = NULL;
	m_cosszagrid = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::DeleteProfiles		2010-2-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::DeleteProfiles()
{
	size_t idx;

	for (idx = 0; idx < m_profiles.size(); ++idx)
	{
		delete m_profiles[idx];
	}
	m_profiles.clear();
	m_pointsindex.clear();
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::AllocateProfiles		2007-12-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::AllocateProfiles( size_t numprofiles )
{
	bool									ok;
	bool									ok1;
	SKTRAN_DiffuseProfileGeometry_V21*		entry;
	size_t									idx;

	ok = (m_profiles.size() == numprofiles);

	if (!ok)
	{
		DeleteProfiles();
		m_profiles.reserve   (numprofiles);
		m_pointsindex.reserve(numprofiles+1);
		ok = true;
		for (idx = 0; idx < numprofiles; idx++)
		{
			entry = new SKTRAN_DiffuseProfileGeometry_V21;
			ok1 = (entry != NULL);
			if (ok1)
			{
				m_profiles.push_back(entry);
				m_pointsindex.push_back(0);		// TODO: fix bug here: upper_bound generates error
			}
			ok = ok && ok1;
			//m_pointsindex.push_back(0);
		}
		m_pointsindex.push_back(0);	
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_TableDiffusePoints_2D_Height_SZA::AllocateProfiles, error allocating space for %Iu diffuse profiles", (size_t)numprofiles);
			DeleteProfiles();
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::FetchGrids		2010-3-29*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::FetchSLONandSZAGrids(const SKTRAN_SpecsInternal_Diffuse_V21*    diffusespecs )
{
	bool								ok;
	const SKTRAN_GridDefCosSZA_V21*		szagrid;
	const SKTRAN_GridDefSLON_V21*		slongrid;	

	szagrid  = diffusespecs->DiffuseProfileCosSZA();
	slongrid = diffusespecs->DiffuseProfileSLON();
	ok =       (szagrid != NULL) && (slongrid != NULL);
	ok = ok && (szagrid->NumAngles() == slongrid->NumAngles());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRANSO_TableDiffusePoints::FetchSLONandSZAGrids, The SZA and SLON are NULL or ifferent sizes, thats a problem");
	}
	else
	{
		szagrid->AddRef();
		slongrid->AddRef();
		if (m_slongrid   != NULL) m_slongrid->Release();
		if (m_cosszagrid != NULL ) m_cosszagrid->Release();
		m_slongrid   = slongrid;
		m_cosszagrid = szagrid;
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::CreateGrid		2007-12-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::CreateGrid_Stage1( std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords, 
													   const SKTRAN_SpecsInternal_RayTracing_V21* raytracingspecs,
													   const SKTRAN_SpecsInternal_Quadrature_V21* quadraturespecs,
													   const SKTRAN_SpecsInternal_Diffuse_V21*    diffusespecs )
{
	size_t									profileidx;
	size_t									ngood;
	bool									ok;
	bool									ok1;
	HELIODETIC_POINT						point;
	HELIODETIC_UNITVECTOR					direction;
	double									sinsza;
	double									cossza;
	double												r;
	double												h;
	double												slon;
	size_t												baseindex;


	h        = raytracingspecs->RayTracingShells()->LowestShell();
	r        = coords->AltitudeToRadius( h );
	ok = FetchSLONandSZAGrids(diffusespecs );											// Fetch the Longitude and latitude of each profile.
	if (ok)
	{
		ngood    = m_cosszagrid->NumAngles();
		ok       = AllocateProfiles( ngood );															// allocate the space for the profiles
		if (ok)																							// If that worked
		{																								// then
			baseindex = 0;
			for (profileidx = 0; profileidx < ngood; ++profileidx )										// repeat the above loop testing for acceptable values
			{																							// so for each potentially good value
				cossza   = m_cosszagrid->CosineSZA(profileidx);
				slon     = m_slongrid->SLON       (profileidx);
				sinsza   = sqrt(1.0-cossza*cossza);
				direction.SetCoords( sinsza*cos(slon), sinsza*sin(slon),  cossza);
				point.Initialize( direction, r, coords.get() );
				ok1	     = m_profiles[profileidx]->ConfigureGeometry_Stage1( point,						// Configure the diffuse profile
																			(SKTRAN_GridIndex)profileidx,
																			coords,
																			RayFactory(),
																			diffusespecs,
																			this  );	
				ok       = ok && ok1;																	// and flag if its ok.
				m_pointsindex[profileidx] = baseindex;
				baseindex += m_profiles[profileidx]->NumPoints();
			}																							// if sza is out of range then ignore
			m_pointsindex[ngood] = baseindex;
		}
	}
	if (!ok) ReleaseResources();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePoints_2D_Height_SZA::ConfigureGeometry		2008-1-13*/
/** Initializes the table. This normally allocates objects and creates rays. After all the tables are created
 *	we should call LinkIncomingRaysToJIndexTables to finish off the wavelength independent procesisng
 *	of the tables.
**/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::ConfigureGeometry(  const SKTRAN_SpecsInternal_V21* specs, SKTRAN_ThreadManager* threadmanager )
{
	const	SKTRAN_SpecsInternal_RayTracing_V21*	raytracingspecs;
	const	SKTRAN_SpecsInternal_Quadrature_V21*	quadraturespecs;
	const	SKTRAN_SpecsInternal_Diffuse_V21*		diffusespecs;
	bool	ok;

	raytracingspecs = specs->RayTracingSpecs();
	quadraturespecs = specs->QuadratureSpecs();
	diffusespecs    = specs->DiffuseSpecs();

	ok =	   SetRayFactory( specs->RayFactoryScattered() );
	ok = ok && CreateGrid_Stage1(  specs->CoordinateSystemObject(), raytracingspecs, quadraturespecs, diffusespecs );
	ok = ok && threadmanager->DiffusePointsTable_ConfigureGeometryStage2	( this );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePoints_2D_Height_SZA::ConfigureGeometry_Stage2MT		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRANSO_TableDiffusePoints::ConfigureGeometry_Stage2MT	( size_t pointindex, SKTRANSO_Quadrature_TLS_V21* threadquadrature )
{
	SKTRAN_DiffusePointGeometry_V21* point;
	bool							ok;

	ok =       LookupDiffusePoint	( pointindex, &point );
	ok = ok && point->ConfigureGeometry_Stage2MT(threadquadrature);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::LookupDiffusePoint		2010-3-29*/
/** Lookup a diffuse point given the index of the point. The index given
 *	is the same as the one given by the table. This is used by the multi-threaded processing
 *	when processing the diffuse points table.
**/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::LookupDiffusePoint( size_t pointindex, SKTRAN_DiffusePointGeometry_V21** point )
{
	size_t							profileindex;
	size_t							heightindex;
	bool							ok;

	ok = PointIndexToProfileAndHeight( pointindex,	&profileindex, &heightindex );
	if (ok)  *point  = m_profiles[profileindex]->AtVar( heightindex );
	else     *point  = NULL;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::LookupDiffusePoint		2010-3-29*/
/** Lookup a diffuse point given the index of the point. The index given
 *	is the same as the one given by the table. This is used by the multi-threaded processing
 *	when processing the diffuse points table.
**/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::PointIndexToProfileAndHeight( size_t pointindex,	size_t* profileindex, size_t* heightindex ) const
{
	std::vector<size_t>::const_iterator	iter;
	bool								ok;

	iter = std::upper_bound( m_pointsindex.begin(), m_pointsindex.end(), pointindex  );		// Find the pointer to the value greater than x
	ok = !(iter == m_pointsindex.end());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_TableDiffusePoints::PointIndexToProfileAndHeight, Error looking up pointindex, That should not happen, its a problem");
		*profileindex = 0;
		*heightindex  = 0;
	}
	else
	{
		*profileindex = (iter - m_pointsindex.begin()) - 1;
		*heightindex  =  pointindex - m_pointsindex[*profileindex];
		NXASSERT(( (*profileindex < m_profiles.size()) && (*heightindex <= m_profiles[*profileindex]->NumPoints()) ));
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::CreateJIndexTables		2010-3-26*/
/**	Create the JIndex tables for the diffuse points in the table. This can only be done after all the tables
 *	in the engine have been initialized with Configure
**/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::CreateJIndexTables_RayIntegral( SKTRAN_ThreadManager* threadmanager  )
{
	bool	ok;

	ok = threadmanager->DiffusePointsTable_CreateJIndices( this );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::ConvertJIndexToRadiancePtr		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_StokesScalar* SKTRANSO_TableDiffusePoints::ConvertJIndexToRadiancePtr( const SKTRANSO_JIndex* entry, ENUM_SKTRAN_JSOURCE jsource ) const
{
	size_t						profileindex;
	size_t						altitudeindex;
	size_t						rayindex;
	const SKTRAN_StokesScalar*	JPtr;

	profileindex  = entry->PositionIndex();
	altitudeindex = entry->HeightIndex();
	rayindex      = entry->VertexIndex();
	JPtr          = m_opticaltable->PointAt( profileindex, altitudeindex)->OutgoingJPtr( rayindex );
	return JPtr;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileGeometry_V21::PointAt		2010-3-31*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_DiffusePointGeometry_V21* SKTRANSO_TableDiffusePoints::PointAt( size_t profileidx, size_t radiusindex ) const
{
	return m_profiles[profileidx]->At(radiusindex);
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::AttachOpticalToGeometry		2010-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::AttachOpticalToGeometry( )
{
	NXASSERT( (m_opticaltable == NULL ));
	m_opticaltable = new SKTRAN_TableDiffusePointsOptical_V21(this);

	return m_opticaltable->AttachToGeometry( );
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::ConfigureOptical		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableDiffusePoints::ConfigureOptical( bool singlescatter, const SKTRAN_TableOpticalProperties_V21*	optprop, SKTRAN_ThreadManager *	threadmanager )
{
	bool	ok;

	ok =  m_opticaltable->ConfigureOptical( singlescatter, optprop, threadmanager );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_TableDiffusePoints::ConfigureOptical, error configuring the wavelength dependent diffuse points table. Thats a problem");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::ScatterIncomingRadiance		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRANSO_TableDiffusePoints::ScatterIncomingRadiance( bool ignorehighaltdiffuse )
{
	return m_opticaltable->ScatterIncomingRadiance(ignorehighaltdiffuse);
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::InitializeFirstOrderIncomingRadiances		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRANSO_TableDiffusePoints::InitializeFirstOrderIncomingRadiances()
{
	return m_opticaltable->InitializeFirstOrderIncomingRadiances();
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::CalculateIncomingRadiance		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRANSO_TableDiffusePoints::CalculateIncomingRadiance( bool ignorehighaltdiffuse )
{
	return m_opticaltable->CalculateIncomingRadiance(ignorehighaltdiffuse);
}
*/
