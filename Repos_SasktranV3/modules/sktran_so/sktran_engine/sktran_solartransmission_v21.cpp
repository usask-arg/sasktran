#include "../sasktranv21_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionProfile_V21::SKTRAN_TableSolarTransmissionProfile_V21		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableSolarTransmissionProfile_V21::SKTRAN_TableSolarTransmissionProfile_V21()
{
	m_cossza  = -9999.0;
	m_slon    = -9999.0;
//	m_heights = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionProfile_V21::~SKTRAN_TableSolarTransmissionProfile_V21		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableSolarTransmissionProfile_V21::~SKTRAN_TableSolarTransmissionProfile_V21()
{
	ReleaseResources();

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionProfile_V21::SKTRAN_TableSolarTransmissionProfile_V21		 2014- 11- 24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableSolarTransmissionProfile_V21::SKTRAN_TableSolarTransmissionProfile_V21	( const SKTRAN_TableSolarTransmissionProfile_V21& moveother)
{
	*this = moveother; //std::move(moveother);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionProfile_V21::operator=		 2014- 11- 24*/
/** The MOVE assignment operator. Required by std::vector when holding unique_ptr objects**/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableSolarTransmissionProfile_V21::operator=( const SKTRAN_TableSolarTransmissionProfile_V21& moveother)
{
	m_cossza = moveother.m_cossza;
	m_slon   = moveother.m_slon;
	m_heights = moveother.m_heights;
	if (moveother.m_rays.size() > 0)
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_TableSolarTransmissionProfile_V21::assignment operator, cannot assign profiles with rays defined as its not yet supported");
		throw("SKTRAN_TableSolarTransmissionProfile_V21::assignment operator error");
	}
	m_rays.resize(0);
	//m_rays.resize( moveother.m_rays.size() );
	//for ( size_t i = 0; i < m_rays.size(); i++)
	//{
	//	m_rays.at(i) = moveother.m_rays.at(i); // std::move( moveother.m_rays.at(i) );
	//}
	//moveother.m_heights = nullptr;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionProfile_V21::ReleaseResources		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableSolarTransmissionProfile_V21::ReleaseResources()
{
//	if (m_heights != NULL) m_heights->Release();
	m_rays.clear();
//	m_heights = NULL;
	m_cossza = -9999.0;
	m_slon   = -9999.0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionProfile_V21::AllocateRays		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableSolarTransmissionProfile_V21::AllocateRays( size_t numrays )
{
	bool	ok;

	m_rays.resize(numrays);
	ok = (m_rays.size() == numrays);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableSolarTransmissionProfile_V21::AllocateRays, Error allocatimng memory for %Iu rays",(size_t)numrays);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionProfile_V21::CreateProfile		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableSolarTransmissionProfile_V21::CreateProfile( size_t profileidx, const SKTRAN_SpecsInternal_V21* specs )
{
	bool ok;

	ReleaseResources();

	m_rayfactory = specs->RayFactoryTransmissionOnly();
	m_heights  = specs->RayTracingSpecs()->RayTracingShells();
	m_cossza   = specs->DiffuseSpecs()->SolarTransmissionCosSZA()->CosineSZA( (SKTRAN_GridIndex)profileidx );
	m_slon     = specs->DiffuseSpecs()->SolarTransmissionSLON()  ->SLON     ( (SKTRAN_GridIndex)profileidx ); 
//	m_heights->AddRef();
	ok         = AllocateRays( m_heights->NumShells() );								// Allocate the points we need for this diffuse profile
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableSolarTransmissionProfile_V21::CreateProfile, Error allocating resources for solar transmission profile %u", (unsigned int)profileidx );
	}
	return ok;
};
	
/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionProfile_V21::RayAtVar		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_RayInternalGeometry* SKTRAN_TableSolarTransmissionProfile_V21::RayAtVar( size_t hidx )
{
//	NXASSERT((hidx < m_numrays));
	return &m_rays.at(hidx);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionProfile_V21::RayAt		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRANSO_RayInternalGeometry* SKTRAN_TableSolarTransmissionProfile_V21::RayAt( size_t hidx ) const
{
//	NXASSERT((hidx < m_numrays));
	return &m_rays.at(hidx);
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::SKTRANSO_TableSolarTransmission		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_TableSolarTransmission::SKTRANSO_TableSolarTransmission()
{
	m_cosszagrid = NULL;
	m_slongrid   = NULL;
	m_numpoints  = 0;
	m_notrequiredforsinglescatter = true;
	m_opticaltable = new SKTRAN_TableSolarTransmissionOptical_V21( this );
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::~SKTRANSO_TableSolarTransmission		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_TableSolarTransmission::~SKTRANSO_TableSolarTransmission()
{
	delete m_opticaltable;
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::ReleaseResources		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_TableSolarTransmission::ReleaseResources()
{
	m_profiles.clear();
	m_pointsindex.clear();
	if ( m_cosszagrid != NULL ) m_cosszagrid->Release();
	if ( m_slongrid   != NULL ) m_slongrid->Release();
	m_cosszagrid = NULL;
	m_slongrid   = NULL;
	m_numpoints = 0;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::ProfileSZAandSLON		2010-5-31*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::ProfileSZAandSLON( size_t profileidx, double* cossza, double* slon) const
{
	bool	ok;

	ok = (profileidx < NumPoints());
	if (ok)
	{
		*cossza = m_profiles[profileidx].CosSZA();
		*slon   = m_profiles[profileidx].SLON();
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_TableSolarTransmission::ProfileSZAandSLON, Error looking up cossza nd SLON of profile index %u", (unsigned int)profileidx);
		*cossza = -9999.0;
		*slon   = -9999.0;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::ProfileHeights		2010-5-31*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_GridDefRayTracingShells_V21* SKTRANSO_TableSolarTransmission::ProfileHeightsPtr( size_t profileidx ) const
{
	bool										ok;
	const SKTRAN_GridDefRayTracingShells_V21*		heights;

	ok      = (profileidx < NumPoints());
	heights = ok ? m_profiles[profileidx].HeightsPtr() : NULL;
	return heights;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::FetchGrids		2010-3-29*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::FetchSLONandSZAGrids(const SKTRAN_SpecsInternal_Diffuse_V21*    diffusespecs )
{
	bool								ok;
	const SKTRAN_GridDefCosSZA_V21*		szagrid;
	const SKTRAN_GridDefSLON_V21*		slongrid;	

	szagrid  = diffusespecs->SolarTransmissionCosSZA();
	slongrid = diffusespecs->SolarTransmissionSLON();
	ok =       (szagrid != NULL) && (slongrid != NULL);
	ok = ok && (szagrid->NumAngles() == slongrid->NumAngles());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRANSO_TableDiffusePoints::FetchSLONandSZAGrids, The SZA and SLON are NULL or different sizes, thats a problem");
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
 *					SKTRANSO_TableDiffusePoints::AllocateProfiles		2007-12-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::AllocateProfiles( size_t numprofiles, const SKTRAN_SpecsInternal_V21* specs )
{
	bool	ok;
	size_t	idx;

	m_profiles.resize   (numprofiles);
	m_pointsindex.resize(numprofiles+1);
	m_numpoints = 0;

	ok = true;

	for (idx = 0; idx < numprofiles; idx++)
	{
		ok = ok && m_profiles[idx].CreateProfile( idx, specs );
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmission_2D_Height_SZA::ConfigureGeometry		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::ConfigureGeometry_Stage1( const SKTRAN_SpecsInternal_V21* specs )
{
	bool ok;

	ReleaseResources();
	m_notrequiredforsinglescatter = (specs->DiffuseSpecs()->LOSSingleScatterTableFactory() != NULL);
	ok =       FetchSLONandSZAGrids( specs->DiffuseSpecs() );
	ok = ok && AllocateProfiles( m_cosszagrid->NumAngles(), specs );
	ok = ok && ConfigureProfileIndexing( );
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "SKTRANSO_TableSolarTransmission::ConfigureGeometry_Stage1, Error configuring the solar geometry");
		ReleaseResources();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::ConfigureProfileIndexing		2010-5-31*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::ConfigureProfileIndexing()
{
	size_t	baseindex;
	size_t	profileidx;

	baseindex = 0;
	for (profileidx = 0; profileidx < m_profiles.size(); ++profileidx )							// repeat the above loop testing for acceptable values
	{																							// so for each potentially good value
		m_pointsindex[profileidx] = baseindex;
		baseindex += m_profiles[profileidx].NumRays();
	}																							// if sza is out of range then ignore
	m_pointsindex[m_profiles.size()] = baseindex;
	m_numpoints = baseindex;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::PointIndexToProfileAndHeight		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::PointIndexToProfileAndHeight	( size_t pointindex, size_t* profileindex, size_t* heightindex) const
{
	std::vector<size_t>::const_iterator	iter;
	bool								ok;

	iter = std::upper_bound( m_pointsindex.begin(), m_pointsindex.end(), pointindex  );		// Find the pointer to the value greater than x
	ok = !(iter == m_pointsindex.end());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_TableSolarTransmission::PointIndexToProfileAndHeight, Error looking up pointindex, That should not happen, its a problem");
		*profileindex = 0;
		*heightindex  = 0;
	}
	else
	{
		*profileindex = (iter - m_pointsindex.begin()) - 1;
		*heightindex  =  pointindex - m_pointsindex[*profileindex];
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::RayAtPoint		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::RayAtPoint( size_t pointindex, const SKTRANSO_RayInternalGeometry** ray ) const
{
	bool	ok;
	size_t	profileindex;
	size_t	hidx ;

	ok = PointIndexToProfileAndHeight( pointindex, &profileindex, &hidx);
	if (ok)
	{
		*ray = m_profiles[profileindex].RayAt(hidx);
	}
	else
	{
		*ray = NULL;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_StokesScalar*SKTRANSO_TableSolarTransmission::ConvertJIndexToRadiancePtr		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_StokesScalar* SKTRANSO_TableSolarTransmission::ConvertJIndexToRadiancePtr( const SKTRANSO_JIndex* entry,  ENUM_SKTRAN_JSOURCE jsource ) const
{
	size_t						pointindex;
	bool						ok;
	const SKTRAN_StokesScalar*	JPtr;

	ok   = ConvertJindexToPosition( entry, &pointindex );
	JPtr = ok ? m_opticaltable->RayAt( pointindex ): NULL;
	return JPtr;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmission_2D_Height_SZA::ConfigureGeometry_Stage2MT		2010-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::ConfigureGeometry( const SKTRAN_SpecsInternal_V21* specs, SKTRAN_ThreadManager* threadmanager )
{

	bool									ok;

	ok = SetRayFactory( specs->RayFactoryTransmissionOnly() );
	ok = ok && ConfigureGeometry_Stage1( specs );
	ok = ok && threadmanager->SolarTransmissionTable_ConfigureGeometryStage2(this, specs );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableSolarTransmission_2D_Height_SZA::ConfigureGeometry, Error configuring the table. Thats a problem" );
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmission_2D_Height_SZA::ConfigureGeometry		2007-11-16*/
/** This will process the point in the Solar Transmisison table. This simply
 *	does teh ray tracing for the ray from the point to the sun. This is called 
 *	from a multi-threaded environment.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::ConfigureGeometry_Stage2MT( size_t pointindex, SKTRANSO_Quadrature_TLS_V21* quadrature, const SKTRAN_SpecsInternal_V21* specs )

{
	bool										ok;
	size_t										profileidx;
	size_t										hidx;
	double										cossza;
	double										sinsza;
	double										slon;
	double										cosslon;
	double										sinslon;
	double										radius;
	double										height;
	HELIODETIC_VECTOR							observer;
	HELIODETIC_UNITVECTOR						look;
	SKTRANSO_RayInternalGeometry*				ray;
	const SKTRAN_SpecsInternal_RayTracing_V21*	raytracingspecs;
	const SKTRAN_CoordinateTransform_V2*		coords;

	raytracingspecs = specs->RayTracingSpecs();
	coords          = specs->CoordinateSystemPtr();

	ok = PointIndexToProfileAndHeight	( pointindex, &profileidx, &hidx) ;
	if (ok)
	{
		std::unique_ptr< SKTRAN_RayOptical_Base>	aray;
		cossza   = m_profiles[profileidx].CosSZA();
		sinsza   = sqrt(1.0-cossza*cossza);
		slon     = m_profiles[profileidx].SLON();
		cosslon  = cos(slon);
		sinslon  = sin(slon);
		height   = m_profiles[profileidx].Height(hidx);
		radius   = coords->AltitudeToRadius( height );											// Get the average/central height of the cell
		ray      = m_profiles[profileidx].RayAtVar(hidx);

		observer.SetCoords( radius*sinsza*cosslon, radius*sinsza*sinslon, radius*cossza );		// Get the radial location  of the observer
		look.SetCoords    ( 0.0, 0.0, 1.0);														// create a look vector directly towards the sun

		ok =      RayFactory()->CreateRayObject( &aray );
		ok = ok && aray->MoveObserver (observer, look );
		ok = ok && aray->TraceRay_NewMethod();
		ray->AssignRay( std::move(aray) );
	}
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "SKTRAN_TableSolarTransmission_2D_Height_SZA::ConfigureGeometry_Stage2MT, Error configuring the solar geometry");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::AttachToGeometry		2010-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::AttachOpticalToGeometry(  )
{
	return m_opticaltable->AttachToGeometry( );
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::ConfigureOptical		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::ConfigureOptical(  bool singlescatter, const SKTRAN_TableOpticalProperties_V21*	optprop, SKTRAN_ThreadManager * threadmanager)
{
	bool	ok;
	ok = threadmanager->SolarTransmissionTable_ConfigureOpticalStage2( this, optprop );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::ConfigureOptical_Stage2MT		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::ConfigureOptical_Stage2MT ( size_t pointindex, SKTRANSO_Quadrature_TLS_V21* threadquadrature )
{
	return m_opticaltable->ConfigureOptical_Stage2( pointindex, threadquadrature );
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission::ConvertJindexToPosition		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableSolarTransmission::ConvertJindexToPosition	( const SKTRANSO_JIndex* entry, size_t* positionindex ) const
{
	size_t	szaindex;
	size_t	heightindex;

	szaindex    = entry->PositionIndex();
	heightindex = entry->HeightIndex();

	NXASSERT(( (szaindex+1) < m_pointsindex.size() ));

	*positionindex = m_pointsindex[szaindex] + heightindex;

	NXASSERT(( *positionindex  < m_pointsindex[szaindex+1] ));

	return true;
}
