#include "include/sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Manager::SKTRAN_HR_Thread_Manager		2013-06-27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Thread_Manager::SKTRAN_HR_Thread_Manager()
{
//	m_coords     = nullptr;
	m_optint     = nullptr;
	m_srcint     = nullptr;
//	m_rayfactory  = nullptr;
	m_solartrans = nullptr;
	m_opttable   = nullptr;
	m_numthreads = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Manager::~SKTRAN_HR_Thread_Manager		2013-06-27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Thread_Manager::~SKTRAN_HR_Thread_Manager()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Manager::SetNumThreads		2013-06-27*/
/** Sets the current number of threads, if numthreads is set to 0 we will 
 *  automatically use the maximum number of threads
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Thread_Manager::SetNumThreads( size_t numthreads )
{
	bool ok = true;
	
	if( numthreads > 0 )
	{
		m_numthreads = numthreads;
	}
	else
	{
#if defined(NXDEBUG)
		m_numthreads = 1;
#else
		m_numthreads = omp_get_max_threads();
#endif
	}
	
	m_threadstore.resize( m_numthreads );
	for( size_t idx = 0; idx < m_numthreads; idx++ )
	{
		ok = ok && m_threadstore[idx].Initialize( m_rayfactory.get() );
	}
	
	omp_set_num_threads( (int)m_numthreads );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Manager::ReleaseResources		2013-06-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Thread_Manager::ReleaseResources()
{
	bool ok = true;

	if( nullptr != m_optint     )  m_optint->     Release();
	if( nullptr != m_srcint     )  m_srcint->     Release();
//	if( nullptr != m_rayfactory )  m_rayfactory-> Release();
	if( nullptr != m_solartrans )  m_solartrans-> Release();
	if( nullptr != m_opttable   )  m_opttable->   Release();

	m_optint     = nullptr;
	m_srcint     = nullptr;
//	m_rayfactory = nullptr;
	m_solartrans = nullptr;
	m_opttable   = nullptr;

	m_threadstore.clear();

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Manager::Initialize		2013-06-27*/
/** Adds references for the internal objects needed by the thread manager 
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Thread_Manager::Initialize( std::shared_ptr< const SKTRAN_CoordinateTransform_V2>&        coords,
										  const SKTRAN_OpticalPropertiesIntegrator_Base&				optint,
										  const SKTRAN_SourceTermIntegrator_Base&						 srcint,
										   std::shared_ptr<const SKTRAN_RayFactory_Base>                  rayfactory,
										  const SKTRAN_SolarTransmission_Base&           solartrans,
										  const SKTRAN_TableOpticalProperties_Base&      opttable )
{
	bool ok = true;

	m_coords =     coords;

	m_optint =     &optint;
	m_optint->     AddRef();
	m_srcint =     &srcint;
	m_srcint->     AddRef();
	m_rayfactory =  rayfactory;
//	m_rayfactory->  AddRef();
	m_solartrans = &solartrans;
	m_solartrans-> AddRef();
	m_opttable =   &opttable;
	m_opttable->   AddRef();

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Manager::CreateDiffuseFirstOrder		2013-06-27*/
/** Two major tasks are done,
 *  1. The first order (single scatter) radiance is calculated for every incoming
 *     diffuse ray
 *  2. Internal diffuse indexes are calculated for every diffuse point, which are
 *     used to make the next order matrix
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Thread_Manager::CreateDiffuseFirstOrder( SKTRAN_HR_Diffuse_Table_CPU* diffusetable, double wlen )
{
	bool ok = true;
	size_t numpoints = diffusetable->NumDiffusePoints() + diffusetable->NumGroundPoints();

	// Want to parallize over every ray instead of every diffuse point, there may only be ~20 diffuse points
	// which can lead to large inefficiencies for higher core count processors
	// Since every diffuse point might have a different number of incoming rays we need to make a vector which maps 
	// each ray index to a (diffuse_point, ray) index
	std::vector<std::pair<size_t, size_t>> indexMap;
	indexMap.reserve(numpoints * diffusetable->DiffusePointAt(0).NumIncomingRays());


	for( int pointidx = 0; pointidx < numpoints; pointidx++ )
	{
		const SKTRAN_HR_Diffuse_Point& point = diffusetable->DiffusePointAt(pointidx);
		for(int rayidx = 0; rayidx < point.NumIncomingRays(); rayidx++)
		{
			indexMap.push_back(std::pair<size_t, size_t>(pointidx, rayidx));
		}
	}

	// Shuffle the indicies so they are processed in random order
	// std::random_shuffle(std::begin(indexMap), std::end(indexMap));

	diffusetable->PreSetup();
	#pragma omp parallel for schedule(dynamic, 1) reduction(&&:ok)
	for( int idx = 0; idx < (int)indexMap.size(); idx++ )
	{
		std::pair<size_t, size_t> index = indexMap[idx];
		int th_id = omp_get_thread_num();
		ok = ok && diffusetable->CalcFirstOrderIncomingRay( index.first, index.second, m_threadstore[th_id].Ray() );
	}

	ok = ok && diffusetable->DeclareFirstOrderInitialized( );

	if( diffusetable->IsSetDiagnostic( 1 ) )
		diffusetable->DumpIncomingRadiances( 1, wlen );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Manager::CreateNextOrderMatrix		2013-06-27*/
/** Creates the next order matrix from the internal diffuse indexes.  After
 *  each diffuse index is addded to the next order matrix the memory is free'd
 **/
/*---------------------------------------------------------------------------*/
/*
bool SKTRAN_HR_Thread_Manager::CreateNextOrderMatrix( SKTRAN_HR_Diffuse_Table& diffusetable )
{
	bool ok = true;

	size_t numpoints = diffusetable.NumDiffusePoints() + diffusetable.NumGroundPoints();
	size_t numrays;

	diffusetable.AllocateNextOrderMatrix();

	#pragma omp parallel for schedule(dynamic,1) private(numrays)
	for( int pointidx = 0; pointidx < numpoints; pointidx++ )
	{
		numrays = diffusetable.DiffusePointAt( pointidx ).NumIncomingRays();
		for( int rayidx = 0; rayidx < numrays; rayidx++ )
		{
			diffusetable.AddRayToNextOrderMatrix( pointidx, rayidx );
		}
	}
	diffusetable.CleanDiffuseIndexes();

	return ok;
}
*/


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Manager::CreateScatteringMatrix		2013-06-20*/
/**  Creates the sparse scattering matrix which scatterings the incoming field
 *   to the outgoing field
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Thread_Manager::CreateScatteringMatrix( SKTRAN_HR_Diffuse_Table_CPU* diffusetable )
{
	bool			ok = true;
	size_t			numpoints = diffusetable->NumDiffusePoints() + diffusetable->NumGroundPoints();



	#pragma omp parallel for schedule(dynamic, 1) reduction(&&:ok)
	for( int pointidx = 0; pointidx < (int)numpoints; pointidx++ )
	{
		ok = ok && diffusetable->CalcScatteringMatrixPoint( pointidx );
	}

	return ok;
}

bool SKTRAN_HR_Thread_Manager::ComputeFieldCPU( SKTRAN_HR_Diffuse_Table_CPU* diffusetable, size_t numorders, double wlen )
{
	bool ok = true;

	if( 1 < numorders )
	{
		size_t numDiffuseScattersToDo = numorders-1;

		ok = ok && CreateDiffuseFirstOrder( diffusetable, wlen );
		ok = ok && CreateScatteringMatrix( diffusetable );

		if( diffusetable->GetDiffuseSourceOrder()<numDiffuseScattersToDo ) ok = ok && ScatterCPU   ( diffusetable );
		if( diffusetable->IsSetDiagnostic( 1 ) ) diffusetable->DumpOutGoingRadiances( 1, wlen );
		for(size_t fullOrderCounter=0; fullOrderCounter<numDiffuseScattersToDo; ++fullOrderCounter){ // Count up to #numorders, even though algorithm should stop on its own first
			if( diffusetable->GetDiffuseSourceOrder()<numDiffuseScattersToDo ) ok = ok && NextOrderCPU ( diffusetable );
			if( diffusetable->IsSetDiagnostic( fullOrderCounter + 2 ) ) diffusetable->DumpIncomingRadiances( fullOrderCounter+2, wlen );

			if( diffusetable->GetDiffuseSourceOrder()<numDiffuseScattersToDo ) ok = ok && ScatterCPU   ( diffusetable );
			if( diffusetable->IsSetDiagnostic( fullOrderCounter + 2 ) ) diffusetable->DumpOutGoingRadiances( fullOrderCounter + 2, wlen );

		}
	}
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_HR_Thread_Manager::ComputeFieldCPU, Something went wrong computing the diffuse field.");
	return ok;
}

bool SKTRAN_HR_Thread_Manager::ScatterCPU( SKTRAN_HR_Diffuse_Table_CPU* diffusetable ) 
{
	bool ok = true;

	size_t numpoints = diffusetable->NumDiffusePoints() + diffusetable->NumGroundPoints();

	//diffusetable->SetOrder ( order );
	#pragma omp parallel for schedule(dynamic, 1) reduction(&&:ok)
	for( int pointidx = 0; pointidx < (int)numpoints; pointidx++ )
	{
		ok = ok && diffusetable->ScatterPoint( pointidx );
	}

	ok = ok && diffusetable->DeclareAllScattered ( );

	return ok;
}

bool SKTRAN_HR_Thread_Manager::NextOrderCPU( SKTRAN_HR_Diffuse_Table_CPU* diffusetable )
{
	bool ok = true;

	size_t numpoints = diffusetable->NumDiffusePoints() + diffusetable->NumGroundPoints();

	#pragma omp parallel for schedule(dynamic, 1) reduction(&&:ok)
	for( int pointidx = 0; pointidx < (int)numpoints; pointidx++ )
	{
		ok = ok && diffusetable->ComputeNextOrderPoint( pointidx );
	}

	ok = ok && diffusetable->DeclareAllIntegrated ( );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Storage::Initialize		2013-06-20*/
/** Allocates the memory needed for the internal thread storage
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Thread_Storage::Initialize( const SKTRAN_RayFactory_Base* rayfactory )
{
	bool ok;

	ok = rayfactory->CreateRayObject( &m_rayopt );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Storage::SKTRAN_HR_Thread_Storage		2013-06-27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Thread_Storage::SKTRAN_HR_Thread_Storage()
{
}

SKTRAN_RayOptical_Base* SKTRAN_HR_Thread_Storage::Ray()
{
	return m_rayopt.get();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Storage::~SKTRAN_HR_Thread_Storage		2013-06-27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Thread_Storage::~SKTRAN_HR_Thread_Storage()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Thread_Storage::ReleaseResources		2013-06-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Thread_Storage::ReleaseResources()
{
	return true;
}
