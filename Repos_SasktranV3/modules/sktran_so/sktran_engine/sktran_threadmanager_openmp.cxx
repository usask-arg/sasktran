#include "../sasktranv21_internals.h"
#include <omp.h>


// This file must be compiled with the /openMP switch on the C++ command line.
// In visyual studio it is set on the "Project" properties of this file C/C++->Language ->OpenMP support.
// Note that pre-compiled headers have to be disabled for just this file as the openMP switch conflicts with
// the pre-compiled headers.

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerOpenMP		2011-5-31*/
/** The OpenMP multi-threaded implementation**/
/*---------------------------------------------------------------------------*/

class SKTRAN_ThreadManagerOpenMP : public SKTRAN_ThreadManager
{


	private:
		size_t										m_numthreads;
		const SKTRAN_TableOpticalProperties_V21*	m_opticalproperties;
		SKTRANSO_Quadrature_TLS_V21**					m_quadrature_tls;						// An array of pointers, one for each thread. This holds thread specific context.
		bool*										m_ok;
		size_t*										m_numpoints;

		static size_t								NumDefaultProcessingThreads	();

	private:
		void										ReleaseResources();
		bool										AllocateQuadratureStorage	( size_t numthreads );
		bool										SetOpticalPropsInThreads	();


	public:
													SKTRAN_ThreadManagerOpenMP	();
		virtual									   ~SKTRAN_ThreadManagerOpenMP	();

	public:
		virtual bool								SetOpticalProps						( const SKTRAN_TableOpticalProperties_V21* optprop);
		virtual bool								CreateThreads						( size_t numthreads, const SKTRAN_SpecsInternal_V21* modelspecifications, const SKTRAN_EngineDiffuseTables* modeltables );
		virtual bool								NotifyWorkerThreadsAndWaitForEnd	( size_t numpoints );
		virtual bool								CloseThreadsAndRelease				( ) { ReleaseResources(); return true;}
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerOpenMP_CreateNewInstance		2011-5-31*/
/** This is a global function that is used to create new instances of the
 *	BoostThreadManager
 **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadManager* SKTRAN_ThreadManagerOpenMP_CreateNewInstance()
{
	SKTRAN_ThreadManager* manager;
	
	manager = new SKTRAN_ThreadManagerOpenMP;
	return manager;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::NumDefaultProcessingThreads		2010-5-26*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_ThreadManagerOpenMP::NumDefaultProcessingThreads()
{
#if defined(NXDEBUG)
	return 1;
#else
	size_t	numcpu;
	numcpu = omp_get_max_threads();			// Get the maximum number of threads as determined by OpenMP
	if (numcpu == 0) numcpu = 1;
	if (numcpu > 150) numcpu = 150;
	return numcpu;
#endif
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::SKTRAN_ThreadManager		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadManagerOpenMP::SKTRAN_ThreadManagerOpenMP()
{
	m_numthreads = 0;
	m_opticalproperties = NULL;
	m_quadrature_tls    = NULL;
	m_ok                = NULL; 
	m_numpoints         = NULL;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::~SKTRAN_ThreadManager		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadManagerOpenMP::~SKTRAN_ThreadManagerOpenMP	()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerOpenMP::ReleaseResources		2011-5-31*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_ThreadManagerOpenMP::ReleaseResources()
{
	size_t	i;

	if (m_quadrature_tls != NULL)
	{
		for (i = 0; i < m_numthreads; i++)
		{
			if ( m_quadrature_tls[i] != NULL)  m_quadrature_tls[i]->Release();
		}
		delete [] m_quadrature_tls;
		m_quadrature_tls = NULL;
	}

	if (m_ok != NULL)        delete [] m_ok;
	if (m_numpoints != NULL) delete [] m_numpoints;

	m_ok = NULL;
	m_numpoints = NULL;
	m_numthreads = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerOpenMP::AlocateQuadratureStorage		2011-5-31*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerOpenMP::AllocateQuadratureStorage( size_t numthreads )
{
	bool	ok;

	NXASSERT(( m_quadrature_tls == NULL ));									// The object should be previously released 

	ok = (numthreads > 0);													// Make sure we are dealing with a positive integer
	if (ok)
	{
		m_quadrature_tls = new SKTRANSO_Quadrature_TLS_V21* [ numthreads ];		// Allcoate the memory
		m_ok             = new bool   [numthreads];
		m_numpoints      = new size_t [numthreads];
		ok = (m_quadrature_tls != NULL) && ( m_ok != NULL) && (m_numpoints != NULL);	
	}
	if (ok)
	{
		m_numthreads = numthreads;
		omp_set_num_threads( (int)m_numthreads );
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_ThreadManagerOpenMP::AllocateQuadratureStorage, Error creating storage for threads. Thats not good");
		ReleaseResources();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::SetOpticalProps		2010-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerOpenMP::SetOpticalProps( const SKTRAN_TableOpticalProperties_V21* optprop)
{
	bool  ok;
	
	m_opticalproperties = optprop;
	ok = SetOpticalPropsInThreads();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerOpenMP::SetOpticalPropsInThreads		2011-6-1*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerOpenMP::SetOpticalPropsInThreads()
{
	size_t	i;

	if (m_opticalproperties != NULL)
	{
		for (i = 0; i < m_numthreads; i++)
		{
			m_quadrature_tls[i]->SetOpticalProps( m_opticalproperties );
		}
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerBoost::CreateThreads		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerOpenMP::CreateThreads(   size_t								numthreads,
										         const SKTRAN_SpecsInternal_V21*		modelspecifications,
										         const SKTRAN_EngineDiffuseTables*		modeltables   )
{
	bool									ok;
	bool									ok1;
	size_t									i;
	SKTRANSO_Quadrature_TLS_V21*				quadrature;
	const SKTRAN_Quadrature_Factory_V21*		quadraturefactory;


	CloseThreadsAndRelease();																				// release any existing resources
	if (numthreads == 0) numthreads = NumDefaultProcessingThreads();
	quadraturefactory = modelspecifications->QuadratureSpecs()->QuadratureFactory();						// Get the factory object to generate thread specific quadrature objects
	ok = AllocateQuadratureStorage( numthreads );															// Allocate space  to hold the quadrature object pointers

	if (ok)
	{
		for (i = 0; i < m_numthreads; ++i)																		// Now create 
		{																										// the quadrature objects for each thread
			ok1  = quadraturefactory->CreateThreadInstance( modelspecifications, modeltables, &quadrature );	// this will return quadrature with an AddRef in place
			m_quadrature_tls[i] = quadrature;																	// assign it to our list of pointers. we will do the release when finished
			m_ok[i]             = true;
			ok = ok && ok1;																						// And get the error status
		}
	}
	ok = ok && SetOpticalPropsInThreads();

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_ThreadManagerOpenMP::CreateThreads, Error creating threads");
		ReleaseResources();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManagerBoost::NotifyWorkerThreadsAndWaitForEnd		2010-3-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManagerOpenMP::NotifyWorkerThreadsAndWaitForEnd( size_t numpoints )
{
	bool	allok;
	size_t	i;
	int		npoints = (int)numpoints;
	size_t		np;

	for ( i = 0; i < m_numthreads; ++i)
	{
		m_ok[i] = true;
		m_numpoints[i] = 0;
	}

// ----- start of OpenMP parallel do loop
	#pragma omp parallel for schedule(dynamic)
		for (int pointindex =0; pointindex < npoints; pointindex++)
		{
			size_t  threadindex;
			bool	ok;

			threadindex = omp_get_thread_num();
			NXASSERT(( threadindex < m_numthreads ));
   			ok = this->ThreadExecuteAction( m_quadrature_tls[threadindex], pointindex);				// Execute the action provided by the manager ( this method MUST BE THREAD SAFE !!!!)
			m_ok[threadindex] = m_ok[threadindex] && ok;
			m_numpoints[threadindex]++;
		}

// ----- end of OpenMP parallel do loop

	allok = true;																					// Get Status of thread actions
	np = 0;
	for ( i = 0; i < m_numthreads; ++i)
	{
		allok = allok && m_ok[i];
		np += m_numpoints[i];
	}
	NXASSERT(( np == numpoints ));
	return allok;
}
