#include "../sasktranv21_internals.h"
//#include <float.h>


/*
#include <process.h>
size_t SKTRAN_ThreadManager::NumDefaultProcessingThreads()
{
	SYSTEM_INFO		info;
	size_t			numwave;
	size_t			numproc;
	
	GetSystemInfo( &info );
	numproc = info.dwNumberOfProcessors;
	numwave = m_wavelengthmanager.NumWavelengths();
	return (numwave > numproc) ? numproc : numwave;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::SKTRAN_ThreadManager		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadManager::SKTRAN_ThreadManager()
{
	m_numpoints = 0;

	m_diffusetable_geometry     = NULL;
	m_diffusetable_optical      = NULL;
	m_solartable_geometry       = NULL;
//	m_solartable_optical        = NULL;
	m_groundpointtable_geometry = NULL;
//	m_groundpointtable_optical  = NULL;
	m_lostable_geometry         = NULL;
	m_lostable_optical          = NULL;
	m_optprop                   = NULL;
	m_observerspecs             = NULL;
	m_specifications            = NULL;
	m_singlescatter             = false;
	m_ignorehighalt				= false;

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::~SKTRAN_ThreadManager		2010-3-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ThreadManager::~SKTRAN_ThreadManager	()
{
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::ThreadExecuteAction		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::ThreadExecuteAction( SKTRANSO_Quadrature_TLS_V21*	threadquadrature, size_t pointindex )
{
	return (this->*m_threadfunction)(threadquadrature, pointindex);
}

