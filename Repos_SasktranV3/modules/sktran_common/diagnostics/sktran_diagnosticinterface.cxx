#include "../sktran_common.h"
#include <nxbase_threads.h>

static std::mutex							g_diagnosticmutex;
//	static nxMutex							g_diagnosticmutex;

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiagnosticInterface::SKTRAN_DiagnosticInterface		2009-1-15*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_DiagnosticInterface::SKTRAN_DiagnosticInterface()
{
	m_engineinterface = NULL;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiagnosticInterface::SKTRAN_DiagnosticInterface		2009-1-15*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_DiagnosticInterface::~SKTRAN_DiagnosticInterface()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiagnosticInterface::Initialize		2009-1-15*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_DiagnosticInterface::Initialize(SKTRAN_Engine_Base* engine)
{
	m_engineinterface = engine;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiagnosticInterface::AddDiagnosticFunction		2009-1-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiagnosticInterface::AddDiagnosticFunction(SKTRAN_DiagnosticFunction func)
{
	m_diagnosticfunctions.push_back( func );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiagnosticInterface::ClearDiagnosticFunctions		2009-1-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiagnosticInterface::ClearDiagnosticFunctions()
{
	m_diagnosticfunctions.clear();
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiagnosticInterface::SetWavelength		2010-4-13*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_DiagnosticInterface::SetWavelength( double wavelenm )
{
	m_wavelennm = wavelenm;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiagnosticInterface::Diagnose		2009-1-12*/
/** Execute user-supplied diagnostic functions. This is invoked from the optical engin **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiagnosticInterface::Diagnose( ENUM_SKTRAN_DIAGNOSTIC_STAGE stage, size_t ordersofscatter , bool processingok )
{
	std::list<SKTRAN_DiagnosticFunction>::iterator			iter;
	SKTRAN_DiagnosticFunction								func;
	bool													ok;
	bool													ok1;
	
	ok = (m_diagnosticfunctions.size() == 0);
	if (!ok)
	{
		std::unique_lock<std::mutex>	lock( g_diagnosticmutex );			// lock this thread while diagnosing the thread, avoids multithread I/O issues
//		nxSingleLock						lock(&m_diagnosticmutex, TRUE);

		ok = (m_engineinterface != NULL);
		NXASSERT(( m_engineinterface != NULL ));		// Must call initialize  in SKTRAN_EningInterface constructor
		for (iter = m_diagnosticfunctions.begin(); !(iter == m_diagnosticfunctions.end()); ++iter)
		{
			func = *iter;
			ok1 = (*func)(  stage, m_wavelennm, ordersofscatter, processingok, m_engineinterface);
			ok = ok && ok1;
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_EngineInterface_V2::ExecuteDiagnosticFunctions, The diagnostic function returned an error. That might be a problem");
			}
		}
	}
	return ok;
}


