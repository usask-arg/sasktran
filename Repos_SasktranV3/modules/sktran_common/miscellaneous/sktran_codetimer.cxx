#include "../sktran_common.h"
#include <time.h>


/*-----------------------------------------------------------------------------
 *					SKTRAN_CodeTimer::SKTRAN_CodeTimer		2008-3-26*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_CodeTimer::SKTRAN_CodeTimer()
{
	Start();
	m_t2 = 0;
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_CodeTimer::Start		2008-3-26*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_CodeTimer::Start( )
{
	m_t1 = (size_t)clock();
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_CodeTimer::Stop		2008-3-26*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_CodeTimer::Stop( )
{
	m_t2 = (size_t)clock();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_CodeTimer::Report		2008-3-26*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_CodeTimer::Report( const char* message )
{
	double	dt;
	double duration;

	if (m_t2 == 0) Stop();
	dt       = (double)(m_t2 - m_t1);
	duration = dt/CLOCKS_PER_SEC;

	nxLog::Verbose(NXLOG_INFO, "Code Timer:: %s = %5.3f seconds\n", (const char*)message, (double)duration);
	NXTRACE(("Code Timer:: %s = %5.3f seconds\n", (const char*)message, (double)duration));
}
