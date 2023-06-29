#include "../sasktranv21_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSightOptical_V21::SKTRAN_TableLinesOfSightOptical_V21		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableLinesOfSightOptical_V21::SKTRAN_TableLinesOfSightOptical_V21(SKTRAN_TableLinesOfSight_V21* geometry)
{
	m_numrays         = 0;
	m_rayoptical      = NULL;
	m_geometry        = geometry;
	m_numraysreserved = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSightOptical_V21::~SKTRAN_TableLinesOfSightOptical_V21		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableLinesOfSightOptical_V21::~SKTRAN_TableLinesOfSightOptical_V21()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSightOptical_V21::ReleaseResources		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableLinesOfSightOptical_V21::ReleaseResources()
{
	if (m_rayoptical != NULL) delete [] m_rayoptical;
	m_rayoptical = NULL;
	m_numrays         = 0;
	m_numraysreserved = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSightOptical_V21::AllocateRays		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableLinesOfSightOptical_V21::AllocateRays( size_t n )
{
	bool	ok;
	size_t	numreserved;

	ok = (n <= m_numraysreserved);
	if (!ok)
	{
		ReleaseResources();
		numreserved = (n*3+1)/2;												// Add quite a bit of slop so we dont keep reallocating.
		m_rayoptical = new SKTRAN_RayLOSOptical_V21[numreserved];
		ok = (m_rayoptical != NULL );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_TableLinesOfSightOptical_V21::AllocateRays, Error allocating space for %Iu observer lines of sight rays. This is a problem", (size_t)n);
		}
		else
		{
			m_numrays = n;
			m_numraysreserved  = numreserved;
		}
	}
	if (!ok) ReleaseResources();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSightOptical_V21::AttachToGeometry		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableLinesOfSightOptical_V21::AttachToGeometry( )
{
	bool											ok;
	bool											ok1;
	size_t											idx;

	NXASSERT(( m_geometry != NULL ));
	ok         = AllocateRays( m_geometry->NumRays() );
	for (idx = 0; idx < m_numrays; idx++)
	{
		ok1 = m_rayoptical[idx].AttachToGeometry( m_geometry->RayAtVar(idx) );
		ok = ok && ok1;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableLinesOfSightOptical_V21::AttachToGeometry, There was an error attaching the optical observer lines of sight to the geometry");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSightOptical_V21::ConfigureOptical		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableLinesOfSightOptical_V21::ConfigureOptical( bool singlescatter, const SKTRAN_TableOpticalProperties_V21*  optprop, SKTRAN_ThreadManager* threadmanager )
{
	bool	ok;

	ok = threadmanager->LinesOfSightTable_ConfigureOpticalStage2( this, singlescatter );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableLinesOfSightOptical_V21::ConfigureOptical, There was an error configuring the optical observer lines of sight");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSightOptical_V21::CalculateObserverIncomingRadiance		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableLinesOfSightOptical_V21::CalculateObserverIncomingRadiance( SKTRAN_StokesScalar* buffer, size_t maxrays)
{
	size_t	i;
	bool	ok1;
	bool	ok;

	ok = (maxrays >= m_numrays);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableLinesOfSightOptical_V21::CalculateObserverIncomingRadiance, user must provide exact space for number of rays");
	}
	else
	{
		NXTRACE_ONCEONLY(firsttime,("SKTRAN_TableLinesOfSightOptical_V21::CalculateObserverIncomingRadiance, it might be possible to multi thread this code\n"));
		for (i = 0; i < m_numrays; ++i )
		{
			ok1 = m_rayoptical[i].CalculateTotalRadianceAtOrigin(buffer+i);
			ok = ok && ok1;
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_TableLinesOfSightOptical_V21::CalculateObserverIncomingRadiance, error calculating incoming radiance");
		}
	}
	return ok;
}



