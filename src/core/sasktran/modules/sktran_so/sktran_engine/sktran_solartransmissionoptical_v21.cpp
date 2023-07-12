#include "../sasktranv21_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionOptical_V21::SKTRAN_TableSolarTransmissionOptical_V21 2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableSolarTransmissionOptical_V21::SKTRAN_TableSolarTransmissionOptical_V21( SKTRANSO_TableSolarTransmission*	solargeom )
{
	m_geometrytable = solargeom;
	m_transmission  = NULL;
	m_numpoints    = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionOptical_V21::~SKTRAN_TableSolarTransmissionOptical_V2_Uniform2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableSolarTransmissionOptical_V21::~SKTRAN_TableSolarTransmissionOptical_V21()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionOptical_V2_Uniform_Uniform::ReleaseResources		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableSolarTransmissionOptical_V21::ReleaseResources()
{
	if (m_transmission != NULL ) delete [] m_transmission;
	m_transmission = NULL;
	m_numpoints    = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionOptical_V2_Uniform_Uniform::AllocateMemory		2009-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableSolarTransmissionOptical_V21::AllocateMemory( size_t numpoints )
{
	bool		ok;

	ok = (numpoints == m_numpoints);
	if (!ok)
	{
		ReleaseResources();
		ok = (numpoints == 0);
		if (!ok)
		{
			m_transmission = new SKTRAN_StokesScalar [numpoints];
			ok = (m_transmission != NULL);
			if (ok)
			{
				m_numpoints = numpoints;
			}
			if (!ok)
			{
				nxLog::Record( NXLOG_WARNING, "SKTRAN_TableSolarTransmissionOptical_V2_Uniform_Uniform::AllocateMemory, Error allocating memory, this is a problem");
				ReleaseResources();
			}
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionOptical_V21::Configure		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableSolarTransmissionOptical_V21::AttachToGeometry()
{
	bool ok;

	ok = AllocateMemory( m_geometrytable->NumPoints() );
	return ok;

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionOptical_V2_Uniform_Uniform::ConfigureOptical		2008-1-4*/
/** Initializes the optical properties of the solar transmission table.
 *	The optical depth of each point in the table from the sun is calculated.
 *	At this point the table is completed.
  **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableSolarTransmissionOptical_V21::ConfigureOptical_Stage2( size_t rayidx, SKTRANSO_Quadrature_TLS_V21* threadquadrature )
{
	bool												ok;
	bool												ok1;
	const SKTRANSO_RayInternalGeometry*					raygeometry;
	SKTRANSO_RayInternalOptical							rayoptical;

	ok = m_geometrytable->RayAtPoint( rayidx, &raygeometry );
	if (ok)
	{
		if (raygeometry->Storage()->GroundIsHit())													// IF this ray hits the ground
		{																				// then
				m_transmission[rayidx] = 0.0;											// we get nothing from the sun
		}																				// and that is that
		else																			// otherwise
		{															
			ok1 =  rayoptical.AttachToGeometry( raygeometry );							// we need to attch 
			ok1 = ok1 && rayoptical.ConfigureOptical( threadquadrature, true, false );	// configure
			m_transmission[rayidx] = rayoptical.TotalTransmissionAlongRay();			// and get the total transmission
			ok  = ok && ok1;
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableSolarTransmissionOptical_V2_Uniform_Uniform::ConfigureOptical, Error configuring the optical component of the Solar Transmission table");
		ReleaseResources();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionOptical_V21::RayAt		2008-2-9*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_StokesScalar* SKTRAN_TableSolarTransmissionOptical_V21::RayAt( size_t position )const
{
	NXASSERT(( position < m_numpoints ));
	return m_transmission + position;
}
