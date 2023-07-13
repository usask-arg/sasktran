#include "../sasktranv21_internals.h"
#include <float.h>

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::SKTRAN_DiffusePointOptical_V21		2007-12-14*/
/** **/
/*---------------------------------------------------------------------------*/


SKTRAN_DiffusePointOptical_V21::SKTRAN_DiffusePointOptical_V21(const SKTRAN_DiffusePointGeometry_V21* geometry )
{
	m_geometry         = geometry;
	m_numlosIn         = 0;
	m_numoutgoingJ     = 0;

	m_scatterobject.SetGeometry( geometry->ScatteringMatrix() );
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::SKTRAN_DiffusePointOptical_V21		2007-12-14*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_DiffusePointOptical_V21::~SKTRAN_DiffusePointOptical_V21()
{
	ReleaseLos();
	ReleaseOutgoingJ();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::ReleaseLos		2007-12-14*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_DiffusePointOptical_V21::ReleaseLos()
{
	m_losIn.clear();
	m_incomingradiance.clear();
	m_numlosIn         = 0;
//	m_losIn            = NULL;
//	m_incomingradiance = NULL;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::ReleaseRadiance		2007-12-14*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_DiffusePointOptical_V21::ReleaseOutgoingJ()
{
	m_outgoingJ.clear(); //if (m_outgoingJ           != NULL) delete [] m_outgoingJ;
	m_numoutgoingJ         = 0;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::AllocateLos		2007-12-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointOptical_V21::AllocateLos()
{
	bool	ok;
	size_t	newlos;

	ok = (m_geometry != NULL);														// do we have a geometry defined
	if (!ok)																		// nope
	{																				// so
		ReleaseLos();																// Release the lines of Sight
		ok = true;																	// and we are done
	}																				// otherwise
	else																			// we do
	{																				// an DiffusePointGeometry
		newlos = m_geometry->NumLOSIn();											// So get the number of lines of sight in the geometry
		ok = (newlos == m_numlosIn);												// See if it equals the number of lines of sight we already have
		if (!ok)																	// If it does not
		{																			// then
			ReleaseLos();															// release the current resources
			ok = (newlos == 0);														// if we need 0 elements
			if (!ok)																// then we are done otherwise lets do the allocation
			{																		// so
				m_losIn.resize(newlos); //           = new SKTRAN_RayInternalDiffuseOptical_V2[newlos];	// Allocate the incoming rays
				m_incomingradiance.resize(newlos);				// Allocate the incoming radiance
				ok = (m_losIn.size() == newlos) && (m_incomingradiance.size() == newlos);				// make sure it worked
				if (!ok)															// if it did not
				{																	// then log a message
					nxLog::Record(NXLOG_WARNING, "SKTRAN_DiffusePointOptical_V21::AllocateLos, Error allocating %Iu line of sight optical elements", (size_t) newlos);
					ReleaseLos();													// and release any memory allcoated
				}																	// and that is that
				else																// otherwise
				{																	// allocation successful
					m_numlosIn = newlos;											// setup the new number of lines of sight
				}																	// and we are done
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::AllocateOutgoingJ		2008-1-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointOptical_V21::AllocateOutgoingJ( const SKTRAN_UnitSphere_V2* unitsphere)
{
	size_t								numpts;
	bool								ok;

	NXASSERT(( m_geometry != NULL ));
	NXASSERT((unitsphere != NULL));
	numpts = unitsphere->NumUnitVectors();
	ok     = (numpts == m_numoutgoingJ);
	if (!ok)
	{
		ReleaseOutgoingJ();
		m_outgoingJ.resize(numpts);
		ok = (m_outgoingJ.size() == numpts);
		if (ok)
		{
			m_numoutgoingJ = numpts;
		}
		else
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_DiffusePointOptical_V21::AllocateOutgoingRadiance, Error allocating %Iu outgoing radiance points", (size_t)numpts);
			ReleaseOutgoingJ();
		}
	}
	ClearOutgoingJ();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::ClearOutgoingRadiances		2008-1-7*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_DiffusePointOptical_V21::ClearOutgoingJ()
{
	size_t	i;

	for (i = 0; i < m_numoutgoingJ; ++i)
	{
		m_outgoingJ[i]       = 0;					// Clear the diffuse scattered terms
	}
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::AttachGeometry		2007-12-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointOptical_V21::AttachToGeometry( )
{
	bool	ok;
	size_t	idx;
	bool	ok1;
//	const	SKTRANSO_RayStorage_InternalJindex*	ray;

	NXASSERT(( m_geometry != NULL));
	ok         =       m_scatterobject.AttachToGeometry();
	ok         = ok && AllocateLos();																			// and configure the internal optical line of sight array
	ok         = ok && AllocateOutgoingJ( m_geometry->ScatteringMatrix()->OutboundUnitSphere() );		// Allocate the outgoing radiance tables

	if (ok)
	{
		for (idx=0; idx < m_numlosIn; ++idx )
		{
//			ray = m_geometry->LOSAt( idx ).Indices();
			ok1 = m_losIn[idx].AttachToGeometry( &m_geometry->LOSAt( idx ) );
			ok  = ok && ok1;
		}
	}

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_DiffusePointOptical_V21::AttachToGeometry, Error attaching optical point to geometry, probably error allocating LOS or outgoing radiances");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::ConfigureOptical		2007-12-14*/
/** Ths is the code that **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointOptical_V21::ConfigureOptical( const SKTRAN_TableOpticalProperties_V21*			optprop,
														     SKTRANSO_Quadrature_TLS_V21*				quadrature)
{
	bool									ok = true;
	bool									ok1;
	size_t									idx;
	size_t									zenidx;
	size_t									aziidx;
	size_t									numzen;
	size_t									numazi;
	size_t									numcells;
	bool									ishomogenous;


	ok = m_scatterobject.ConfigureOptical	( optprop );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_DiffusePointOptical_V21::ConfigureOptical, There was an error configuring the scatterin matrix object");
	}
	else
	{
		ishomogenous = optprop->IsOptionTrue( SKTRAN_TableOpticalProperties_V21::OPTION_IS_HORIZONTALLY_UNIFORM ) ;		// Are the optical properties height profiles the same at all points

	//	NXASSERT(( (m_scatterobject.NumIncoming() == m_numlosIn) && (m_scatterobject.NumOutgoing() == m_numoutgoingJ) ));

		numzen  = m_geometry->IncomingUnitVectors()->NumZenith();

		for( zenidx = 0; zenidx < numzen; ++zenidx )																// Slowly scan over zenith angles
		{																														// For each zenith angle
			ok =       m_geometry->IncomingUnitVectors()->GetZenithVertexIndexAndNumVertexAzimuth( zenidx, &idx, &numazi );			// Get the starting vertex index for this zenith index as well as the number of azimuth points at this zenith angle
			ok = ok && m_losIn[idx].ConfigureOptical( quadrature, false, false,false);								// Configure the optical, diffuse points are never single scatter, reset the transmission and cell factor caches
			numcells = m_losIn[idx].GeometryRay()->Storage()->NumCells();
			idx++;
			if (!ok) nxLog::Record(NXLOG_WARNING,"SKTRAN_DiffusePointOptical_V21::ConfigureOptical, Error configuring zenith index %d", (int)zenidx);
			for( aziidx = 1; aziidx < numazi; ++aziidx )															// for each azimuth
			{																										// we can use homoegenous atmospheric symmetry to optimize the configuration
				ok1 = (!ishomogenous || (numcells == m_losIn[idx].GeometryRay()->Storage()->NumCells()));						// Make sure homogenous cells have the same number of elements
				NXASSERT(ok1);																						// and assert in debug mode
				if (!ok1) nxLog::Record(NXLOG_WARNING,"SKTRAN_DiffusePointOptical_V21::ConfigureOptical, ishomogenous = %d, numcells = %d, %d, azimuth index %d at zenith index %d at altitude = %12.5f", (int)ishomogenous, (int)numcells, (int)m_losIn[idx].GeometryRay()->Storage()->NumCells(), (int) aziidx, (int)zenidx, (double)Geometry()->Location().Altitude());
				ok1 = ok1 && m_losIn[idx++].ConfigureOptical( quadrature, false, ishomogenous, ishomogenous);				// as we can reuse the extiction calculation paths along rays at the same zenith angle
				if (ok && !ok1) nxLog::Record(NXLOG_WARNING,"SKTRAN_DiffusePointOptical_V21::ConfigureOptical, Error configuring azimuth index %d at zenith index %d", (int) aziidx, (int)zenidx);
				ok  = ok && ok1;
			}
		}
		NXASSERT(( idx == m_numlosIn));
	}
	ClearOutgoingJ();						// Clear the outbound radiances
//	ConfigureOutboundDiffuseForwardScatter();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_DiffusePointOptical_V21::ConfigureOptical, There was an error configuring this diffuse optical point. Thats not good.");
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::InitializeFirstOrderIncomingRadiances		2008-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointOptical_V21::InitializeFirstOrderIncomingRadiances()
{
	size_t	i;

	for (i = 0; i < m_numlosIn; ++i )
	{
		m_incomingradiance[i] = m_losIn[i].GetFirstOrderIncomingRadiance();
		NXASSERT(( NXFINITE( m_incomingradiance[i] ) ));
	}
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21::CalculateIncomingRadiation		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointOptical_V21::CalculateIncomingRadiances( bool ignorehighaltdiffuse )
{
	bool	ok = true;
	bool	ok1;
	size_t	i;

	ok = (ignorehighaltdiffuse && m_geometry->IsHighAltitude() );
	if (!ok)
	{
		ok = true;
		for (i = 0; i < m_numlosIn; ++i )
		{
			ok1 = m_losIn[i].CalculateTotalRadianceAtOrigin( &m_incomingradiance.at(i));
			ok  = ok && ok1;
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21::ScatterIncomingRadiance		2008-2-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointOptical_V21::ScatterIncomingRadiance( bool ignorehighaltitude )
{
	bool	ok;

	ok = (ignorehighaltitude && m_geometry->IsHighAltitude() );
	if (!ok) 
	{
		ok = m_scatterobject.ScatterIncomingRays( this );
	}
	return ok;
}

