#include "../sasktranv21_internals.h"



/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayLOSGeometry_V21::SKTRANSO_RayLOSGeometry_V21		2010-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_RayLOSGeometry_V21::SKTRANSO_RayLOSGeometry_V21()
{
	m_singlescattertable = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayLOSGeometry_V21::~SKTRANSO_RayLOSGeometry_V21		2010-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_RayLOSGeometry_V21::~SKTRANSO_RayLOSGeometry_V21()
{
	if (m_singlescattertable != NULL) m_singlescattertable->Release();
	m_singlescattertable = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayLOSGeometry_V21::SetINternalSolarTransmissionTable		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_RayLOSGeometry_V21::SetInternalSolarTransmissionTable	(SKTRANSO_TableRayLOS* table )
{
	if (table != NULL) table->AddRef();
	if (m_singlescattertable != NULL) m_singlescattertable->Release();
	m_singlescattertable = table;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayLOSGeometry_V21::ConfigureInternalSolarTransmissionTable		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_RayLOSGeometry_V21::ConfigureInternalSolarTransmissionTableGeometry(SKTRANSO_Quadrature_TLS_V21* quadrature)
{
	bool	ok;

	ok = (m_singlescattertable == NULL);
	if (!ok)
	{
		ok = m_singlescattertable->ConfigureGeometry( quadrature, this );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRANSO_RayLOSGeometry_V21::ConfigureInternalSolarTransmissionTableGeometry, Error configuring internal single scatter table");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayLOSGeometry_V21::ConfigureInternalSolarTransmissionTableOptical		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_RayLOSGeometry_V21::ConfigureInternalSolarTransmissionTableOptical(SKTRANSO_Quadrature_TLS_V21* quadrature)
{
	bool	ok;

	ok = (m_singlescattertable == NULL);
	if (!ok)
	{
		ok = m_singlescattertable->ConfigureOptical( quadrature );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRANSO_RayLOSGeometry_V21::ConfigureInternalSolarTransmissionTableOptical, Error configuring internal single scatter table");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayLOSGeometry_V21::InternalSolarTransmissionTable		 2014- 11- 26*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRANSO_JindexTableBase*	 SKTRANSO_RayLOSGeometry_V21::InternalSolarTransmissionTable() const
{
	return m_singlescattertable;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayLOSOptical_V21::SKTRAN_RayLOSOptical_V21		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayLOSOptical_V21::SKTRAN_RayLOSOptical_V21()
{
	m_losgeometryray = NULL;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayLOSOptical_V21::~SKTRAN_RayLOSOptical_V21		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayLOSOptical_V21::~SKTRAN_RayLOSOptical_V21()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayLOSOptical_V21::AttachToGeometry		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayLOSOptical_V21::AttachToGeometry( SKTRANSO_RayLOSGeometry_V21* los)
{
	bool	ok;

	m_losgeometryray = los;
	ok = SKTRAN_RayInternalDiffuseOptical_V2::AttachToGeometry(los);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayLOSOptical_V21::ConfigureOptical		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayLOSOptical_V21::ConfigureOptical( SKTRANSO_Quadrature_TLS_V21* quadrature, bool singlescatter, bool usecachedtransmission, bool usecachedcellfactors )
{
	bool	ok;

	ok = (m_losgeometryray != NULL);
	if (ok)
	{
		ok =       m_losgeometryray->ConfigureInternalSolarTransmissionTableOptical( quadrature );
		ok = ok && SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical( quadrature, singlescatter, usecachedtransmission, usecachedcellfactors );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayLOSOptical_V21::ConfigureOptical, There was an error configuring the line of sight");
	}
	return ok;


};



















/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS::SKTRANSO_TableRayLOS		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_TableRayLOS::SKTRANSO_TableRayLOS( std::weak_ptr< const SKTRAN_RayFactory_Base> rayfactory)
{
	m_losrayfactory         = rayfactory;
	m_ray                   = NULL;
	m_distancefromobserver.SetStatic();

}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS::~SKTRANSO_TableRayLOS		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_TableRayLOS::~SKTRANSO_TableRayLOS()
{
	ReleaseStorage();
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS::InterpolateTable		2010-6-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableRayLOS::InterpolateTable( const HELIODETIC_POINT&				location,
												  const HELIODETIC_UNITVECTOR&	look,
												  SKTRANSO_JIndex*			vertexdescriptortable,
												  size_t						maxpts,
												  size_t*						npts,
												  double						userweight
												) const
{
	double				s;
	SKTRAN_GridIndex	indices[2];
	double				weight[2];
	size_t				numpts;
	size_t				gidx;
	size_t				idx;
	bool				ok;
	HELIODETIC_VECTOR	o;
	HELIODETIC_VECTOR	p;


	o = m_ray->Ray()->GetObserver();
	p = location.Vector();
	s = (p-o).Magnitude();
	
	/* 2014-07-29 djz828, changed interpolation mode from OUTOFBOUND_ZERO to OUTOFBOUND_TRUNCATE 
	 * OUTOFBOUND_ZERO was causing problems when the line of sight ray was looking at the ground,
	 * occasionally s would be larger than the last element in m_distancefromobserver (from rounding errors), and 
	 * InterpolateTable would return 0.0.
	 */
	numpts = m_distancefromobserver.FindingBoundingIndices( s, SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, &(indices[0]), &(weight[0]), 2);
	gidx = 0;
	ok   = true;
	for (idx = 0; idx < numpts; idx++)
	{
		if (weight[idx] != 0.0)
		{
			NXASSERT(( weight[idx] > 0.0 ));
			ok = ok && (gidx < maxpts);
			if (ok)
			{
				vertexdescriptortable[gidx].ConfigureLOSSolarTransmissionIndex( indices[idx], weight[idx]*userweight );
				++gidx;
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableRayLOS_Legacy::InterpolateTable, There was an error interpolating the Line of Sight Solar transmission table. Thats not good");
		gidx = 0;
	}
	*npts = gidx;
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS::ConvertJIndexToRadiancePtr		2010-6-18*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_StokesScalar* SKTRANSO_TableRayLOS::ConvertJIndexToRadiancePtr( const SKTRANSO_JIndex*		entry,
																				 ENUM_SKTRAN_JSOURCE			jsource
																		 ) const
{
	size_t	idx;

	idx = entry->HeightIndex();
	return (&m_incomingsolarradiance.at(idx));
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS::ReleaseStorage		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_TableRayLOS::ReleaseStorage()
{
	m_distancefromobserver.erase();
	m_raystosun.clear();
	m_incomingsolarradiance.clear();
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS::AllocateStorage		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableRayLOS::AllocateStorage( size_t numpoints)
{
	bool	ok;

	ok = (numpoints == NumRaysToSun());
	if (!ok)
	{
		ReleaseStorage();
		ok = (numpoints == 0);
		if (!ok)
		{
			m_raystosun.resize( numpoints ); 
			m_incomingsolarradiance.resize( numpoints);
			m_distancefromobserver.AllocateGridArray(numpoints);
			ok =    ( m_raystosun.size() == numpoints) && ( m_incomingsolarradiance.size() == numpoints) && (m_distancefromobserver.size() == numpoints);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRANSO_TableRayLOS, Error allocating internal memory for single scatter table");
				ReleaseStorage();
			}
		}
	}
	return ok;
		
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableRayLOS_Legacy::ConfigureEntry		2010-6-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableRayLOS::ConfigureEntry( size_t index, double distancefromorigin, SKTRANSO_Quadrature_TLS_V21* quadrature )
{
	HELIODETIC_VECTOR								point;
	HELIODETIC_UNITVECTOR							look;
	bool											ok;
	HELIODETIC_VECTOR								origin;;
	std::unique_ptr< SKTRAN_RayOptical_Base>		raytosun;
	std::shared_ptr< const SKTRAN_RayFactory_Base>	factory;


//	static bool firsttime = true;
//	if (firsttime)
//	{
//		nxLog::Record(NXLOG_INFO, "SKTRANSO_TableRayLOS::ConfigureEntry, **** TODO **** Tweak up code so it works for curved and straight rays ");
//		firsttime = false;
//	}
	
	factory = LOSRayFactory().lock();

	NXASSERT((m_ray != NULL));
	look.SetCoords(0.0, 0.0, 1.0);

	point.SetCoords( m_ray->Ray()->LookVector(), distancefromorigin);
	point += m_ray->Ray()->GetObserver();

	ok =       factory->CreateRayObject( &raytosun );
	ok = ok && raytosun->MoveObserver( point, look);
	ok = ok && raytosun->TraceRay_NewMethod();

	m_raystosun.at(index).AssignRay( std::move( raytosun ) );
	m_distancefromobserver.AtVar(index) = distancefromorigin;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS::ConfigureRayToSunLocations		2010-6-18*/
/** A helper function for derived classes. The derived class passes in a list
 *	of locations along the ray where it wants single scatter solar transmission
 *	entries to be placed.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableRayLOS::ConfigureRayToSunLocations( std::vector< SKTRAN_Distance>& shelllocations, SKTRANSO_Quadrature_TLS_V21* quadrature,	 SKTRANSO_RayLOSGeometry_V21* losray )
{
	bool								ok;
	bool								ok1;
	size_t								idx;
	size_t								numpoints;

	m_ray = losray;
	numpoints = shelllocations.size();;														
	ok       = AllocateStorage( numpoints );
	for (idx = 0; idx < numpoints; idx++)
	{
		ok1 = ConfigureEntry( idx, shelllocations[idx], quadrature);
		ok  = ok && ok1;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS::ConfigureOptical		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableRayLOS::ConfigureOptical( SKTRANSO_Quadrature_TLS_V21 *	quadrature )
{

	bool									ok = true;
	bool									ok1;
	const SKTRAN_RayStorage_Base*			raygeometry;
	SKTRANSO_RayInternalOptical				rayoptical;
	size_t									idx;

	NXTRACE_ONCEONLY(firsttime, (" SKTRANSO_TableRayLOS::ConfigureOptical, We could configure rayoptical so it does not continuosly allocate and deallocate geometry settings\n"));

	for (idx = 0; idx < NumRaysToSun(); ++idx)
	{
		raygeometry = m_raystosun.at(idx).Storage();
		if (raygeometry->GroundIsHit())													// IF this ray hits the ground
		{																				// then
			m_incomingsolarradiance[idx] = 0.0;											// we get nothing from the sun
		}																				// and that is that
		else																			// otherwise
		{															
			ok1 =        rayoptical.AttachToGeometry( &m_raystosun.at(idx) );					// we need to attch 
			ok1 = ok1 && rayoptical.ConfigureOptical( quadrature, true, false );		// configure
			m_incomingsolarradiance[idx] = rayoptical.TotalTransmissionAlongRay();		// and get the total transmission
			ok  = ok && ok1;
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_TableRayLOS::ConfigureOptical, Error configuring the optical component of the Solar Transmission table");
	}
	return ok;
}
