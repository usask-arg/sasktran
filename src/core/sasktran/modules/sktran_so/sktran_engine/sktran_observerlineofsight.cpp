#include "../sasktranv21_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSight_V21::SKTRAN_TableLinesOfSight_V21		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableLinesOfSight_V21::SKTRAN_TableLinesOfSight_V21()
{
	m_opticaltable         = new SKTRAN_TableLinesOfSightOptical_V21(this);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSight_V21::~SKTRAN_TableLinesOfSight_V21		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableLinesOfSight_V21::~SKTRAN_TableLinesOfSight_V21()
{
	ReleaseResources();
	delete m_opticaltable;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSight_V21::ReleaseResources		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableLinesOfSight_V21::ReleaseResources()
{
	m_raygeometry.clear();
	m_opticaltable->ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSight_V21::AllocateRays		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableLinesOfSight_V21::AllocateRays( size_t n, const SKTRAN_TableRayLOSFactory* singlescatterfactory )
{
	bool					ok;
	SKTRANSO_TableRayLOS*		lostable;
	size_t					idx;
	bool					ok1;

	ReleaseResources();
	ok = (n == 0);
	if (!ok)
	{
		m_raygeometry.resize(n);
		for (size_t i = 0; i < n; i++)
		{
			std::unique_ptr< SKTRAN_RayOptical_Base	>				ray;
			m_rayfactory->CreateRayObject( &ray );
			m_raygeometry.at(i).AssignRay( std::move(ray) );
		}
		ok = m_raygeometry.size() == n;
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_TableLinesOfSight_V21::AllocateRays, Error allocating space for %Iu observer line of sight rays", (size_t)n);
		}
		else
		{
			for (idx =0; idx < NumRays(); ++idx)
			{
				if (singlescatterfactory != NULL)
				{
					ok1 =     singlescatterfactory->CreateInternalSingleScatterTable( &lostable, m_solarrayfactory );
				}
				else
				{
					ok1 = true;
					lostable = NULL;
				}
				ok1 = ok1 && m_raygeometry.at(idx).SetInternalSolarTransmissionTable( lostable );
				if (lostable != NULL) lostable->Release();
				ok = ok && ok1;
			}
		}
	}
	if (!ok) ReleaseResources();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSight_V21::TranslateLimbViewingGeometryToOsculatingSphere		2014-4-8*/
/** Converts an observer**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableLinesOfSight_V21::TranslateLimbViewingGeometryToOsculatingSphere( SKTRAN_LineOfSightEntry_V2* entry, const SKTRAN_CoordinateTransform_V2* coordinates)
{
	double								latitude;
	double								longitude;
	double								height;
	nxVector							offsetobserver;
	nxVector							offsettangentpoint;
	nxVector							look;
	nxGeodetic							earthgeoid 	= coordinates->TrueGeoid();
	nxGeodetic							oscgeoid    = coordinates->OsculatingGeoid();


	earthgeoid.FromTangentPointLocation( entry->Observer(), entry->Look() );										// Get the tangent point locations in the True Earth system
	latitude  = earthgeoid.GeodeticLatitude();																		// Get the latitude, longitude and height
	longitude = earthgeoid.GeodeticLongitude();																		// of the tangent point
	height    = earthgeoid.Height();
	NXTRACE(("Original              line of sight tangent point (lat = %12.7f, lng = %12.7f, height= %15.7f)\n", (double)latitude, (double)longitude, (double)height));

	offsetobserver = coordinates->TranslateGeoidToOsculatingSphere( entry->Observer() );							// Get the observers geocentric location in the osculating sphere system

	#if defined(NXDEBUG)
		oscgeoid.FromTangentPointLocation(offsetobserver, entry->Look() );
		double alatitude  = oscgeoid.GeodeticLatitude();
		double alongitude = oscgeoid.GeodeticLongitude();
		double aheight    = oscgeoid.Height();
	NXTRACE(("Old Osculating sphere line of sight tangent point (lat = %12.7f, lng = %12.7f, height= %15.7f)\n", (double)alatitude, (double)alongitude, (double)aheight));

	#endif
	oscgeoid.FromGeodetic( latitude, longitude, height );															// Get the tangent point location
	offsettangentpoint = oscgeoid.Location();																		// in the osculating sphere location
	look = (offsettangentpoint-offsetobserver).UnitVector();
	entry->Configure( offsetobserver, look, entry->Mjd() );

	#if defined(NXDEBUG)
		oscgeoid.FromTangentPointLocation( entry->Observer(), entry->Look());
		latitude  = oscgeoid.GeodeticLatitude();
		longitude = oscgeoid.GeodeticLongitude();
		height    = oscgeoid.Height();
		NXTRACE(("New Osculating sphere line of sight tangent point (lat = %12.7f, lng = %12.7f, height= %15.7f)\n", (double)latitude, (double)longitude, (double)height));
	#endif
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSight_V21::SetLinesOfSight		2010-4-13*/
/** Copies the users lines of sight and stores them internally.
 *
 *	\par DEPRECATED
 *	The lines of sight are adjusted so they have the same tangent altitude in the osculating
 *	sphere as they do in the real Earth geoid. As of 2014-04-08 we stopped this as it
 *	required changes in the ray direction which would change the solar scattering angle.
 *	The technique may be re-instated at a later date but for now we removed the correction
 *	as it may do more harm than good.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableLinesOfSight_V21::SetLinesOfSight( const SKTRAN_LineOfSightArray_V21& observerlinesofsight, const SKTRAN_CoordinateTransform_V2* coordinates,  SKTRAN_RayTracingRegionManager*	raytracingregionmanager )
{
//	size_t								idx;
	bool								ok;
//	bool								ok1;
//	SKTRAN_LineOfSightEntry_V2*			entry;

	ReleaseResources();
	
	ok = m_observerlinesofsight.DeepCopy( observerlinesofsight );			// Copy over the rays.
//	nxLog::Record(NXLOG_INFO,"SKTRAN_TableLinesOfSight_V21::SetLinesOfSight, We should determine if we should adjust rays as we go from true Earth to Osculating sphere");
/*	for (idx = 0; idx < m_observerlinesofsight.NumRays(); idx++ )
	{
		ok1 = m_observerlinesofsight.GetRayVar( idx, &entry );
		if (ok1 )
		{
			earthgeoid.FromTangentPointLocation( entry->Observer(), entry->Look() );
			latitude  = earthgeoid.GeodeticLatitude();
			longitude = earthgeoid.GeodeticLongitude();
			height    = earthgeoid.Height();
			NXTRACE(("Original              line of sight tangent point (lat = %12.7f, lng = %12.7f, height= %15.7f)\n", (double)latitude, (double)longitude, (double)height));

			offsetobserver = coordinates->TranslateGeoidToOsculatingSphere( entry->Observer() );

	#if defined(NXDEBUG)
			oscgeoid.FromTangentPointLocation(offsetobserver, entry->Look() );
			double alatitude  = oscgeoid.GeodeticLatitude();
			double alongitude = oscgeoid.GeodeticLongitude();
			double aheight    = oscgeoid.Height();
			NXTRACE(("Old Osculating sphere line of sight tangent point (lat = %12.7f, lng = %12.7f, height= %15.7f)\n", (double)alatitude, (double)alongitude, (double)aheight));

	#endif
			oscgeoid.FromGeodetic( latitude, longitude, height );
			offsettangentpoint = oscgeoid.Location();
			look = (offsettangentpoint-offsetobserver).UnitVector();
			entry->Configure( offsetobserver, look, entry->Mjd() );

	#if defined(NXDEBUG)
			oscgeoid.FromTangentPointLocation( entry->Observer(), entry->Look());
			latitude  = oscgeoid.GeodeticLatitude();
			longitude = oscgeoid.GeodeticLongitude();
			height    = oscgeoid.Height();
			NXTRACE(("New Osculating sphere line of sight tangent point (lat = %12.7f, lng = %12.7f, height= %15.7f)\n", (double)latitude, (double)longitude, (double)height));
	#endif
		}
		ok = ok && ok1;

	}
*/
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableLinesOfSight_V21::SetLinesOfSight, There were errors converting the lines of sight to osculating sphere system. Thats not good.");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSight_V21::ConfigureGeometry		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableLinesOfSight_V21::ConfigureGeometry( bool										 singlescatter,
													  const SKTRAN_SpecsInternal_V21*			 modelspecs,
													  SKTRAN_ThreadManager*						 threadmanager)
{
	bool										 ok;
	size_t										 n;

	n            =       m_observerlinesofsight.NumRays( );
	m_rayfactory = modelspecs->RayFactoryLOS();
	m_solarrayfactory = modelspecs->RayFactoryTransmissionOnly();
	ok = (m_rayfactory.get() != nullptr);
	ok = ok && AllocateRays( n, modelspecs->DiffuseSpecs()->LOSSingleScatterTableFactory() );
	ok = ok && m_opticaltable->AttachToGeometry();
	ok = ok && threadmanager->LinesOfSightTable_ConfigureGeometryStage2( this, &m_observerlinesofsight, modelspecs, singlescatter );
	if (!ok)
	{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_TableLinesOfSight_V21::ConfigureGeometry, There were errors creating the observer line of sight rays and linking the Jindex tables. This is going to cause problems");
			ReleaseResources();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableLinesOfSight_V21::ConfigureOptical		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableLinesOfSight_V21::ConfigureOptical( bool singlescatter,  const SKTRAN_TableOpticalProperties_V21*  optprop,  SKTRAN_ThreadManager* threadmanager )
{
	return m_opticaltable->ConfigureOptical( singlescatter, optprop, threadmanager );
}


