#include "../sasktranv21_internals.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::SKTRAN_DiffusePointGeometry_V21		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_DiffusePointGeometry_V21::SKTRAN_DiffusePointGeometry_V21()
{
	m_diffuseheightindex    = (size_t)(-1);
	m_diffuseprofileindex   = (size_t)(-1);
//	m_losGeometryIn         = NULL;
//	m_numlosIn              = 0;
//	m_raytracingspecs       = NULL;
//	m_coords                = NULL;
	m_incomingunitvectors   = NULL;
//	m_zenithsIn             = NULL;
//	m_azimuthsIn            = NULL;
	m_ishighaltitudediffuse = false;
	m_opticalpoint          = NULL;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::CreateOpticalPoint		2012-9-14*/
/** We used to create the Optical point inside the constructor but this
 *	created problems when the points were put inside std::vector using vector.resize
 *	as the pointer was copied in the resizing.  The allocation is now done after we
 *	have allocated the std::vector.
 **/ 
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointGeometry_V21::CreateOpticalPoint()
{
	bool	ok;
	ok = (m_opticalpoint == NULL);
	if (ok)
	{
		m_opticalpoint  = new SKTRAN_DiffusePointOptical_V21(this);
		m_opticalpoint->AddRef();

	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_DiffusePointGeometry_V21::CreateOpticalPoint, OpticalPoint already defined, this function should only be called once per instance as part of the object construction");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::~SKTRAN_DiffusePointGeometry_V21		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_DiffusePointGeometry_V21::~SKTRAN_DiffusePointGeometry_V21()
{
	ReleaseResources();
	if (m_opticalpoint != NULL) m_opticalpoint->Release();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::ReleaseResources		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_DiffusePointGeometry_V21::ReleaseResources()
{
	m_losGeometryIn.clear();
//	if ( m_scatteringmatrix != NULL ) m_scatteringmatrix ->Release();
//	if ( m_azimuthsIn       != NULL ) m_azimuthsIn->Release();
//	if ( m_zenithsIn        != NULL ) m_zenithsIn->Release();
	if ( m_incomingunitvectors != NULL ) m_incomingunitvectors->Release();
//	if ( m_raytracingspecs  != NULL ) m_raytracingspecs->Release();
//	if ( m_coords           != NULL ) m_coords->Release();

	m_diffuseheightindex    = (size_t)(-1);
	m_diffuseprofileindex   = (size_t)(-1);
//	m_losGeometryIn         = NULL;
//	m_numlosIn              = 0;
//	m_raytracingspecs       = NULL;
//	m_zenithsIn             = NULL;
//	m_azimuthsIn            = NULL;
	m_incomingunitvectors   = NULL;
//	m_scatteringmatrix      = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::GeographicZenith		2008-3-20*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_UNITVECTOR SKTRAN_DiffusePointGeometry_V21::LocalZenith() const
{
	return  m_location.LocalZenith();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::AllocateLOSArray		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointGeometry_V21::AllocateLOSArray()
{
	bool	ok;
	size_t	numlosIn;

	NXASSERT(( ( m_incomingunitvectors != NULL)  ));		// Make sure the grids are defined
	//NXASSERT(( m_losGeometryIn == NULL ));									// Make sure the LOS array is undefined

	numlosIn   = m_incomingunitvectors->NumUnitVectors();					// Get the number of unit vectors 
	ok = (numlosIn > 0);													// Make sure we have something to do
	if (!ok)																// if not then log
	{																		// a warning
		nxLog::Record(NXLOG_WARNING, "SKTRAN_DiffusePointGeometry_V2::AllocateLOSArray, No lines of sight to allocate" );
	}																		// and that is that
	else																	// otherwise
	{																		// do the allocation
		m_losGeometryIn.clear();
		m_losGeometryIn.resize(numlosIn);
//		for (size_t i = 0; i < numlosIn; i++)
//		{
//			std::unique_ptr<SKTRAN_RayStorage_Base>		trajectorystorage;
//			m_losGeometryIn.push_back( std::move(std::unique_ptr< SKTRANSO_RayInternalGeometry> ( new SKTRANSO_RayInternalGeometry( std::move(trajectorystorage), nullptr) ) ));
//		}
		ok = (m_losGeometryIn.size() == numlosIn);
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING, "SKTRAN_DiffusePointGeometry_V21::AllocateLOSArray, Error allocating memory for %Iu line of sight rays", (size_t)numlosIn );
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::ConfigureLOSArray		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointGeometry_V21::ConfigureLOSArray( SKTRANSO_Quadrature_TLS_V21* quadraturespecs )
{
	bool										ok;
	bool										ok1;
	double										x,y,z;
	HELIODETIC_UNITVECTOR						xprime;
	HELIODETIC_UNITVECTOR						yprime;
	HELIODETIC_UNITVECTOR						zprime;
	size_t										numzen;
	size_t										numazi;
	size_t										izen;
	size_t										iazi;
//	double										zenith;
//	double										azimuth;
//	double										coszen;
//	double										sinzen;
///	double										cosazi;
//	double										sinazi;
	nxVector									v;
	size_t										idx;
	HELIODETIC_UNITVECTOR						look;
	HELIODETIC_UNITVECTOR						unitvectors[3];

	m_location.LocalUnitVectors( &unitvectors[0], 3);		// Get the profile's local x',y',z' unit vectors from the parent profile
	xprime = unitvectors[0];
	yprime = unitvectors[1];
	zprime = unitvectors[2];

	ok = AllocateLOSArray();
	if (ok)
	{
		NXTRACE_ONCEONLY(firsta,("**** 2008-12-23 ***** SKTRAN_DiffusePointGeometry_V2::ConfigureLOSArray, Potential bug if sun is directly parallel to ZPrime\n"));
		numzen  = m_incomingunitvectors->NumZenith();											// Get the number of zeniths in the line of sight table
		for( izen = 0; izen < numzen; ++izen )											// Slowly scan over zenith angles
		{																				// For each zenith angle
			ok1     = m_incomingunitvectors->GetZenithVertexIndexAndNumVertexAzimuth( izen, &idx, &numazi );

			for( iazi = 0; iazi < numazi; ++iazi )										// quickly scan over azimuth
			{																			// for each azimuth

				v = m_incomingunitvectors->UnitVectorAt( idx );
				x = v.X();		// sinzen*cosazi;										// NOte that azimuth is zero degrees at the due south position
				y = v.Y();      // sinzen*sinazi;										// and is 90 degrees or pi/2 at the due East direction (in the system with sun at pole)
				z = v.Z();		// coszen;
				look.SetCoords( x*xprime.X() + y*yprime.X() + z*zprime.X(),				// look = x*xprime + y*yprime + z*zprime
					            x*xprime.Y() + y*yprime.Y() + z*zprime.Y(),
								x*xprime.Z() + y*yprime.Z() + z*zprime.Z() );			// Get the look vector in heliodetic coordinates

//				ok1     = m_losGeometryIn[idx]->CreateRay( m_raytracingspecs->RayTracingShells(), quadraturespecs, m_location.Vector(), look, &z );
				std::unique_ptr< SKTRAN_RayOptical_Base>		ray;

				ok1     = m_diffuserayfactory->CreateRayObject     ( &ray      ); 
				ok1     = ok1 && ray->MoveObserver( m_location.Vector(), look  );
				ok1     = ok1 && ray->TraceRay_NewMethod();
				m_losGeometryIn[idx].AssignRay( std::move(ray) );
				ok      = ok && ok1;
				idx++;
			}
		}
		NXASSERT(( idx == m_incomingunitvectors->NumUnitVectors() ));
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING,"SKTRAN_DiffusePointGeometry_V21::ConfigureLOSArray, Error configuring the line of sight array"); 
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::ZenithAndAzimuthOfIncomingLos		2008-4-8*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_DiffusePointGeometry_V2::ZenithAndAzimuthOfIncomingLos( size_t idxoflosin, double* zenith, double*azimuth )	const
{
	bool	ok;
	size_t	izen;
	size_t	iazi;
	size_t	numzen;
	size_t	numazi;

	ok = (m_numlosIn > 0) && ( idxoflosin < m_numlosIn);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_DiffusePointGeometry_V21::ZenithAndAzimuthOfIncomingLos, You are indexing out of the valid bounds of the array");
		*zenith  = 9999.0;
		*azimuth = 9999.0;
	}
	else
	{
		numzen   = m_zenithsIn->NumAngles();											// Get the number of zeniths in the line of sight table
		numazi   = m_azimuthsIn->NumAngles();											// Get the number of azimuths in the line of sight table
		izen     = idxoflosin/numazi;
		iazi     = idxoflosin%numazi;
		*zenith  = m_zenithsIn->At (izen);
		*azimuth = m_azimuthsIn->At(iazi);
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::ZenithAndAzimuthOfIncomingLos		2008-4-8*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_DiffusePointGeometry_V2::ZenithAndAzimuthOfOutboundSource( size_t idxofsrc, double* zenith, double*azimuth )	const
{
	bool										ok;
	nxVector									v;
	const SKTRAN_UnitSphere_V2*					outboundsphere;

	outboundsphere = m_scatteringmatrix.OutboundUnitSphere();
	ok             = ( outboundsphere != NULL) && (idxofsrc < outboundsphere->NumUnitVectors() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_DiffusePointGeometry_V21::ZenithAndAzimuthOfOutboundSource, You are indexing out of the valid bounds of the array");
		*zenith  = 9999.0;
		*azimuth = 9999.0;
	}
	else
	{
		v        = outboundsphere->UnitVectorAt(idxofsrc);
		*zenith  = 90.0 - v.Latitude();
		*azimuth = v.Longitude();
	}
	return ok;
}
*/



/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::ConfigureGeometry		2007-12-12*/
/** The first stage of the point configuration. This configuration is
 *	performed in a single threaded envoironment and there are no
 *	syncronization issues.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointGeometry_V21::ConfigureGeometry_Stage1( std::shared_ptr<const SKTRAN_RayFactory_Base>&	  diffuserayfactory,
															    //const SKTRAN_SpecsInternal_RayTracing_V21*		  raytracingspecs,
																const  SKTRAN_SpecsInternal_Diffuse_V21*		  diffusespecs,
																const HELIODETIC_POINT&							  point,
																SKTRAN_GridIndex								  diffuseprofileindex,
																SKTRAN_GridIndex								  diffuseheightindex
															  )
{
	bool ok;
	const SKTRAN_UnitSphere_V2*	unitsphere;

	ReleaseResources();												// Release the current resources

	m_location              = point;
	m_diffuseheightindex    = diffuseheightindex;
	m_diffuseprofileindex   = diffuseprofileindex;
	m_diffuserayfactory     = diffuserayfactory;
	m_incomingunitvectors   = diffusespecs->IncomingUnitSphere( diffuseprofileindex, diffuseheightindex );
	m_ishighaltitudediffuse = (diffusespecs->MaxDiffuseAltitude(diffuseprofileindex) <= point.Altitude() );	
	unitsphere              = diffusespecs->OutboundUnitSphere(diffuseprofileindex, diffuseheightindex) ;

	ok =     ( m_incomingunitvectors != NULL) && (unitsphere != NULL);
	if (ok)
	{
		m_incomingunitvectors->AddRef();								// Hold onto the inocoming unit sphere
		m_scatteringmatrix.ConfigureGeometry_Stage1( m_location, m_incomingunitvectors, unitsphere );
	}
	else																// If it failed
	{																	// then log a message
		nxLog::Record(NXLOG_WARNING, "SKTRAN_DiffusePointGeometry_V21::ConfigureGeometry_Stage1, Error configuring zenith angles or LOS array");
		ReleaseResources();												// then release everything rathe rthan have things partially defined
	}																	// finally if the initialization failed
	return ok;															// return the status
}																		// and that is that.

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::ConfigureGeometry_Stage2MT		2007-12-12*/
/** The second stage of configuring the diffuse point. This is performed in a
 *	multi-threaded environment. This second stage allocates and configures
 *	all of the rays associated with the point. This code must be thread safe.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointGeometry_V21::ConfigureGeometry_Stage2MT( SKTRANSO_Quadrature_TLS_V21* quadrature )
{
	bool ok;

//	ok  =      (m_scatteringmatrix != NULL);
	ok  =       m_scatteringmatrix.ConfigureGeometry_Stage2MT();
	ok  = ok && ConfigureLOSArray(quadrature);							// Configure the corresponding line of sight arrays
	if (!ok)															// If it failed
	{																	// then log a message
		nxLog::Record(NXLOG_WARNING, "SKTRAN_DiffusePointGeometry_V21::ConfigureGeometry, Error configuring scattering matrix, zenith angles or LOS array");
		ReleaseResources();												// then release everything rathe rthan have things partially defined
	}																	// finally if the initialization failed
	return ok;															// return the status
}																		// and that is that.

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::IncomingSolidAngleArray		2009-1-21*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_DiffusePointGeometry_V2::IncomingSolidAngleArray( const double** domega, size_t* numangles )const
{

	*domega    = m_scatteringmatrix.SolidAngleArray();
	*numangles = m_scatteringmatrix.InboundUnitSphere()->NumUnitVectors();
	return true;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointGeometry_V21::CreateJIndexTables_MT		2008-2-6*/
/** This function creates all of the JIndexTables for each ray at this point
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffusePointGeometry_V21::CreateJIndexTables_MT( SKTRANSO_Quadrature_TLS_V21* quadrature  )
{
	bool									ok;
	bool									ok1;
	size_t									losCtr;
	SKTRANSO_RayInternalGeometry*			los;

	ok = true;
	for (losCtr = 0; losCtr < NumLOSIn(); ++losCtr)						// for all of the incoming
	{																	// lines of sight at this point
		los       = LOSAtVar(losCtr);									// Get the line of sight
		ok1       = los->IndicesVar()->CreateJIndexTables( quadrature, los, false );		// Link the tables, diffuse points are never single scatter
		ok        = ok && ok1;
	}																	// and scan over all elements and all lines of sight
	return ok;
}
