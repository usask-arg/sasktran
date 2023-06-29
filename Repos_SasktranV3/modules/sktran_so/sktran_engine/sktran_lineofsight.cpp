#include "../sasktranv21_internals.h"
#include <float.h>


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayInternalOptical::SKTRANSO_RayInternalOptical		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_RayInternalOptical::SKTRANSO_RayInternalOptical()
{
	m_geometryray       = NULL;
	m_totaltransmission = 0.0;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayInternalOptical::~SKTRANSO_RayInternalOptical		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_RayInternalOptical::~SKTRANSO_RayInternalOptical()
{
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayInternalOptical::AttachToGeometry		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_RayInternalOptical::AttachToGeometry	( const SKTRANSO_RayInternalGeometry* geometry )
{
	m_geometryray       = geometry;
	m_totaltransmission = 0.0;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayInternalDiffuseOptical_V2::LinkSSTables		2007-11-23*/
/** Link to the single scatter tables**/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_RayInternalOptical::ConfigureOptical( SKTRANSO_Quadrature_TLS_V21* quadrature, bool totaltransmissiononly, bool usecachedtransmission  )
{
	bool	ok;

	ok = quadrature->CalculateRayTransmission( this, &m_totaltransmission, totaltransmissiononly, usecachedtransmission);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayInternalDiffuseOptical_V2::SKTRAN_RayInternalDiffuseOptical_V2		2007-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayInternalDiffuseOptical_V2::SKTRAN_RayInternalDiffuseOptical_V2(  )
{
	m_internaldiffusegeometry = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayInternalDiffuseOptical_V2::~SKTRAN_RayInternalDiffuseOptical_V2		2007-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayInternalDiffuseOptical_V2::~SKTRAN_RayInternalDiffuseOptical_V2( )
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayInternalDiffuseOptical_V2::AttachToGeometry		2008-2-1*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayInternalDiffuseOptical_V2::AttachToGeometry( const SKTRANSO_RayInternalGeometry* los )
{
	bool	ok;

	NXTRACE_ONCEONLY( firsttime,("**** 2008-12-23 ***** SKTRAN_RayInternalDiffuseOptical_V2::AttachToGeometry, the solar transmission geometry can be fully completed inside configure otpical, saves a lot of space!!\n"));

	m_internaldiffusegeometry  = los;
	ok                         = SKTRANSO_RayInternalOptical::AttachToGeometry( los );
	m_singlescatterJ           = SKTRAN_DBL_TO_STOKES_SCALAR(0.0);

	ok = ok && m_diffuseJ.AttachToGeometry          ( los->Indices()->DiffuseJIndex() );
	ok = ok && m_groundpointMSJ.AttachToGeometry    ( los->Indices()->GroundPointMSJIndex() );
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR,"SKTRAN_RayInternalDiffuseOptical_V2::Allocate, Error allocating memory for optical ray");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical		2008-2-8*/
/** The reference ray (if not NULL) is a ray at the same zenith angle but different azimuth.
 *	Under these circulatances teh reference ray has pre-calculated the
 *	shell transmission buffer and the cell transmission buffer as we can take advantage
 *	of spherical symmetry.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical( SKTRANSO_Quadrature_TLS_V21*				quadrature,
														    bool									singlescatteronly,
															bool									usecachedtransmission,
															bool									usecachedcellfactors )
{
	bool										ok;
	bool										ok1;
	bool										groundishit;
	SKTRAN_StokesScalar							groundsignalss;
	SKTRAN_StokesScalar							atmosemission;
	SKTRAN_JValueTable_V21						Jtable;


	groundishit        = GeometryRay()->Ray()->Storage()->GroundIsHit();																	// We are good if we dont hit the ground.

	ok  =       SKTRANSO_RayInternalOptical::ConfigureOptical( quadrature,															// Configure the total transmission along teh ray
															false,																// request that we want more than just the total transmisison,
															usecachedtransmission);												// Let user decide if we can use pre-cached values for shell and total transmission
	if (!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical, FAILURE inside SKTRANSO_RayInternalOptical::ConfigureOptical");
	// --- evaluate the atmospheric single scatter terms
	ok1  = quadrature->CreateJValueTable_AtmosphericSingleScatter( this,														// Set up the atmospheric scatter term
																   m_internaldiffusegeometry->Indices()->SingleScatterJIndex(),			// convert the atmospheric single scatter JIndex
																   &Jtable,											// to a Jvalue array
																   true,														// We can definitely use cached transmission factors.
																   usecachedcellfactors);										// Let user decide if we can use cached cell factors

	if (!ok1) nxLog::Record(NXLOG_ERROR,"SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical, FAILURE inside quadrature->CreateJValueTable_AtmosphericSingleScatter");
	ok  = ok && ok1;
	m_singlescatterJ  = Jtable.Evaluate();																			// Immediately evaluate the single scatter term
	NXASSERT(( m_singlescatterJ >= 0.0 ));																						// and check for weird stuff

	// --- evaluate the atmospheric emissions
	ok1  = quadrature->CreateJValueTable_AtmosphericEmissions( this,															// Set up the atmospheric scatter term
															   m_internaldiffusegeometry->Indices()->AtmosphericEmissionJIndex(),			// convert the atmospheric single scatter JIndex
																&Jtable,											// to a Jvalue array
																true,														// We can definitely use cached transmission factors.
																true);														// We can definitely use the cached cell factorsLet user decide if we can use cached cell factors

	if (!ok1) nxLog::Record(NXLOG_ERROR,"SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical, FAILURE inside quadrature->CreateJValueTable_AtmosphericEmissions");
	ok  = ok && ok1;
	atmosemission = Jtable.Evaluate();																			// Immediately evaluate the single scatter term
	m_singlescatterJ += atmosemission;																			// Immediately evaluate the single scatter term
	NXASSERT(( atmosemission   >= 0.0 ));																						// and check for weird stuff
	NXASSERT(( m_singlescatterJ >= 0.0 ));																						// and check for weird stuff
	
	
	if (groundishit)																											// if we hit the ground
	{																															// then
	// --- evaluate the ground albedo, single scatterterms

		ok1  = ok && quadrature->CreateJValueTable_GroundSingleScatter(  this,													// convert the 
																		m_internaldiffusegeometry->Indices()->GroundPointSSJIndex(),		// the single scatter ground point table
																		&Jtable,												// into a JValue table
																		true);													// and we can use cached results
		if (!ok1) nxLog::Record(NXLOG_ERROR,"SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical, FAILURE inside quadrature->CreateJValueTable_GroundSingleScatter");

		groundsignalss    = Jtable.Evaluate();																					// Evaluate immediately
		m_singlescatterJ += groundsignalss;																						// and add it onto the single scatter term
		ok = ok && ok1 ;																										// keep track of the status
		NXASSERT(( groundsignalss >= 0.0 ));																					// Check for weird stuff
		NXASSERT(( m_singlescatterJ >= 0.0 ));

	// --- evaluate the ground emissions term
		ok1  = ok && quadrature->CreateJValueTable_GroundEmissions(  this,														// convert the 
																	 m_internaldiffusegeometry->Indices()->GroundPointEmissionJIndex(),		// the single scatter ground point table
																	 &Jtable,													// into a JValue table
																	 true);														// and we can use cached results
		if (!ok1) nxLog::Record(NXLOG_ERROR,"SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical, FAILURE inside quadrature->CreateJValueTable_GroundEmissions");

		groundsignalss    = Jtable.Evaluate();																					// Evaluate immediately
		m_singlescatterJ += groundsignalss;																						// and add it onto the single scatter term
		NXASSERT(( groundsignalss >= 0.0 ));																					// Check for weird stuff
		NXASSERT(( m_singlescatterJ >= 0.0 ));
		ok = ok && ok1 ;																										// keep track of the status
	}																										

    // ---- single scatter is flagged in the geometry as blank Jindex tables
	if (singlescatteronly)																											// If we have more than just single scatter
	{
		m_diffuseJ.Clear();
		m_groundpointMSJ.Clear();
	}
	else
	{																															// then 
		ok1  = quadrature->CreateJValueTable_AtmosphericDiffuseScatter( this,													// convert the
																		m_internaldiffusegeometry->Indices()->DiffuseJIndex(),				// diffuse JIndex
																		&m_diffuseJ,											// to a JValue table
																		true,													// Use cached ray transmissions buffer
																		true);													// Use cached cell transmissions

		if (!ok1) nxLog::Record(NXLOG_ERROR,"SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical, FAILURE inside quadrature->CreateJValueTable_AtmosphericDiffuseScatter");
		ok = ok && ok1;																											// and update status

		if (groundishit)																										// if we hit the ground
		{																														// then convert
			ok1  = quadrature->CreateJValueTable_InterpolateGroundDiffuseScatter( this,													// the multiple scatter 
																	   m_internaldiffusegeometry->Indices()->GroundPointMSJIndex(),		// off the ground		
																	   &m_groundpointMSJ,										// to a JValue table
																	   true);													// use cached ray transmissions
			if (!ok1) nxLog::Record(NXLOG_ERROR,"SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical, FAILURE inside quadrature->CreateJValueTable_InterpolateGroundDiffuseScatter");

			ok = ok && ok1;																										// and we are done
		}																														
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING," SKTRAN_RayInternalDiffuseOptical_V2::ConfigureOptical, Error configuring the JValue tables. Thats a problem");
		m_diffuseJ.Clear();
		m_groundpointMSJ.Clear();
		m_singlescatterJ = SKTRAN_DBL_TO_STOKES_SCALAR(-9.0E28);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayInternalDiffuseOptical_V2::CalculateTotalRadianceAtOrigin		2008-2-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayInternalDiffuseOptical_V2::CalculateTotalRadianceAtOrigin( SKTRAN_StokesScalar* radiance )
{
	SKTRAN_StokesScalar				diffuse;
	SKTRAN_StokesScalar				ground;

	diffuse    = m_diffuseJ.Evaluate();
	ground     = m_groundpointMSJ.Evaluate();
	NXASSERT(( diffuse >= 0.0 ));
	NXASSERT(( ground  >= 0.0 ));
	*radiance  = m_singlescatterJ + diffuse + ground;
	NXASSERT(( NXFINITE( *radiance ) ));

	return true;
}

