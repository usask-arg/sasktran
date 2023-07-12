#include "../sasktranv21_internals.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointGeometry_V21::SKTRAN_ScatteringMatrixPointGeometry_V21		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ScatteringMatrixPointGeometry_V21::SKTRAN_ScatteringMatrixPointGeometry_V21()
{
//	m_incomingzenith   = NULL;
//	m_incomingazimuth  = NULL;
	m_unitsphere         = NULL;
	m_incomingunitsphere = NULL;
//	m_numincominglos   = 0;
//	m_dOmega           = NULL;
//	m_losazimuth       = NULL;
//	m_loszenith        = NULL;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointGeometry_V21::ReleaseResources		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_ScatteringMatrixPointGeometry_V21::ReleaseResources()
{
//	if (m_incomingzenith  != NULL) m_incomingzenith->Release();
//	if (m_incomingazimuth != NULL) m_incomingazimuth->Release();
	if (m_unitsphere         != NULL) m_unitsphere->Release();
	if (m_incomingunitsphere != NULL) m_incomingunitsphere->Release();
//	if (m_dOmega             != NULL) delete [] m_dOmega;
//	if (m_losazimuth      != NULL) delete [] m_losazimuth;
//	if (m_loszenith       != NULL) delete [] m_loszenith;
//	m_incomingzenith   = NULL;
//	m_incomingazimuth  = NULL;
	m_unitsphere       = NULL;
//	m_numincominglos   = 0;
//	m_dOmega           = NULL;
//	m_losazimuth       = NULL;
//	m_loszenith        = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointGeometry_V21::~SKTRAN_ScatteringMatrixPointGeometry_V21		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ScatteringMatrixPointGeometry_V21::~SKTRAN_ScatteringMatrixPointGeometry_V21()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointGeometry_V21::ConfigureGeometry		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ScatteringMatrixPointGeometry_V21::ConfigureGeometry_Stage1( const HELIODETIC_POINT&                           location,
																		const SKTRAN_UnitSphereLatLonGrid*				  incomingunitsphere,
																		const SKTRAN_UnitSphere_V2*						  outboundunitsphere )

{
	bool ok;


	ReleaseResources();																					// Release the current resources

	m_location        = location;
//	m_incomingzenith  = incomingzenith;
//	m_incomingazimuth = incomingazimuth;
	m_unitsphere         = outboundunitsphere;
	m_incomingunitsphere = incomingunitsphere;
	ok = (m_incomingunitsphere != NULL)  && (m_unitsphere != NULL );
	if (ok)
	{
		m_incomingunitsphere->AddRef();
		m_unitsphere->AddRef();
	}
	else															// If it failed
	{																// then log a message
		nxLog::Record(NXLOG_WARNING, "SKTRAN_DiffusePointGeometry_V21::ConfigureGeometry, Error setting incoming and outbound angles, pointers are NULL, thats not good.");
		ReleaseResources();											// then release everything rathe rthan have things partially defined
	}																// finally if the initialization failed
	return ok;														// return the status
}																	// and that is that.

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointGeometry_V21::ConfigureGeometry		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ScatteringMatrixPointGeometry_V21::ConfigureGeometry_Stage2MT( )
{
//	bool	ok;

//	ok =       ConfigureSolidAngles();
//	ok = ok && ConfigureSinesAndCosines(  );

/*	if (!ok)														// If it failed
	{																// then log a message
		nxLog::Record(NXLOG_WARNING, "SKTRAN_DiffusePointGeometry_V21::ConfigureGeometry_Stage2, Error configuring zenith angles or solid angle array");
		ReleaseResources();											// then release everything rathe rthan have things partially defined
	}																// finally if the initialization failed
	return ok;														// return the status
*/
	return true;
}																	// and that is that.



/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointGeometry_V21::ConfigureSinesAndCosines		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_ScatteringMatrixPointGeometry_V2::ConfigureSinesAndCosines(  )
{
	size_t				numzenith;
	size_t				numazimuth;
	size_t				idx;
	bool				ok;

	// --- Calculate all of the the signals coming into this point.

	numzenith   = m_incomingzenith->NumAngles();				// The number of incoming zenith angles
	numazimuth  = m_incomingazimuth->NumAngles();				// The number of incoming azimuths
	ok          =       m_coszenith.SetSize ( numzenith  );
	ok          = ok && m_sinzenith.SetSize ( numzenith  );
	ok          = ok && m_cosazimuth.SetSize( numazimuth );
	ok          = ok && m_sinazimuth.SetSize( numazimuth );

	if (ok)
	{
		for (idx = 0; idx < numzenith; idx++ )
		{
			m_coszenith.At(idx) = cos(m_incomingzenith->At(idx) );
			m_sinzenith.At(idx) = sin(m_incomingzenith->At(idx) );
		}

		for (idx = 0; idx < numazimuth; idx++ )
		{
			m_cosazimuth.At(idx) = cos(m_incomingazimuth->At(idx) );
			m_sinazimuth.At(idx) = sin(m_incomingazimuth->At(idx) );
		}
	}
	else
	{
		m_coszenith.erase();
		m_sinzenith.erase();
		m_cosazimuth.erase();
		m_sinazimuth.erase();
		nxLog::Record(NXLOG_WARNING,"SKTRAN_ScatteringMatrixPointGeometry_V21::ConfigureSinesAndCosines(  ), Error allocating memory for sines and cosines, thats a problem");
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointGeometry_V21::ConfigureSolidAngles		2008-1-10*/
/** Caches the array of solid angles associated with the incoming rays.
 *	The solid angle is a quadrature weight applied to each incoming ray in
 *	addition to the scattering coefficients. The solid angle is based upon the
 *	area surrounding a ray: the ray is centered in the middle of the area.
 *
 *	The incoming rays should not be in the zenith direction as this will result
 *	in a solid angle weight of 0.0, implying this ray has no impact upon the 
 *	the scattering calculation.  
 *
 *	I do calculate the solid angle weights so they add up to exactly 2 pi for ground points
 *	and exactly 4 pi for atmospheric points.  This is done by using mu space (= cos(zenithh angle))
 *	space rather than zenith angle and making sure we cover the entire sphere by
 *	padding the ends of the azimuth and zenith arrays with proper end point
 *	values. I do assume the zenith and azimuth arrays are in sorted, ascending
 *	order (ie zenith/azimuth start at angle close to 0 and finish at angles closes to pi/2pi )
 **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_ScatteringMatrixPointGeometry_V2::ConfigureSolidAngles(  )
{
	size_t				numincominglos;
	size_t				i;

	numincominglos = m_incomingunitsphere->NumUnitVectors();
	m_dOmega       = new double [numincominglos];
	for (i = 0; i < numincominglos; i++)
	{
		m_dOmega[i] = m_incomingunitsphere->CubatureWeightAt(i);
	}
	return true;
}
*/


/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::SKTRAN_ScatteringMatrixPointOptical_V21		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ScatteringMatrixPointOptical_V21::SKTRAN_ScatteringMatrixPointOptical_V21()
{
	m_scatmatrixgeometry = NULL;
	m_scatterarray       = NULL;
	m_numincoming        = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::ReleaseResources		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_ScatteringMatrixPointOptical_V21::ReleaseResources()
{
	if (m_scatterarray != NULL) delete [] m_scatterarray;
	m_scatterarray = NULL;
	m_numincoming = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::~SKTRAN_ScatteringMatrixPointOptical_V21		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_ScatteringMatrixPointOptical_V21::~SKTRAN_ScatteringMatrixPointOptical_V21()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::SetGeometry		2010-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_ScatteringMatrixPointOptical_V21::SetGeometry( const SKTRAN_ScatteringMatrixPointGeometry_V21* geometry)
{
	m_scatmatrixgeometry = geometry;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::NumIncoming		2011-6-29*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_ScatteringMatrixPointOptical_V21::NumIncoming	() const
{
	return m_scatmatrixgeometry->InboundUnitSphere()->NumUnitVectors();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::NumOutgoing		2011-6-29*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_ScatteringMatrixPointOptical_V21::NumOutgoing() const
{
	return m_scatmatrixgeometry->OutboundUnitSphere()->NumUnitVectors();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::AttachToGeometry		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ScatteringMatrixPointOptical_V21::AttachToGeometry()
{
	size_t	numscat;
	bool	ok;
	size_t	numincoming;
	size_t	numoutgoing;

	NXASSERT(( m_scatmatrixgeometry != NULL ));

	ReleaseResources();
	numincoming   = m_scatmatrixgeometry->InboundUnitSphere()->NumUnitVectors();
	numoutgoing   = m_scatmatrixgeometry->OutboundUnitSphere()->NumUnitVectors();
	numscat       = numincoming * numoutgoing;
	ok = (numscat == 0);				// See if there is anything to do.
	if (!ok)							// There is so allocate memory
	{
		m_scatterarray = new SKTRAN_PhaseMatrixScalar[numscat];
		ok = (m_scatterarray != NULL);
		if (!ok)
		{
				nxLog::Record(NXLOG_WARNING, "SKTRAN_ScatteringMatrixPointOptical_V21::AttachToGeometry, Error allocating memory for %Iu elements in scatter matrix.  You might be using too much memeory", (size_t)numscat );
		}
	}
	if (!ok) ReleaseResources();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::ComputeMultipliersAndAdjustScatterArray		2013-1-10*/
/** Check conservation of photons.  For each incoming ray, check total outbound
 *	radiance.  Note the ratio of integrated outbound radiance for each incoming
 *	rad.  If multipliers are outisde the range [1e-1,1e1], compute the mean phase
 *	matrix over the incoming ray that is causing the trouble.  This provides an 
 *	order-of-magnitude improvement in the multiplier for these 'grazing' rays.**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ScatteringMatrixPointOptical_V21::ComputeMultipliersAndAdjustScatterArray( const SKTRAN_TableOpticalProperties_V21* opticalprops )
{
	bool										ok = true;
	double										w;
	size_t										inidx;
	size_t										outidx;
	SKTRAN_StokesScalar							ratio,ratio1;					// factor to check conservation of photons 
	const SKTRAN_UnitSphere_V2*					unitsphere;
	size_t										numincoming;
	size_t										numoutgoing;
	double										scatextinction;

	unitsphere	  = m_scatmatrixgeometry->OutboundUnitSphere();
	numincoming   = m_scatmatrixgeometry->InboundUnitSphere()->NumUnitVectors();
	numoutgoing   = m_scatmatrixgeometry->OutboundUnitSphere()->NumUnitVectors();
	scatextinction= opticalprops->ScatteringExtinctionPerCM( m_scatmatrixgeometry->Location() );

	#if defined(NXDEBUG)
		double				wsum;
	#endif

	NXASSERT(( scatextinction > 0.0 ));
	for( inidx = 0; inidx < numincoming; ++inidx )							// check photon conservation on scatter:
	{																			// for each incoming ray,
		ratio   = 0.0;															// integrate the (scattered) source functions (vol. emission)
		#if defined(NXDEBUG)
			wsum = 0.0;
		#endif
		for( outidx = 0; outidx < numoutgoing; ++outidx )									// from every outgoing point
		{																					// and compare it to the incoming radiance's vol. emission
			w	   = unitsphere->CubatureWeightAt(outidx);									// 
			ratio += m_scatterarray[inidx+numincoming*outidx]*w;							// if photons are perfectly conserved,
			#if defined(NXDEBUG)
				wsum += w;
			#endif
		}
		ratio /= m_scatmatrixgeometry->InboundUnitSphere()->CubatureWeightAt(inidx);		// divide out theincoming solid angle array 
		ratio /= scatextinction;															// this ratio: integrated o/b to this i/c rad

		if ( (ratio < pow(10,-0.75)) || (ratio > pow(10,0.75)) )
		{
			ratio1 = ComputeMeanPhaseFunctionAndAdjustScatterArray( inidx,opticalprops );
			static int loopcount = 0;
			if (loopcount < 10)
			{
				nxLog::Record(NXLOG_WARNING, "SKTRAN_ScatteringMatrixPointOptical_V2::ComputeMultipliersAndAdjustScatterArray, There is a big scattering ratio adjustment for photon conservation (%e). Computed mean phase function, changing this to improved ratio (%e).", (double)ratio,(double)ratio1);
				loopcount++;
			}
			ratio = ratio1;											// update ratio also after taking mean
		}
		NXASSERT(( ratio > 0.0 ));
		NXASSERT(( (wsum > 12.566) && (wsum < 12.567) ));
		for( outidx = 0; outidx < numoutgoing; ++outidx )						
		{												
			m_scatterarray[inidx+numincoming*outidx] /= ratio;
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V2::ConfigureOptical		2008-1-10*/
/** Evaluates the scattering matrix geometry. Gets the scattering cross-section
 *	times the solid angle for each incoming ray direction
 **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_ScatteringMatrixPointOptical_V21::ConfigureOptical	( const SKTRAN_TableOpticalProperties_V21* opticalprops )
{
	bool												ok = true;
	size_t												outidx;
	size_t												inidx;
	size_t												numincoming;
	double												cosangle;
	double												lx,ly,lz;
	double												x,y,z;
	nxVector											outlook;
	size_t												scatidx;
	const SKTRAN_UnitSphere_V2*							unitsphere;
	const SKTRAN_UnitSphereLatLonGrid*					incomingsphere;
	SKTRANSO_JIndex									spatialindex[8];
	size_t												numspatial;
	const nxVector*										v;
	double												scatextinction;
	size_t												numoutgoing;

//	NXASSERT(( (m_scatterarray != NULL) ));		// Make sure the user has made a successful call to AttachGeometry 

	unitsphere      = m_scatmatrixgeometry->OutboundUnitSphere();
	incomingsphere  = m_scatmatrixgeometry->InboundUnitSphere();
	numincoming     = incomingsphere->NumUnitVectors(); 
	numoutgoing     = unitsphere->NumUnitVectors();

	scatextinction	= opticalprops->ScatteringExtinctionPerCM( m_scatmatrixgeometry->Location() );

	numspatial = N_ELEMENTS(spatialindex);
	opticalprops->GetBoundingSpatialPoints(m_scatmatrixgeometry->Location(), spatialindex, &numspatial );

	scatidx = 0;																// reset the loop counter for the scattering coefficients array
	for( outidx = 0; outidx < numoutgoing; ++outidx )							// loop over all of the outbound directions
	{																			// for each outbound direction
		outlook  = unitsphere->UnitVectorAt(outidx);									// Get the direction of the outbound vector in local coordinates
		lx       = outlook.X();															// and store locally as lx, ly, lz
		ly       = outlook.Y();
		lz       = outlook.Z();
		for(  inidx = 0; inidx < numincoming; ++inidx )										// loop slowly over zenith angle
		{																				// for each zenith
			v         = &incomingsphere->UnitVectorAt(inidx);
			x         = v->X();															// Get the X component of the incoming ray in local coords
			y         = v->Y();															// Get the Y component of the incoming ray in local coords
			z         = v->Z();
            cosangle  = -x*lx - y*ly - z*lz;											// Get the dot product to get the cosine of scattering angle
			ok        = ok && opticalprops->GetScatteringCoefficientCM2ForJindex( &spatialindex[0], numspatial, cosangle, &m_scatterarray[scatidx] );
			m_scatterarray[scatidx] *= SKTRAN_DBL_TO_STOKES_SCALAR( incomingsphere->CubatureWeightAt(inidx) );	// Multiply the phase matrix by the solid angle
			scatidx++;
		}																				// and point to next scattering entry and solid angle entry
		NXASSERT(( inidx == numincoming ));
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_ScatteringMatrixPointOptical_V21::ConfigureOptical, There was an error configuring the optical properties");
	}
	ComputeMultipliersAndAdjustScatterArray( opticalprops );										// get the factors: 

	NXASSERT(( scatidx == (numincoming * numoutgoing) ));
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::ScatterIncomingRays		2008-1-10*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_ScatteringMatrixPointOptical_V21::ScatterIncomingRays( SKTRAN_DiffusePointOptical_V21* point )
{
	size_t						inidx;
	size_t						scatidx;
	size_t						outidx;
	SKTRAN_StokesScalar			scatter;
	SKTRAN_StokesScalar*		scatptr;
	const SKTRAN_StokesScalar*	incomingradiance;
	size_t						numincoming;
	size_t						numoutgoing;

	numincoming   = m_scatmatrixgeometry->InboundUnitSphere()->NumUnitVectors();
	numoutgoing   = m_scatmatrixgeometry->OutboundUnitSphere()->NumUnitVectors();

	NXASSERT(( numincoming == point->NumIncomingLos() ));
	NXASSERT(( numoutgoing == point->NumOutgoingRad() ));

	incomingradiance = point->IncomingRadianceArray();

	NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_ScatteringMatrixPointOptical_V21::ScatterIncomingRays, think about loop unrolling the tight innner loop\n"));

	scatidx = 0;																// reset the loop counter for the scattering coefficients array
	for( outidx = 0; outidx < numoutgoing; ++outidx )							// loop over all of the outbound directions
	{																			// for each outbound direction
		scatter = 0.0;
		for(  inidx = 0; inidx < numincoming; ++inidx )							// loop slowly over zenith angle
		{																		// for each zenith
			scatter +=  incomingradiance[inidx] * m_scatterarray[scatidx];		// multiply incoming radiance by scattering coefficient and solid angle and then sum
			scatidx++;
		}																		// do all of the incoming
		NXASSERT((scatter >= 0.0));
		scatptr = point->OutgoingJPtr(outidx);
		*scatptr = scatter;
	}
	return true;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21::ComputeMeanPhaseFunctionAndAdjustScatterArray		2013-1-10*/
/**	For the specified 'in' direction, set up high-resolution points within its
 *	associated solid angle and, for each outbound direction, compute a mean value 
 *	to use for the phase function using all high-res incoming points. **/
/*---------------------------------------------------------------------------*/

SKTRAN_StokesScalar SKTRAN_ScatteringMatrixPointOptical_V21::ComputeMeanPhaseFunctionAndAdjustScatterArray( size_t inidx,const SKTRAN_TableOpticalProperties_V21* opticalprops )
{
	SKTRAN_StokesScalar						ratio;
	const SKTRAN_UnitSphere_V2*				unitsphere;
	const SKTRAN_UnitSphereLatLonGrid*		incomingsphere;
	SKTRAN_PhaseMatrixScalar				scatter;
	SKTRAN_PhaseMatrixScalar				totalscatter;
	SKTRANSO_JIndex						spatialindex[8];
	size_t									numspatial,numoutgoing,numincoming,idxzen,idxazi,outidx,numsmallzen,numsmallazi;
	double									zen,azi,cosangle,cosangle1;
	double									w,cosanglethresh,scatextinction,inzen,inazi,domega;
	double									dtheta,dphi,domegaEst,deltaAngle,x,y,z,xp,yp,zp,hp;
	const nxVector*							v;
	nxVector								outlook;
	bool									ok(true);

	incomingsphere	= m_scatmatrixgeometry->InboundUnitSphere();
	unitsphere		= m_scatmatrixgeometry->OutboundUnitSphere();
	numincoming     = incomingsphere->NumUnitVectors(); 
	numoutgoing     = unitsphere->NumUnitVectors();
	scatextinction	= opticalprops->ScatteringExtinctionPerCM( m_scatmatrixgeometry->Location() );
	numspatial		= N_ELEMENTS(spatialindex);
	opticalprops->GetBoundingSpatialPoints( m_scatmatrixgeometry->Location(),spatialindex,&numspatial );
	
	cosanglethresh	= nxmath::cosd( 3.0 );									// define 'small angle' threshold
	deltaAngle		= 0.007*nxmath::Pi/180.0;								// angular resolution for surface averaging
	v = &incomingsphere->UnitVectorAt(inidx);

	incomingsphere->LocalLookToAziZen( *v,&inazi,&inzen );					// incoming vector info 
	domega			= incomingsphere->CubatureWeightAt( inidx );			// 
	dtheta			= nxmath::Pi/double(incomingsphere->NumZenith());		// approx zenith size of domega
	dphi			= 2*nxmath::Pi/double(incomingsphere->NumAzimuthVertex());// approx azi size of domega
	if (inzen < 1e-4 )
	{
		domegaEst		= sin(inzen+0.5*dtheta)*dphi*dtheta;				// make est. solid angle finite
		numsmallzen		= (size_t)floor(0.5*dtheta/deltaAngle);				// polar: smaller zen size
	}
	else if (fabs(inzen-nxmath::Pi)<1e-4)
	{
		domegaEst		= sin(inzen-0.5*dtheta)*dphi*dtheta;
		numsmallzen		= (size_t)floor(0.5*dtheta/deltaAngle);				// polar: smaller zen size
	}
	else
	{
		domegaEst		= sin(inzen)*dphi*dtheta;							// est. solid angle for this inc ray
		numsmallzen		= (size_t)floor(dtheta/deltaAngle);					// number of zenith rays for interp.
	}
	numsmallazi		= (size_t)floor(dphi/deltaAngle);				// number of azimuth rays for interp.
	ratio			= 0.0;
	for ( outidx=0;outidx<numoutgoing;outidx++ )
	{
		outlook  = unitsphere->UnitVectorAt(outidx);
		x	= outlook.X();
		y	= outlook.Y();
		z	= outlook.Z();
		cosangle = x*(-(*v).X()) + y*(-(*v).Y()) + z*(-(*v).Z());
		if (cosangle >= cosanglethresh)						// if this scattering angle is small,
		{													// do high-res points
			totalscatter	= 0.0;
			for ( idxzen=0;idxzen<numsmallzen;idxzen++ )			// loop slowly over zen surrounding
			{														// the incoming vector
				zen	= inzen - 0.5*dtheta + double(idxzen)*deltaAngle;// 
				zp	= cos( zen );
				hp	= sin( zen );
				for ( idxazi=0;idxazi<numsmallazi;idxazi++ )		// loop quickly over azimuths
				{
					azi	= inazi - 0.5*dphi + double(idxazi)*deltaAngle;
					xp	= cos( azi )*hp;
					yp	= sin( azi )*hp;
					cosangle1	= x*(-xp) + y*(-yp) + z*(-zp);		// eval phase funct at point
					ok = ok && opticalprops->GetScatteringCoefficientCM2ForJindex( spatialindex,numspatial,cosangle1,&scatter );
					totalscatter	+= scatter;						// estimate mean phase funct
				}
			}
			if (numsmallzen*numsmallazi!=0)
			{
				totalscatter /= double(numsmallzen*numsmallazi);
			}
			m_scatterarray[inidx+numincoming*outidx] = totalscatter*SKTRAN_DBL_TO_STOKES_SCALAR( domega*(domega/domegaEst) );// set to adjusted (mean) value
		}
		w	   = unitsphere->CubatureWeightAt(outidx);			// get cub weight and find 
		ratio += m_scatterarray[inidx+numincoming*outidx]*w;			// new estimated ratio
	}
	ratio /= (domega*scatextinction);							// divide out theincoming solid angle array 
	return ratio;
}

bool SKTRAN_ScatteringMatrixPointOptical_V21::ScatterIncomingRays( SKTRAN_DiffusePointOptical_V21* point )
{

	size_t							numindices;
	size_t							numloops;
	size_t							numrem;
	size_t							idx;
	size_t							scatidx;
	size_t							outidx;
	size_t							l;

	const SKTRAN_StokesScalar*		incomingradiance;
	SKTRAN_PhaseMatrixScalar*		scattarrayptr;
	SKTRAN_PhaseMatrixScalar*		s;
	SKTRAN_StokesScalar				sum;
	const SKTRAN_StokesScalar*		r;
	SKTRAN_StokesScalar*			scatptr;
	size_t				  		    numincoming;
	size_t							numoutgoing;


//	NXASSERT(( m_numincoming == point->NumIncomingLos() ));
//	NXASSERT(( m_numoutgoing == point->NumOutgoingRad() ));

	incomingradiance = point->IncomingRadianceArray();
	numincoming = NumIncoming();
	numoutgoing = NumOutgoing();
	numindices = numincoming;								// So get the number of incoming radainces .
	numloops   = numincoming/16;							// break it into chunks of 16
	numrem     = numincoming%16;							// but dont forget the remainder

	scatidx = 0;											// reset the loop counter for the scattering coefficients array
	for( outidx = 0; outidx < numoutgoing; ++outidx )		// loop over all of the outbound directions
	{														// for each outbound direction
		sum           = 0.0;
		scattarrayptr =  m_scatterarray +scatidx;			// Get the array of scattering coefficient
		idx = 0;											// set the index into the radiances and weights
		for (l = 0; l < numloops; l++ )						// for the twelve point loop
		{													// loop over the 12 points
			r    = incomingradiance + idx;		// get the radiance pointer
			s    = scattarrayptr    + idx;		// get the weight points
			sum +=   (r[ 0])*s[ 0]				// add the product of the radiance times the weight
				+    (r[ 1])*s[ 1]				// Do all twelve. This is essentially loop unrolling
				+    (r[ 2])*s[ 2]				// and helps the compiler speed up the processor
				+    (r[ 3])*s[ 3]				// by avoiding excessive loops.
				+    (r[ 4])*s[ 4]				// We might even be able to use SSE/MMX instructions in this region
				+    (r[ 5])*s[ 5]
				+    (r[ 6])*s[ 6]
				+    (r[ 7])*s[ 7]
				+    (r[ 8])*s[ 8]
				+    (r[ 9])*s[ 9]
				+    (r[10])*s[10]
				+    (r[11])*s[11]
				+    (r[12])*s[12]
				+    (r[13])*s[13]
				+    (r[14])*s[14]
				+    (r[15])*s[15];				// do all 16
			idx +=  16;
		}
		
		for (l = 0; l < numrem; l++ )									// For the last few
		{																// just
			sum +=  (incomingradiance[idx])*scattarrayptr[idx];			// do the calculation
			idx++;														// and that is that
		}

		scatptr  = point->OutgoingJPtr(outidx);
		*scatptr = sum;

		NXASSERT(( idx == numincoming ));
		scatidx += numincoming;										// Point to the next set of scattering matrices in teh square array
	}																		// do all of the incoming
	return true;
}

