#include "../sasktranv21_internals.h"
#include "../sasktran_legacy21.h"

/*-----------------------------------------------------------------------------
 *					class SKTRAN_Quadrature_Factory_Legacy_V21		2010-5-17*/
/** A factory object used to create new instances of quadrature obejcts
 *	for each thread.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_Quadrature_Factory_Legacy_V21 : public SKTRAN_Quadrature_Factory_V21
{
	public:
		virtual			   ~SKTRAN_Quadrature_Factory_Legacy_V21(){};
		virtual bool		CreateThreadInstance( const SKTRAN_SpecsInternal_V21*		modelspecifications,
												  const SKTRAN_EngineDiffuseTables*		modeltables,
												  SKTRANSO_Quadrature_TLS_V21**			instance ) const; 
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsInternal_Quadrature_V21_Legacy 2010-5-17*/
/** The internal quadrature specifications used by the engine. The quadrature
 *	specs just keeps a pointer to the factory object that can create a new
 *	quadrature instance for each processing thread within the sasktran engine.
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_Quadrature_V21_Legacy : public SKTRAN_SpecsInternal_Quadrature_V21
{
	private:
		SKTRAN_Quadrature_Factory_Legacy_V21*					m_factory;

	public:
																SKTRAN_SpecsInternal_Quadrature_V21_Legacy	();
		virtual												   ~SKTRAN_SpecsInternal_Quadrature_V21_Legacy	();
		virtual const SKTRAN_Quadrature_Factory_V21*				QuadratureFactory							() const { return m_factory; }

};

#include "sktran_quadrature_tls_v2_legacy_header.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Quadrature_V21_Legacy::CreateInternalSpecs		2010-5-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_Quadrature_V21_Legacy::CreateInternalSpecs( const SKTRAN_SpecsInternal_Quadrature_V21** internalspecs ) const
{
	bool	ok;

	*internalspecs = new SKTRAN_SpecsInternal_Quadrature_V21_Legacy();
	ok = (*internalspecs != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_SpecsUser_Quadrature_V21_Legacy::CreateInternalSpecs, Error allocating memory");
	}
	else
	{
		(*internalspecs)->AddRef();
	}
	return ok;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Quadrature_V21_Legacy::SKTRAN_SpecsInternal_Quadrature_V21_Legacy		2010-5-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_Quadrature_V21_Legacy::SKTRAN_SpecsInternal_Quadrature_V21_Legacy()
{
	m_factory = new SKTRAN_Quadrature_Factory_Legacy_V21;
	m_factory->AddRef();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Quadrature_V21_Legacy::~SKTRAN_SpecsInternal_Quadrature_V21_Legacy		2010-5-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_Quadrature_V21_Legacy::~SKTRAN_SpecsInternal_Quadrature_V21_Legacy()
{
	m_factory->Release();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_Quadrature_Factory_Legacy_V21::CreateThreadInstance		2010-5-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Quadrature_Factory_Legacy_V21::CreateThreadInstance( const SKTRAN_SpecsInternal_V21*		modelspecifications,
																 const SKTRAN_EngineDiffuseTables*		modeltables,
																 SKTRANSO_Quadrature_TLS_V21**				instance ) const
{
	bool ok;
	SKTRANSO_Quadrature_TLS_V2_Legacy*	quadrature;

	quadrature = new SKTRANSO_Quadrature_TLS_V2_Legacy( modelspecifications, modeltables );
	ok         = (quadrature != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_Quadrature_Factory_Legacy_V21::CreateThreadInstance, Error allocating thread instance, thats a problem");
	}
	else
	{
		quadrature->AddRef();
	}
	*instance  = quadrature;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::SKTRANSO_Quadrature_TLS_V2_Legacy		2010-3-8*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_Quadrature_TLS_V2_Legacy::SKTRANSO_Quadrature_TLS_V2_Legacy( const SKTRAN_SpecsInternal_V21*		modelspecifications,
																  const SKTRAN_EngineDiffuseTables*	    modeltables)
{
	bool	ok;

	m_cellfactors                 = NULL;
	m_shelltransmissionsbuffer    = NULL;
	m_scratchbuffer               = NULL;
	m_maxshellsalongray           = 0;
	m_maxquadpointsalongray       = 0;
	m_opticalprops                = NULL;
	m_albedobrdf                  = NULL;
	m_enginetables                = modeltables;
	ok = Initialize( modelspecifications);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRANSO_Quadrature_TLS_V2_Legacy::Constructor, Error initializing the thread dependent quadrature object, That might be a problem");
	}

}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::~SKTRANSO_Quadrature_TLS_V2_Legacy		2010-3-8*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_Quadrature_TLS_V2_Legacy::~SKTRANSO_Quadrature_TLS_V2_Legacy()
{
	ReleaseTemporaryWorkAreas();
	if (m_opticalprops != NULL ) m_opticalprops->Release();
	if (m_albedobrdf   != NULL ) m_albedobrdf->Release();

}



/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineOptical_V2::ReleaseTemporaryWorkAreas		2009-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_Quadrature_TLS_V2_Legacy::ReleaseTemporaryWorkAreas()
{
	if ( m_cellfactors               != NULL ) delete [] m_cellfactors;
	if ( m_shelltransmissionsbuffer  != NULL ) delete [] m_shelltransmissionsbuffer;
	if ( m_scratchbuffer             != NULL ) delete [] m_scratchbuffer;

	m_cellfactors                 = NULL;
	m_shelltransmissionsbuffer    = NULL;
	m_scratchbuffer               = NULL;
	m_maxshellsalongray           = 0;
	m_maxquadpointsalongray       = 0;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::MaxQuadraturePointsAlongRay		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

size_t  SKTRANSO_Quadrature_TLS_V2_Legacy::MaxQuadraturePointsAlongRay( const SKTRAN_SpecsInternal_RayTracing_V21* raytracingspecs)
{
	size_t	maxshells;
	size_t	quadpoints;

	maxshells   = raytracingspecs->MaxShellsAlongRay();
	quadpoints  = MaxQuadraturePointsPerCell();
	return (maxshells*quadpoints);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineOptical_V2::AllocateTemporaryWorkAreas		2009-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::Initialize( const SKTRAN_SpecsInternal_V21* modelspecifications )

{
	bool	ok;
	size_t											maxshellsalongray;
	size_t											maxQ;
	size_t											maxpoints;
	size_t											maxindices;
	const SKTRAN_SpecsInternal_RayTracing_V21*		raytracingspecs; 
	const SKTRAN_AlbedoBRDF_V2*						albedobrdf;

	albedobrdf       = modelspecifications->GroundPointSpecs()->AlbedoBRDF();
	raytracingspecs  = modelspecifications->RayTracingSpecs();
	m_coordinates      = modelspecifications->CoordinateSystemObject();


	if (albedobrdf   != NULL) albedobrdf->AddRef();
	if (m_albedobrdf != NULL) m_albedobrdf->Release();
	m_albedobrdf          = albedobrdf;



	maxshellsalongray     = raytracingspecs->MaxShellsAlongRay();
	maxQ                  = MaxQuadraturePointsAlongRay( raytracingspecs );
	maxpoints             = m_enginetables->DiffusePointsTable()->MaxInterpolatedPointsPerJ();
	maxindices            = maxQ*maxpoints;


	ok =    ( maxshellsalongray <= m_maxshellsalongray )						// We dont need to re-allocate if it is less than we need
		 && ( maxQ              <= m_maxquadpointsalongray);				// check both the 
	if (!ok)
	{
		ReleaseTemporaryWorkAreas();
		maxshellsalongray     += 50;
		maxQ                 += 100;

		m_shelltransmissionsbuffer  = new double [maxshellsalongray];
		m_scratchbuffer             = new double [maxQ];
		m_cellfactors               = new double [maxQ];

		ok =    ( m_shelltransmissionsbuffer != NULL )
			 && ( m_scratchbuffer            != NULL )
			 && ( m_cellfactors              != NULL );
		if (ok)
		{
			m_maxshellsalongray     = maxshellsalongray;
			m_maxquadpointsalongray = maxQ;
		}
	}
	ok  = ok && m_jindexworkspace.AllocateMaximumStorage( maxshellsalongray, maxQ, maxindices );
//	ok  = ok && m_jtableworkarea.ReserveStorage         ( maxindices );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Quadrature_TLS_V2_Legacy::ConfigureThreadLocalStore, Error allocating storage for %u shells and %u quad points along ray, This is a problem", (size_t)maxshellsalongray, (size_t)maxQ);
		ReleaseTemporaryWorkAreas();
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::ConfigureOptical		2010-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::SetOpticalProps( const SKTRAN_TableOpticalProperties_V21*	optprop)					//!< The optical properties table;
{
	bool	ok;

	if (optprop        != NULL ) optprop->AddRef();
	if (m_opticalprops != NULL ) m_opticalprops->Release();
	m_opticalprops  = optprop;

	ok = (m_opticalprops != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Quadrature_TLS_V2_Legacy::ConfigureOptical, There was an error during ConfigureOptical. Thats a problem");
	}
	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_QuadratureSourceFunction_HomogenousCell_InternalDiffuse_V2::GetCellQuadraturePointRadiusSzaAndLookTowardsObserver		2008-2-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::GetCellQuadraturePointAltitudeSzaAndLookTowardsObserver( const SKTRAN_RayStorage_Base* ray, size_t cellidx, HELIODETIC_POINT* point, HELIODETIC_UNITVECTOR* look, double* celllengthmeters ) const
{
	bool								ok;

	ok = ray->CellMidPoint( cellidx, point );
	*celllengthmeters = ray->CellLength( cellidx);
	*look             =  ray->AverageLookVectorTowardsObserver(cellidx);

	if (!ok)
	{
		point->Clear();
		look->SetCoords(0,0,0);
		*celllengthmeters = 0.0;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_QuadratureSourceFunction_HomogenousCell_InternalDiffuse_V2::GeometryJIndexTableFromQuadraturePointsAlongRay		2008-2-9*/
/** Calculate the Jindex that represents the source functions evaluated along
 *	the ray at the quadrature points. This calculation is common to both
 *	single scatter and diffuse terms and is given by
 *	
 *	\f$ \sum_{k=0}^{11} w_{k}J_{D,k} \f$
**/
/*---------------------------------------------------------------------------*/
bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTableForQuadraturePointsAlongRay( const SKTRANSO_JindexTableBase* tablesource, const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex )
{
	bool								ok;
	bool								ok1;
	size_t								numquadrature;
	size_t								quadidx;
	SKTRANSO_JIndex						jdescriptor[100];
	size_t								nj;
	HELIODETIC_UNITVECTOR				look;
	double								length;
	HELIODETIC_POINT					heliolocation;

	m_jindexworkspace.ResetCounters();																			// Reset the counters in this temporary storage area
	numquadrature = ray->Storage()->NumCells();
	ok = (numquadrature == 0);
	if (!ok)
	{
		ok = true;
		for (quadidx=0; quadidx < numquadrature; quadidx++ )
		{
			ok1 = GetCellQuadraturePointAltitudeSzaAndLookTowardsObserver( ray->Storage(), quadidx, &heliolocation, &look, &length);			// Get info about the one quadrature point for this cell
			ok1 = ok1 && tablesource->InterpolateTable( heliolocation, look, jdescriptor, N_ELEMENTS(jdescriptor), &nj, 1.0);		// Look up the J index descriptors for this quadrature point
			ok1 = ok1 && m_jindexworkspace.InsertQuadraturePointEntries( jdescriptor, nj);											// Add these J index descriptots for the quadrature
			ok  = ok  && ok1;
			m_jindexworkspace.FinishCellEntries();
		}
	}
	ok = ok && jindex->DeepCopy( m_jindexworkspace );
	if (!ok)
	{
		 jindex->Clear();
		 nxLog::Record(NXLOG_WARNING, "SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTableForQuadraturePointsAlongRay, error creating Jindex table. Thats a problem");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_QuadratureSourceFunctionHomogenousCell_V2::GeometryJIndexTableFromGroundPoint		2008-2-9*/
/** Indexes either the table of diffuse groundpoints for multiple scatter ground
 *	points or the solar transmission table for sinngle scatter points**/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTableForGroundPoint( const SKTRANSO_JindexTableBase* tablesource, const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex, bool isSingleScatterTerm )
{
	bool								ok;
	HELIODETIC_UNITVECTOR				look;
	SKTRANSO_JIndex					jdescriptor[100];
	size_t								nj;
	double								coszen;
	HELIODETIC_POINT					location;
//	HELIODETIC_VECTOR					groundgeo;
	double								R;				// The reflectance function
	double								cossza;

	m_jindexworkspace.ResetCounters();																	// Reset the counters in this temporary storage area
	ok = (!ray->Storage()->GroundIsHit());																	// it only makes sense if ground is hit;
	if (!ok)
	{																										// if it does then
		look        = ray->Storage()->AverageLookVectorAwayFromObserver(ray->Storage()->NumCells()-1);		// Get the look vector of teh last segment of the ray (away from point)
		ok          = ray->Storage()->LocationOfPoint( ray->Storage()->NumCells(), &location);															// Get the ground point in geographic coordinates
		coszen      = -(location.LocalZenith() & look );										// Get the cosine of the zenith angle, -1 as ray is coming into ground point.
		R           = m_albedobrdf->ReflectanceFunction(coszen);								// Get the Reflectance function for this zenith angle
		if (isSingleScatterTerm)																// For Single scatter we also add in the cosine of zenith angle
		{																						// so
			cossza = location.CosSZA();															// get the cosine of the solar zenith angle
			if (cossza < 0) cossza = 0.0;  /* ndl303 added this check 2010-11-23 */				// if the sun is below the horizon at this point then set the value to zero
			R *= cossza;																		// Modify the reflectance by the incoming solar zenith angle factor.
		}

//		NXASSERT(( location.Altitude() == ray->LowestShellInRayTracingGrid() ));
		NXASSERT(( R >= 0.0 ));

		ok    =       tablesource->InterpolateTable( location, look, jdescriptor, N_ELEMENTS(jdescriptor), &nj, R);
		ok    = ok && m_jindexworkspace.InsertQuadraturePointEntries( jdescriptor, nj);		// Each ground point only uses 
	}																								// and we are done
	m_jindexworkspace.FinishCellEntries();
	ok = ok && jindex->DeepCopy( m_jindexworkspace );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_QuadratureSourceFunction_HomogenousCell_InternalDiffuse_V2::CellTransmissionFactor		2008-2-10*/
/** Calculates the exponent integral across a pathlength. The integral is 
 *	made in centimeters rather than meters.
 **/
/*---------------------------------------------------------------------------*/

double SKTRANSO_Quadrature_TLS_V2_Legacy::CellTransmissionFactor( const HELIODETIC_POINT& point, double celllengthmeters ) 
{
	double rho;
	double	factor;
	double	tau;
	double	celllengthcm;

	celllengthcm = celllengthmeters*100.0;
	rho          = m_opticalprops->TotalExtinctionPerCM( point );
	tau          = rho*celllengthcm;
	if (tau < 1.0E-10)
	{
		factor = celllengthcm;
	}
	else
	{
		factor  = (1.0 - exp( -tau ))/rho;
	}
	return factor;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_QuadratureSourceFunction_HomogenousCell_InternalDiffuse_V2::FillCellTransmissionBuffer		2008-5-6*/
/** This fills the m_cellfactors array with the cell factor multiplied by the
 *	optical depth of the shell from the observer.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::FillCellFactorBuffer(  const SKTRANSO_RayInternalOptical* ray, size_t* numweights, bool usecache )
{
	size_t								numcells;
	double								celllengthmeters;
	const SKTRANSO_RayInternalGeometry*	geomray;
	size_t								cellidx;
	HELIODETIC_POINT					point;
	bool								ok;

	numcells  = ray->GeometryRay()->Storage()->NumCells();
	if (!usecache)
	{
		geomray  = ray->GeometryRay();
		for (cellidx = 0; cellidx < numcells; cellidx++)
		{
			ok                            = geomray->Storage()->CellMidPoint( cellidx, &point );
			celllengthmeters              = geomray->Storage()->CellLength( cellidx);
			m_cellfactors[cellidx]        = CellTransmissionFactor( point, celllengthmeters)*m_shelltransmissionsbuffer[cellidx];
		}
	}
	*numweights = numcells;
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_Quadrature_SS_LookupTable_Uniform::QuadratureWeightsForAtmosphericSingleScatter		2010-3-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::QuadratureWeightsForAtmosphericSS( const SKTRANSO_RayInternalOptical* ray, double* weights, size_t* numweights)
{
	size_t									numcells;
	size_t									cellidx;
	const SKTRAN_RayStorage_Base*			storage;
	SKTRAN_PhaseMatrixScalar				scatteringcm2;
	double									cosscatangle;
	bool									ok = true;
	HELIODETIC_POINT						point;

	
	storage           = ray->GeometryRay()->Storage();
	numcells          = storage->NumCells();
	*numweights       = numcells;

	NXASSERT(( numcells <= m_maxshellsalongray ));

	for (cellidx = 0; cellidx < numcells; cellidx++)
	{
		ok               = storage->CellMidPoint( cellidx, &point );
		cosscatangle     = storage->AverageLookVectorAwayFromObserver(cellidx).Z();											// Get the cosine of the scattering angle (use From observer as Z is in opposite direction of vector from Sun)
		ok               = ok && m_opticalprops->GetScatteringCoefficientCM2( point, cosscatangle,  &scatteringcm2 );	// Get the scattering coefficient for this angle
		weights[cellidx] = m_cellfactors[cellidx]*scatteringcm2;
	}
	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::QuadratureWeightsForAtmosphericSS		 2015- 3- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::QuadratureWeightsForAtmosphericEmission( const SKTRANSO_RayInternalOptical* ray, double* weights, size_t* numweights)
{
	size_t									numcells;
	size_t									cellidx;
	const SKTRAN_RayStorage_Base*			storage;
	bool									ok = true;

	
	storage           = ray->GeometryRay()->Storage();
	numcells          = storage->NumCells();
	*numweights       = numcells;

	NXASSERT(( numcells <= m_maxshellsalongray ));

	for (cellidx = 0; cellidx < numcells; cellidx++)
	{
		weights[cellidx] = m_cellfactors[cellidx];
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_QuadratureOpticalDepth_V2::OpticalDepthOfCell		2008-1-31*/
/** Calculates the optical depth of a segment of a ray from startintercept
 *	to exitintercept where startinrcept is closer to the observer.  The
 *	algorithm uses linear interpolation of the extinction per unit length
 *	with altitude.
 *
 *	This code implicitly assumes that there are no rays that start on one
 *	shell, go through the tangent point and come out the other side. The user
 *	must ensure 
**/
/*---------------------------------------------------------------------------*/

double SKTRANSO_Quadrature_TLS_V2_Legacy::OpticalDepthOfSegment( size_t cellidx, const SKTRAN_RayStorage_Base* storage) 
{
	bool								ok;
	double								tau;
//	double								r0,t0;
//	double								h1,h0;
//	double								r1,t1;
//	double								rt;
//	double								sigmak;
//	double								sigmaf;
//	double								dl;
//	SKTRAN_Distance			centershell;
//	HELIODETIC_POINT					location;
	HELIODETIC_POINT								startpt;
	HELIODETIC_POINT								endpt;
	SKTRAN_OpticalDepthCalculator_LinearWithHeight	opticaldepthcalculator;
	double											sigma0;
	double											sigma1;


//	NXASSERT(( startintercept <= exitintercept ));

/*
	rt = ray->GetTangentRadius ( );
	NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_QuadratureOpticalDepth_V2::OpticalDepthOfSegment, Rewrite code so we dont call GetInterceptAltitudeTangentAndZenith twice as its a big inefficiency hit, possibly cache r0 and r1 info in ray or Jindex tables\n"));

	r0 = ray->ShellRadius( cellidx );
	r1 = ray->ShellRadius( cellidx + 1);

	centershell = 0.5*( startintercept + exitintercept );
	ray->LocationAlongRayAsPoint( centershell, &location  );			// Get the nominal heliodetic location of this ray segment.
*/
	storage->LocationOfPoint( cellidx,     &startpt );			// Get the nominal heliodetic location of this ray segment.
	storage->LocationOfPoint( cellidx+1 ,  &endpt   );			// Get the nominal heliodetic location of this ray segment.

	opticaldepthcalculator.ConfigureQuadratureCoefficients( startpt, endpt);
	sigma0 = m_opticalprops->TotalExtinctionPerCM(startpt);
	sigma1 = m_opticalprops->TotalExtinctionPerCM(endpt);
	tau     = opticaldepthcalculator.OpticalDepthFromStartToEnd( sigma0, sigma1);
		
/*
	ok =       ray->GetInterceptAltitudeTangentAndZenith( startintercept, &h0, &t0, r0);
	ok = ok && ray->GetInterceptAltitudeTangentAndZenith( exitintercept,  &h1, &t1, r1);
	ok = ok && m_opticalprops->GetEffectiveExtinctionPerCMWithHeight( location, h0,r0, h1,r1, &sigmak, &sigmaf );
    sigmak *= 100.0;															// Convert extictions per CM to extinction per meter
	sigmaf *= 100.0;
	if (ok)
	{
		dl      = fabs( (double)(exitintercept - startintercept) );				// Watch out for path lengths that are very small
		if (dl < 0.0001)														// They ar eusually due to rounding issues
		{																		// if we have one
			tau = 0.0;															// then set the optical depth to zero
		}																		// and we are done.
		else
		{
			if ( t1 >= t0 )										// Case 1, Upward rays, or very short cell segments at the tangent point
			{
				tau = sigmak*(t1-t0) + sigmaf*0.5*( (r1*t1-r0*t0) + rt*rt*log( (r1+t1)/(r0+t0) ) );
			}
			else 												// case 2. downward ray, passes from upper shell to lower shell
			{
				tau = sigmak*(t0-t1) + sigmaf*0.5*( (r0*t0-r1*t1) + rt*rt*log( (r0+t0)/(r1+t1) ) );
			}
//	#if defined(_DEBUG)
			/*
			* This is included to check that rays that start and finish at the same height are really very short path length
			* objects defined near to the tangent point. They are often the result of floating point limitations.
			* Paths that start on one side of the tangent point and finish on the other need to break the calculation
			* down into two sections or use the formula
			*
			* 2011-05-31 ndl303. I have had problems with this check for the internal solar transmission table used for each ray.
			* In this table there are occassional ponts
			*/
/*			double				cossza;
			HELIODETIC_VECTOR	pt;

			if (h1 == h0 )
			{
				HELIODETIC_UNITVECTOR look ;

				ray->LocationAlongRayAsVector( startintercept, &pt);
				look = ray->AverageLookVector();
				cossza = (pt.X()*look.X() + pt.Y()*look.Y() +pt.Z()*look.Z())/pt.Magnitude();
				ok     =  ray->GetInterceptAltitudeTangentAndZenith( startintercept, &h0, &t0, r0 );
				NXASSERT(( (cossza > -0.00001) && (cossza < 0.00001) ));
	//			NXASSERT(( fabs( exitintercept - startintercept) < 3.0000 ));
			}

//	#endif
		}
	}
#if defined(NXDEBUG)
	if (tau > 0.0)
	{
		double frac = fabs( (tau-od)/tau );
		NXASSERT(( frac < 0.01 ));
	}
	else
	{
		NXASSERT(( od == 0.0));
	}
#endif
*/
	ok = (tau >= 0.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_QuadratureOpticalDepth_V2::OpticalDepthOfSegment, Error looking up optical depth of a segment");
		tau = 0;
	}
	return tau;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_QuadratureOpticalDepth_V2::OpticalDepthOfCell		2008-1-31*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRANSO_Quadrature_TLS_V2_Legacy::OpticalDepthOfCell(  const SKTRAN_RayStorage_Base* storage, size_t cellidx )
{
	return OpticalDepthOfSegment( cellidx, storage );
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CalculateRayTransmission		2008-2-8*/
/** Calculates the total transmission of a ray and, if required, the
 *	transmission from each shell along the ray to the ray origin.
 *
 * \par Theory
 *	This implementation assumes that the extinction, \f$\sigma\f$, varies linearly with height.
 *	Note that this is not the same as varying linearly along the ray. The optical depth between
 *	one end of a ray segment at radius \f$r_0\f$ and another segment at \f$r_1\f$
 *  is given by,
 *
 *	\f$ \Delta\tau = \sigma_K\left( t_1 - t_0 \right) + \frac{\sigma_F}{2}\left( (r_1t_1 - r_0t_0) + r_T^2\ln\left[\frac{r_1+t_1}{r_0+t_0}\right]\right) \f$
 *
 *	where
 *		- \f$ \Delta\tau\f$  is the optical depth of ray travelling between a shell at radius \f$r_1\f$ and another shell at radius \f$r_0\f$.
 *		- \f$t_1\f$ is the radius of one end point o fthe ray.
 *	    - \f$t_0\f$ is the radius o fthe other end point of the ray.
 *		- \f$t_1\f$ is the distance of the end point from the rays tangent point.
 *	    - \f$t_0\f$ is the distance of the other end point from the ray's tangent point ,
 *		- \f$\sigma(r)\f$ is the extinction at radius \f$r\f$ derived from a linear interpolation of the tabulated extinction: \f$\sigma_K + \sigma_{F}\,r\f$ 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CalculateRayTransmission( SKTRANSO_RayInternalOptical* ray, double* usertransmission, bool totaltransmissiononly, bool usecachedtransmission )
{
	size_t								numcells;
	size_t								cellidx;
	double								tau;
	double								dtau;
	const SKTRAN_RayStorage_Base*		geometryray;
	double								transmission;
	bool								ok;

	transmission = 1.0;
	tau          = 0.0;
	geometryray = ray->GeometryRay()->Storage();
	numcells    = geometryray->NumCells();

	if ( totaltransmissiononly )																			// The rays in the solar transmission table only need
	{																										// the total transmission
		for (cellidx = 0; cellidx < numcells; cellidx++ )													// so we can do that a little more efficiently
		{																									// for each cell
			dtau   = OpticalDepthOfCell( geometryray, cellidx );								// get the optical depth of the cell
			tau   -= dtau;																					// get the total optical depth
			NXASSERT(( dtau >= 0.0 ));																		// make sure everything is going as planned
		}																									// do all the cells
		transmission = exp(tau);																			// and get the total transmission
	}
	else																									// Other rays also need to store the optical depth of each segment boundary
	{																										// from the observer, this is the purpose of the m_shelltransmissionbuffer
		if ( !usecachedtransmission)																		// if we are not useing a pre-calculated transmission
		{																									// then 
			for (cellidx = 0; cellidx < numcells; cellidx++ )												// we calculate the table for the entire ray
			{																									// so for each dcell
				m_shelltransmissionsbuffer[cellidx] = exp(tau);													// get the updated transmission of the "closest" point from the observer
				dtau                               = OpticalDepthOfCell( geometryray, cellidx );				// Get the optical depth of the cell
				tau                               -= dtau;														// increment the optical depth
				NXASSERT(( dtau >= 0.0 ));																		// and make sure weird stuff is not happening
			}																								// do all of the cells
			m_shelltransmissionsbuffer[numcells]  = exp(tau);												// and get the total transmission of all of the cells
		}																									// nothing to do if using cached result
		transmission =  m_shelltransmissionsbuffer[numcells];												// 
	}

	ok = ( tau  <= 0 );
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING,"SKTRAN_QuadratureOpticalDepth_V2::CalculateRayTransmission, This ray has a negative optical depth (tau = %e). That cannot be good!", (double)-tau );
	}
	*usertransmission = transmission;
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::RaySolarTransmissionTable		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRANSO_JindexTableBase*	SKTRANSO_Quadrature_TLS_V2_Legacy::RaySolarTransmissionTable( const SKTRANSO_RayInternalGeometry* ray)
{
	const SKTRANSO_JindexTableBase*	solartransmissiontable;

	solartransmissiontable = ray->InternalSolarTransmissionTable();														// Ask the ray what solar transmission table it wants to use
	if (solartransmissiontable == NULL) solartransmissiontable = m_enginetables->SolarTransmissionTable();		// If it returns NULL then use the default
	NXASSERT(( solartransmissiontable != NULL ));
	return solartransmissiontable;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_AtmosphericSingleScatter		2010-3-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_AtmosphericSingleScatter( const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex )
{
	return CreateJIndexTableForQuadraturePointsAlongRay( RaySolarTransmissionTable(ray), ray, jindex);
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_Emissions		 2015- 3- 4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_AtmosphericEmissions ( const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex)
{
	bool	ok;

	ok = m_enginetables->EmissionsTable()->IsEmpty();
	if (ok) jindex->Clear();
	else
	{
		ok = CreateJIndexTableForQuadraturePointsAlongRay( m_enginetables->EmissionsTable(), ray, jindex );
	}
	return ok; 
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_AtmosphericSingleScatter		2010-3-9*/
/** Calculate the single scatter signal along a ray using the Solar Transmission
 *	Table. The JIndex table is configured so that each point is the incoming
 *	solar irradiance at each quadrature point along the ray. linear combination o fthe 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_AtmosphericSingleScatter(	SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																					const SKTRANSO_JIndexArray&			jindex,
																					SKTRAN_JValueTable_V21*					jvalue,
																					bool									usecachedtransmission,
																					bool									usecachedcellfactors )
{
	bool	ok;
	size_t	numweights;
	double	transmission;

	ok = jindex.IsEmpty();
	if (ok) jvalue->Clear();
	else
	{
		ok =       jvalue->AttachToGeometry              ( jindex );																// Attach the Jindex to a Jvalue table
		ok = ok && jvalue->SetWeightsAndRadiancePtrs     ( RaySolarTransmissionTable(ray->InternalDiffuseGeometry()), SKTRAN_JSINGLESCATTER );	// Convert the indices into pointers
		ok = ok && CalculateRayTransmission              ( ray, &transmission, false, usecachedtransmission );						// Calculate the ray transmisison from each shell to the observer
		ok = ok && FillCellFactorBuffer					 ( ray, &numweights, usecachedcellfactors );
		ok = ok && QuadratureWeightsForAtmosphericSS     ( ray, m_scratchbuffer, &numweights);						// Get the weights for each quadrature point
		ok = ok && jvalue->AdjustWeightsAndTrim          ( m_scratchbuffer,  numweights );							// Apply the weights to the jvalue
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_Emissions		 2015- 3- 4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_AtmosphericEmissions( SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																				const SKTRANSO_JIndexArray&				jindex,
																				SKTRAN_JValueTable_V21*					jvalue,
																				bool									usecachedtransmission,
																				bool									usecachedcellfactors)
{
	bool	ok;
	size_t	numweights;
	double	transmission;

	ok = jindex.IsEmpty();
	if (ok) jvalue->Clear();
	else
	{
		ok =       jvalue->AttachToGeometry					( jindex );												// Attach the Jindex to a Jvalue table
		ok = ok && jvalue->SetWeightsAndRadiancePtrs		(m_enginetables->EmissionsTable(), SKTRAN_JEMISSION );		// Convert the indices into pointers
		ok = ok && CalculateRayTransmission					( ray, &transmission, false, usecachedtransmission );		// Calculate the ray transmisison from each shell to the observer
		ok = ok && FillCellFactorBuffer						( ray, &numweights, usecachedcellfactors );
		ok = ok && QuadratureWeightsForAtmosphericEmission  ( ray, m_scratchbuffer, &numweights);						// Get the weights for each quadrature point
		ok = ok && jvalue->AdjustWeightsAndTrim				( m_scratchbuffer,  numweights );							// Apply the weights to the jvalue
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_AtmosphericDiffuseScatter		2010-3-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_AtmosphericDiffuseScatter( const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex )
{
	return CreateJIndexTableForQuadraturePointsAlongRay( m_enginetables->DiffusePointsTable(), ray, jindex );
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_AtmosphericDiffuseScatter		2010-3-10*/
/** Calculate the weights to be applied to the source function at each
 *	quadradture point along the ray. The weight for the homoegenous cell is given by
 *	$f$w_i= e^{-\tau(s_0)}\left[ \frac{ 1 - e^{-\bar{\sigma}L}}{\bar{\sigma}}\right]$f$ 
 *
 *	\param shelltransmissions
 *	An input array [numshells] of transmissions from each of the shell boundaries to the ray orgin/observer.
 *
 *	\param weights
 *	An output array [numquadpoints] of weights to be applied to each quadrature point along the
 *	ray.
 *	
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_AtmosphericDiffuseScatter(  SKTRANSO_RayInternalOptical*		ray,
																					const SKTRANSO_JIndexArray&	jindex,
																					SKTRAN_JValueTable_V21*			jvalue,
																					bool							usecachedtransmission,
																					bool						    usecachedcell_cellfactors
																				  )
{
	bool	ok;
	size_t	numweights;
	double	transmission;

	ok =       jvalue->AttachToGeometry				( jindex );													// Attach the Jindex to a Jvalue table
	ok = ok && jvalue->SetWeightsAndRadiancePtrs	( m_enginetables->DiffusePointsTable(), SKTRAN_JDIFFUSE );		// Convert the indices into pointers using the Diffuse Points Table
	ok = ok && CalculateRayTransmission             ( ray, &transmission, false, usecachedtransmission );			// Calculate the ray transmisison from each shell to the observer
	ok = ok && FillCellFactorBuffer					( ray, &numweights, usecachedcell_cellfactors );				// Get the factors to be applied to each source function at each quadrature point
	ok = ok && jvalue->AdjustWeightsAndTrim			( m_cellfactors,  numweights );									// Apply the weights to the jvalue
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_GroundDiffuseScatter		2010-3-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_InterpolateGroundDiffuseScatter( const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex )
{
	return CreateJIndexTableForGroundPoint( m_enginetables->DiffuseGroundPointTable(), ray, jindex, false );
}




/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::QuadratureWeightsForMSGroundPoint		2010-3-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::QuadratureWeightsForMSGroundPoint( const SKTRANSO_RayInternalOptical*			ray,
																					  double*									weights,
																					  size_t*									numweights )
{
	if( ray->GeometryRay()->Storage()->GroundIsHit())														// it only makes sense if ground is hit;
	{																							// if it does then
		*numweights = 1;
		*weights    = ray->TotalTransmissionAlongRay();
	}
	else
	{
		*numweights = 0;
	}
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_GroundDiffuseScatter		2010-3-11*/
/** Creates the JValue table that is used to evaluate the contribution of
 *	the reflection of the diffuse field off the ground into the end point of the ray.
 *	The ground point table stores the "reflected" upward flux. The Jindex has the 
 *	albedo angular component already built in so all this code has to do is propagate
 *	the upward flux to the end point of the ray.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_InterpolateGroundDiffuseScatter(         SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																						   const SKTRANSO_JIndexArray&					jindex,
																							     SKTRAN_JValueTable_V21*					jvalue,
																							     bool									usecachedtransmission
																						  )
{
	bool		ok;
	double		transmission;
	size_t		numweights;

	ok =       jvalue->AttachToGeometry				( jindex );															// Attach the Jindex to a Jvalue table
	ok = ok && jvalue->SetWeightsAndRadiancePtrs	( m_enginetables->DiffuseGroundPointTable(),SKTRAN_JGROUND_UPFLUX );	// Convert the indices into pointers using the Diffuse Points Table
	ok = ok && CalculateRayTransmission             ( ray, &transmission, false, usecachedtransmission );					// Calculate the ray transmisison from each shell to the observer
	ok = ok && QuadratureWeightsForMSGroundPoint	( ray, m_scratchbuffer, &numweights);						// Get the factors to be applied to each multiple scatter source function
	ok = ok && jvalue->AdjustWeightsAndTrim			( m_scratchbuffer,  numweights );										// Apply the weights to the jvalue
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_GroundSingleScatter		2010-3-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_GroundSingleScatter( const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex )
{
	return CreateJIndexTableForGroundPoint( RaySolarTransmissionTable(ray), ray, jindex, true );
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_GroundEmissions		 2015- 3- 5*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_GroundEmissions( const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex)
{
	bool								ok;
	HELIODETIC_UNITVECTOR				look;
	SKTRANSO_JIndex						jdescriptor[100];
	size_t								nj;
	HELIODETIC_POINT					location;
	

	ok = (!ray->Storage()->GroundIsHit()) || ( m_enginetables->EmissionsTable()->IsEmpty());
	if (ok)
	{
		jindex->Clear();
	}
	else
	{
																										// if it does then
		double	nan =  std::numeric_limits<double>::quiet_NaN();

		m_jindexworkspace.ResetCounters();																	// Reset the counters in this temporary storage area
		look.SetCoords( nan, nan, nan);											
		ok    =       m_enginetables->EmissionsTable()->InterpolateTable( location, look, jdescriptor, N_ELEMENTS(jdescriptor), &nj, 1.0);
		ok    = ok && m_jindexworkspace.InsertQuadraturePointEntries( jdescriptor, nj);	
		m_jindexworkspace.FinishCellEntries();
		ok = ok && jindex->DeepCopy( m_jindexworkspace );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_GroundEmissions		 2015- 3- 5*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_GroundEmissions(		 SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																		   const SKTRANSO_JIndexArray&					jindex,
																				 SKTRAN_JValueTable_V21*				jvalue,
																				 bool									usecachedtransmission)
{
	bool		ok;
	double		transmission;

	ok = (jindex.IsEmpty()) ||  (m_enginetables->EmissionsTable()->IsEmpty());
	if (ok)
	{
		jvalue->Clear();
	}
	else
	{
		ok =       jvalue->AttachToGeometry				( jindex );																// Attach the Jindex to a Jvalue table
		ok = ok && jvalue->SetWeightsAndRadiancePtrs	(  m_enginetables->EmissionsTable(),SKTRAN_JEMISSION );					// Convert the indices into pointers using the Solar Transmission Table
		ok = ok && CalculateRayTransmission             ( ray, &transmission, false, usecachedtransmission );					// Calculate the ray transmisison from each shell to the observer
		ok = ok && jvalue->AdjustWeightsByConstantFactorAndTrim( transmission);					// Get the factors to be applied to each multiple scatter source function
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::QuadratureWeightsForSSGroundPoint		2009-1-20*/
/** Create the JValue table that will evaluate the contribution of single
 *	scatter from the ground to the observer at the ray origin. The Jindex table
 *	stores the cos(sza) which converts the downward irradiance to a downward flux
 *	 and it also stores the angular dependence of the albedo. This code must apply the
 *	wavelength dependent albedo and must propagate the upward flux to the observer
 *	at the origin.
 *
 *	This code finds the downard irradiance by interpolating the points in the
 *	Solar Transmission Table. This this quadrature class is not ideal for Single Scatter
 *	calculations as it creates the entire solar transmission table.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::QuadratureWeightsForSSGroundPoint( const SKTRANSO_RayInternalOptical*			ray,
																					  const SKTRAN_TableOpticalProperties_V21*	opticalprops,
																					  double*									weights,
																					  size_t*									numweights)
{
	const SKTRAN_RayStorage_Base*			storage;
	double									albedo;
	HELIODETIC_POINT						groundpoint;
	bool									ok = true;

	
	storage   = ray->GeometryRay()->Storage();
	ok   = !storage->GroundIsHit();
																								// it only makes sense if ground is hit;
	if (!ok)
	{																							// if it does then
		storage->LocationOfPoint( storage->NumCells(), &groundpoint );
		ok           = opticalprops->Get_AlbedoForDeprecatedLegacyCode( groundpoint, &albedo );				// Get the optical properties at this point
		*weights     = albedo*ray->TotalTransmissionAlongRay();
		*numweights  = 1;
	}
	else
	{
		*weights = 0.0;
		*numweights = 0;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_GroundSingleScatter		2010-3-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_GroundSingleScatter(  SKTRAN_RayInternalDiffuseOptical_V2*			ray,
																			  const SKTRANSO_JIndexArray&		jindex,
																			  SKTRAN_JValueTable_V21*			jvalue,
																			  bool								usecachedtransmission
																				     	)
{
	bool		ok;
	double		transmission;
	size_t		numweights;

	ok =       jvalue->AttachToGeometry				( jindex );																// Attach the Jindex to a Jvalue table
	ok = ok && jvalue->SetWeightsAndRadiancePtrs	( RaySolarTransmissionTable(ray->InternalDiffuseGeometry()),SKTRAN_JGROUNDSS );		// Convert the indices into pointers using the Solar Transmission Table
	ok = ok && CalculateRayTransmission             ( ray, &transmission, false, usecachedtransmission );					// Calculate the ray transmisison from each shell to the observer
	ok = ok && QuadratureWeightsForSSGroundPoint	( ray, m_opticalprops, m_scratchbuffer, &numweights);					// Get the factors to be applied to each multiple scatter source function
	ok = ok && jvalue->AdjustWeightsAndTrim			( m_scratchbuffer,  numweights );										// Apply the weights to the jvalue
	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_GroundPointIncomingRaysUpwardFlux		2010-5-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJIndexTable_GroundPointIncomingRaysUpwardFlux	( SKTRAN_GroundPointDiffuseGeometry_V21* groundpoint )
{
	bool	ok;

	ok = groundpoint->CreateIncomingJIndexTableForDownwardFlux( );
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_GroundPointIncomingRaysUpwardFlux		2010-5-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Quadrature_TLS_V2_Legacy::CreateJValueTable_GroundPointIncomingRaysUpwardFlux	( SKTRANSO_GroundPointDiffuseOptical*	groundpoint, const SKTRAN_TableOpticalProperties_V21* optprop )
{
	bool	ok;

	ok = groundpoint->ConfigureOptical( optprop, this);
	return ok;
}


