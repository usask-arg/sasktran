#include "../sktran_common.h"
#include <omp.h>

class SKTRAN_MCPhoton_Base;

/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::SKTRAN_SolarTransmission_NoTable		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SolarTransmission_NoTable::SKTRAN_SolarTransmission_NoTable()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::~SKTRAN_SolarTransmission_NoTable		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SolarTransmission_NoTable::~SKTRAN_SolarTransmission_NoTable()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::ReleaseResources		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_SolarTransmission_NoTable::ReleaseResources()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::MakeThreadSafeFor		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable::MakeThreadSafeFor( size_t numthreads )
{
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::TransmissionAtVector		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable::TransmissionAtVector	(	const HELIODETIC_VECTOR&						vec,
																double&											transmission) const
{
	bool ok = true;

	SKTRAN_RayOptical_Base*							rayopt;
	HELIODETIC_UNITVECTOR							sununit;
	std::unique_ptr<SKTRAN_RayOptical_Base>			ray;
	
	sununit.SetCoords(0,0,1);
	ok = ok && RayFactory()->CreateRayObject( &ray);
	ok = ok && ray->MoveObserver( vec, sununit);
	ok = ok && ray->TraceRay_NewMethod();
	if(ok)
	{
		rayopt = ray.get();
		if(rayopt->Storage()->GroundIsHit())
		{
			transmission = 0.0;
		} else
		{
			//ok = ok && m_rayintegrator->CalculateRayTransmission_withMinContainer(rayopt,&transmission,true,true);
			ok = ok && Integrator()->CalculateRayScalarTransmission(rayopt,&transmission,true,true);
			transmission = exp(-1.0*rayopt->TotalOpticalDepth());
		}
	}
	if (!ok) transmission = std::numeric_limits<double>::quiet_NaN();
	return ok;
}

bool SKTRAN_SolarTransmission_NoTable::TransmissionAtVector(const double&									wavelength,
															const HELIODETIC_VECTOR&						vec,
															double&											transmission) const
{
	bool ok = true;

	SKTRAN_RayOptical_Base*							rayopt;
	HELIODETIC_UNITVECTOR							sununit;
	std::unique_ptr<SKTRAN_RayOptical_Base>			ray;

	sununit.SetCoords(0, 0, 1);
	ok = ok && RayFactory()->CreateRayObject(&ray);
	ray->SetWavelength(wavelength);
	ok = ok && ray->MoveObserver(vec, sununit);
	ok = ok && ray->TraceRay_NewMethod();
	if (ok)
	{
		rayopt = ray.get();
		if (rayopt->Storage()->GroundIsHit())
		{
			transmission = 0.0;
		}
		else
		{
			//ok = ok && m_rayintegrator->CalculateRayTransmission_withMinContainer(rayopt,&transmission,true,true);
			ok = ok && Integrator()->CalculateRayScalarTransmission(rayopt, &transmission, true, true);
			transmission = exp(-1.0*rayopt->TotalOpticalDepth());
		}
	}
	if (!ok) transmission = std::numeric_limits<double>::quiet_NaN();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::FillTable_ClassSpecific		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable::FillTable_ClassSpecific( )
{
	// There is no actual table for this class
	return true;
}
 


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::TransmissionAtPoint		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable::TransmissionAtPoint	(	const HELIODETIC_POINT&							point, 
																double&											transmission) const
{	
	return TransmissionAtVector(point.Vector(),transmission);
}

bool SKTRAN_SolarTransmission_NoTable::TransmissionAtPoint  (   const double&									wavelength, 
																const HELIODETIC_POINT&							point,
																double&											transmission) const
{
	return TransmissionAtVector(wavelength, point.Vector(), transmission);
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::SourceTermAtPoint		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable::SourceTermAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source ) const
{
	bool ok = true;

	double transmission;
	double kscatt;

	ok = ok && TransmissionAtPoint( qobj.GetPoint(), transmission );
	ok = ok && Integrator()->GetOpticalProps()->GetScatteringCoefficientCM2( qobj.GetPoint(), CosAngleToSource(qobj.GetLookAway(), &qobj.GetPoint()), &kscatt );
	kscatt *= 100;
	NXASSERT((NXFINITE(kscatt)));
	*source = transmission * kscatt;

	return ok;
}

bool SKTRAN_SolarTransmission_NoTable::SourceTermAtPoint(const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, double* source) const
{
	bool ok = true;

	double transmission;
	double kscatt;

	ok = ok && TransmissionAtPoint(wavelength, qobj.GetPoint(), transmission);
	ok = ok && Integrator()->GetOpticalProps()->GetScatteringCoefficientCM2(wavelength, qobj.GetPoint(), CosAngleToSource(qobj.GetLookAway(), &qobj.GetPoint()), &kscatt);
	kscatt *= 100;
	NXASSERT((NXFINITE(kscatt)));
	*source = transmission * kscatt;

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::GroundSourceAtPoint		 2016- 12- 22*/
/** **/
/*---------------------------------------------------------------------------*/

//bool SKTRAN_SolarTransmission_NoTable::GroundSourceAtPoint( const SKTRAN_SourceTermQueryObject_Base& outboundray,
//													        double*									 source ) const
//{
//	bool ok = true;
//	double brdf;
//	double transmission;
//	double mu_in;
//	double	mu_out;
//	double cosdphi;
//	const HELIODETIC_POINT&        point = outboundray.GetPoint();
//	const HELIODETIC_UNITVECTOR&   lookv = outboundray.GetLookAway();
//	HELIODETIC_UNITVECTOR		   sunvector;
//
//	Sun()->UpdateSun();
//	Sun()->SunUnitVector( &sunvector);
//	ok = ok && TransmissionAtPoint( point, transmission);
//	SKTRAN_TableOpticalProperties_Base::GroundBRDFAngles( point, sunvector, lookv, &mu_in, &mu_out, &cosdphi);
//
//	ok = ok && Integrator()->GetOpticalProps()->GetBRDF( point, mu_in, mu_out, cosdphi, &brdf);
//
//	*source = transmission*mu_in*brdf;
//
//	return ok;
//}
//

/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::MonteCarlo_SingleScatteredRadianceAtPoint		 2016- 12- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable::MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj,
                                                                         double&                              radiance ) const
{
	bool ok = true;
	double transmission;
	double scattCoeff_solar;
	double kscatt;				
	size_t threadid = omp_get_thread_num();

	Sun()->UpdateSun();
	ok = ok && TransmissionAtPoint( qobj.GetPoint(), transmission );
	ok = ok && Integrator()->GetOpticalProps()->GetScatteringCoefficientCM2( qobj.GetPoint(), Sun()->CosAngleToSun(qobj.GetLookAway()), &scattCoeff_solar);	// New look is heliodetic --> cosScatAngle = dot(-zhat,newLook) = -newLook.z
	kscatt = Integrator()->GetOpticalProps()->ScatteringExtinctionPerCM(qobj.GetPoint());
	scattCoeff_solar = 0.0==kscatt ? 0.0 : scattCoeff_solar / kscatt;			// GetScatteringCoefficientCM2 returns P_normalized*kscatt

	radiance = transmission * scattCoeff_solar;

	return ok;
 }


bool SKTRAN_SolarTransmission_NoTable::MonteCarlo_SingleScatteredRadianceAtPoint( const double& wavelength, 
																				   const SKTRAN_SourceTermQueryObject_Base& qobj,
																				   double&                              radiance) const
{
	bool ok = true;
	double transmission;
	double scattCoeff_solar;
	double kscatt;
	size_t threadid = omp_get_thread_num();

	Sun()->UpdateSun();
	ok = ok && TransmissionAtPoint(wavelength, qobj.GetPoint(), transmission);
	ok = ok && Integrator()->GetOpticalProps()->GetScatteringCoefficientCM2(wavelength, qobj.GetPoint(), Sun()->CosAngleToSun(qobj.GetLookAway()), &scattCoeff_solar);	// New look is heliodetic --> cosScatAngle = dot(-zhat,newLook) = -newLook.z
	kscatt = Integrator()->GetOpticalProps()->ScatteringExtinctionPerCM(wavelength, qobj.GetPoint());
	scattCoeff_solar = 0.0 == kscatt ? 0.0 : scattCoeff_solar / kscatt;			// GetScatteringCoefficientCM2 returns P_normalized*kscatt
	radiance = transmission * scattCoeff_solar;
	return ok;
 }

/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable::MonteCarlo_GroundScatteredRadianceAtPoint		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable::MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, 
                                                                        double&                              radiance ) const
{
	bool ok = true;

	double transmission;
	double scattCoeff_solar;

	Sun()->UpdateSun();
	ok = ok && TransmissionAtPoint( qobj.GetPoint(), transmission);
	scattCoeff_solar = Sun()->CosAngleToSun(qobj.GetPoint().UnitVector())/nxmath::Pi;					// Lambertian scatter --> intensity of light at lambertian scatter point ~ cos(angleBetweenSurfaceNormalAndSolarDirection)  ~=  scatterVector_{heliodetic}.dot(solarVector_{heliodetic})  =  scatterVector_{heliodetic}.dot([0 0 1]) = scatterVector_{heliodetic}.Z(); factor \pi comes from normalization over rotation about \hat{r}; This should really be done inside a Ground_Lambertian class
	if(scattCoeff_solar < 0.0) scattCoeff_solar = 0.0;

	radiance = transmission * scattCoeff_solar;

	return ok;
}

 
bool SKTRAN_SolarTransmission_NoTable::MonteCarlo_GroundScatteredRadianceAtPoint(  const double&							 wavelength, 
																				const SKTRAN_SourceTermQueryObject_Base& qobj,
																					double&									 radiance) const
{
	bool ok = true;

	double transmission;
	double scattCoeff_solar;

	Sun()->UpdateSun();
	ok = ok && TransmissionAtPoint(wavelength, qobj.GetPoint(), transmission);
	scattCoeff_solar = Sun()->CosAngleToSun(qobj.GetPoint().UnitVector()) / nxmath::Pi;					// Lambertian scatter --> intensity of light at lambertian scatter point ~ cos(angleBetweenSurfaceNormalAndSolarDirection)  ~=  scatterVector_{heliodetic}.dot(solarVector_{heliodetic})  =  scatterVector_{heliodetic}.dot([0 0 1]) = scatterVector_{heliodetic}.Z(); factor \pi comes from normalization over rotation about \hat{r}; This should really be done inside a Ground_Lambertian class
	if (scattCoeff_solar < 0.0) scattCoeff_solar = 0.0;

	radiance = transmission * scattCoeff_solar;

	return ok;
}

/**************************************************
 * SKTRAN_SolarTransmission_NoTable_reuseRays code
 **************************************************/

SKTRAN_SolarTransmission_NoTable_reuseRays::SKTRAN_SolarTransmission_NoTable_reuseRays()
{
	m_rayopt1.SetThreadStorageInitializer( std::bind( &SKTRAN_SolarTransmission_NoTable_reuseRays::InitializeThreadSafeStorageEntry, this, std::placeholders::_1) );   //Configure the function used to initialize thread storage data. 
	// *** NOTE _1 MAY NOT WORK AS IS IT REFERS TO std::placeholders::_1 and WE MUST USE std::placeholders::_1				  //Configure the function used to initialize thread storage data. 
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable_reuseRays::~SKTRAN_SolarTransmission_NoTable_reuseRays		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SolarTransmission_NoTable_reuseRays::~SKTRAN_SolarTransmission_NoTable_reuseRays(){
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable_reuseRays::InitializeThreadSafeStorageEntry		 2014- 11- 12*/
/**	The function used to initialize newly created thread storage data. This
 *	function must be thread safe as it will often be running in any one of
 *	a number of worker threads. This code simply uses the RayFactory to create a new 
 *	ray object on the heap with auto destruct capability.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable_reuseRays::InitializeThreadSafeStorageEntry( std::unique_ptr<SKTRAN_RayOptical_Base>* ptr ) const
{
		return RayFactory()->CreateRayObject( ptr );
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable_reuseRays::Allocate		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable_reuseRays::Allocate(size_t numThreads){

	bool ok = true;
	
//	m_sununit.SetCoords(0.0,0.0,1.0);
	DeallocateInternalRays( );
//	m_rayopt.resize(numThreads);
//	for(size_t ridx=0; ridx < numThreads; ridx++)
//	{
//		ok = ok && RayFactory()->CreateRayObject( &m_rayopt[ridx] );
//	}
//
//	if(!ok) DeallocateInternalRays( );


	return ok;	
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable_reuseRays::ReleaseResources		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_SolarTransmission_NoTable_reuseRays::ReleaseResources(){

	DeallocateInternalRays( );
	SKTRAN_SolarTransmission_NoTable::ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable_reuseRays::DeallocateInternalRays		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_SolarTransmission_NoTable_reuseRays::DeallocateInternalRays( )
{
	m_rayopt1.Clear();
	m_numThreadsSafeFor = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable_reuseRays::MakeThreadSafeFor		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable_reuseRays::MakeThreadSafeFor(size_t numThreads){

	bool ok;

	ok = Allocate(numThreads);
	m_numThreadsSafeFor = ok ? numThreads : 0;

	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable_reuseRays::TransmissionAtVector		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable_reuseRays::TransmissionAtVector	(	const HELIODETIC_VECTOR&					vec,
																			double&										transmission) const
{
	std::unique_ptr<SKTRAN_RayOptical_Base>*	raytls;

	bool ok = true;

//	m_rayopt[threadid]->Initialize(coords, vec, m_sun->GetSunUnit(threadid));
//	m_rayopt[threadid]->StorageAccessVar()->SetCoords( &coords );
//  ok = ok && RayFactory()->CreateRay( &(m_rayopt[threadid]));
	ok =       m_rayopt1.LookupUpThreadData( &raytls );
	ok = ok && (*raytls)->MoveObserver(  vec,  Sun()->GetSunUnit() );
	ok = ok && (*raytls)->TraceRay_NewMethod();
	if(ok)
		if( (*raytls)->Storage()->GroundIsHit())
		{
			transmission = 0.0;
		}
		else
		{
			ok = Integrator()->CalculateRayScalarTransmission_withMinContainer(raytls->get(),&transmission,true,true);
			transmission = std::exp(- (*raytls)->TotalOpticalDepth());
		}
	return ok;
}

bool SKTRAN_SolarTransmission_NoTable_reuseRays::TransmissionAtVector(	const double&								wavelength,
																		const HELIODETIC_VECTOR&					vec,
																		double&										transmission) const
{
	std::unique_ptr<SKTRAN_RayOptical_Base>*	raytls;

	bool ok = true;

	//	m_rayopt[threadid]->Initialize(coords, vec, m_sun->GetSunUnit(threadid));
	//	m_rayopt[threadid]->StorageAccessVar()->SetCoords( &coords );
	//  ok = ok && RayFactory()->CreateRay( &(m_rayopt[threadid]));
	ok = m_rayopt1.LookupUpThreadData(&raytls);
	ok = ok && (*raytls)->MoveObserver(vec, Sun()->GetSunUnit());
	(*raytls)->SetWavelength(wavelength);
	ok = ok && (*raytls)->TraceRay_NewMethod();
	
	if (ok)
		if ((*raytls)->Storage()->GroundIsHit())
		{
			transmission = 0.0;
		}
		else
		{
			ok = Integrator()->CalculateRayScalarTransmission_withMinContainer(raytls->get(), &transmission, true, true);
			transmission = std::exp(-(*raytls)->TotalOpticalDepth());
		}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_NoTable_reuseRays::TransmissionAtPoint		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_NoTable_reuseRays::TransmissionAtPoint(	const HELIODETIC_POINT&	point, 
																		double&					transmission) const
{

	return TransmissionAtVector(point.Vector(), transmission);
}

bool SKTRAN_SolarTransmission_NoTable_reuseRays::TransmissionAtPoint(	const double& wavelength, 
																		const HELIODETIC_POINT&	point,
																		double&					transmission) const
{

	return TransmissionAtVector(wavelength, point.Vector(), transmission);
}

bool SKTRAN_SolarTransmission_NoTable_NoLOSSource::SourceTermAtPoint(const SKTRAN_SourceTermQueryObject_Base& qobj, double* source) const
{
	*source = 0.0;

	return true;
}