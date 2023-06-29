#include "include/sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Second_Order_Source::SKTRAN_HR_Diffuse_Second_Order_Source		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Diffuse_Second_Order_Source::SKTRAN_HR_Diffuse_Second_Order_Source()
{
	m_incomingsphere = NULL;
	m_integrator     = NULL;
	m_opticaltable   = NULL;
	m_rayfactory     = NULL;
	m_solartable     = NULL;
    m_emissiontable  = nullptr;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Second_Order_Source::~SKTRAN_HR_Diffuse_Second_Order_Source		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Diffuse_Second_Order_Source::~SKTRAN_HR_Diffuse_Second_Order_Source()
{
	if( NULL != m_incomingsphere ) m_incomingsphere->Release();
	if( NULL != m_integrator ) m_integrator->Release();
	if( NULL != m_opticaltable ) m_opticaltable->Release();
//	if( NULL != m_rayfactory ) m_rayfactory->Release();
	if( NULL != m_solartable ) m_solartable->Release();
    if( nullptr != m_emissiontable ) m_emissiontable->Release();

	m_incomingsphere = NULL;
	m_integrator = NULL;
	m_opticaltable = NULL;
	m_rayfactory   = NULL;
	m_solartable = NULL;
    m_emissiontable = nullptr;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Second_Order_Source::RotateIncomingRay		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Second_Order_Source::RotateIncomingRay( const nxVector& inray, const HELIODETIC_POINT& pt, HELIODETIC_UNITVECTOR& outlook ) const
{
	bool ok = true;

	HELIODETIC_UNITVECTOR				xprime,yprime,zprime;
	double								x,y,z;
	HELIODETIC_UNITVECTOR				unitvectors[3];

	pt.LocalUnitVectors( &unitvectors[0], 3 );
	xprime = unitvectors[0];
	yprime = unitvectors[1];
	zprime = unitvectors[2];

	x = inray.X();
	y = inray.Y();
	z = inray.Z();

	outlook.SetCoords( x*xprime.X() + y*yprime.X() + z*zprime.X(),				// look = x*xprime + y*yprime + z*zprime
					   x*xprime.Y() + y*yprime.Y() + z*zprime.Y(),
					   x*xprime.Z() + y*yprime.Z() + z*zprime.Z() );			// Get the look vector in heliodetic coordinates


	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Second_Order_Source::RadianceFromLook		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Second_Order_Source::RadianceFromLook( const HELIODETIC_UNITVECTOR& look, const HELIODETIC_POINT& pt, double& radiance ) const
{
	bool											ok;
	std::unique_ptr<SKTRAN_RayOptical_Base>			ray;
    double                                          temp = 0.0;

	ok =       RayFactoryPtr()->CreateRayObject( &ray );
	ok = ok && ray->MoveObserver( pt.Vector(), look);
	ok = ok && ray->TraceRay_NewMethod();
	ok = ok && m_integrator->CalculateRayScalarTransmission( ray.get(), NULL, false, false );

    // Should change this to loop over all sources... 
    ok = ok && m_integrator->SingleScatterRadiance   ( ray.get(), temp, *m_solartable );
    radiance = 0.0;
    radiance += temp;
    ok = ok && m_integrator->SingleScatterRadiance   ( ray.get(), temp, *m_emissiontable );
    radiance += temp;
     
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Diffuse_Second_Order_Source::RadianceFromLook, Error calculating single scatter radiance. Thats not good.");
		radiance = 0;
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Second_Order_Source::SourceTermAtPoint		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Second_Order_Source::SourceTermAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source ) const
{
	bool ok = true;
	
	HELIODETIC_UNITVECTOR outlook;
	double				  radiance;
	double				  scattcoeff;
	double				  cosangle;

	*source = 0.0;
	for( size_t rayidx = 0; rayidx < m_incomingsphere->NumUnitVectors(); rayidx++ )
	{
		ok = ok && RotateIncomingRay( m_incomingsphere->UnitVectorAt( rayidx ), qobj.GetPoint(), outlook );
		ok = ok && RadianceFromLook( outlook, qobj.GetPoint(), radiance );
		cosangle = -1.0*(qobj.GetLookAway().X() * outlook.X() + qobj.GetLookAway().Y() * outlook.Y() + qobj.GetLookAway().Z() * outlook.Z());
		ok = ok && m_opticaltable->GetScatteringCoefficientCM2( qobj.GetPoint(), cosangle, &scattcoeff );
		*source += radiance * scattcoeff * m_incomingsphere->CubatureWeightAt( rayidx );

	}
	*source *= 100;		// cm to m
	return ok;
}

bool SKTRAN_HR_Diffuse_Second_Order_Source::SourceTermAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source ) const
{
	bool ok = true;
	// No polarization for now
	source->SetTo(0.0);
	SKRTFLOAT temp_scalar;
	ok = ok && SourceTermAtPoint( qobj, &temp_scalar );
	source->Assign_I( temp_scalar );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Second_Order_Source::SetUnitSphere		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Second_Order_Source::SetUnitSphere( const SKTRAN_UnitSphere_V2& unitsphere )
{
	bool ok = true;

	m_incomingsphere = &unitsphere;
	m_incomingsphere->AddRef();

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Second_Order_Source::SetOptical		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Second_Order_Source::SetOptical( std::shared_ptr< const SKTRAN_RayFactory_Base>  rayfactory,
														const SKTRAN_OpticalPropertiesIntegrator_Base&  integrator,
														const SKTRAN_TableOpticalProperties_Base&       opticaltable,
														const SKTRAN_SolarTransmission_Base&            solartable, 
                                                        const SKTRAN_EmissionTable_Base&                emissiontable )
{
	bool ok = true;

//	if (rayfactory   != nullptr)   rayfactory->AddRef();
//	if (m_rayfactory != nullptr) m_rayfactory->Release();
	m_rayfactory = rayfactory;


	m_integrator = &integrator;
	m_integrator->AddRef();
	m_opticaltable = &opticaltable;
	m_opticaltable->AddRef();
	m_solartable = &solartable;
	m_solartable->AddRef();
    m_emissiontable = &emissiontable;
    m_emissiontable->AddRef();
    

	return ok;
}
