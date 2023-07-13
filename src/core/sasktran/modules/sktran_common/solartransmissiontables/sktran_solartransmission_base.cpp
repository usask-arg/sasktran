#include "../sktran_common.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base::SKTRAN_SolarTransmission_Base		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SolarTransmission_Base::SKTRAN_SolarTransmission_Base( )
{
	m_sun        = nullptr;
	m_integrator = nullptr;
//	m_rayfactory  = nullptr;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base::~SKTRAN_SolarTransmission_Base		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SolarTransmission_Base::~SKTRAN_SolarTransmission_Base( )
{
	ReleaseResources( );
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base::ReleaseResources		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_SolarTransmission_Base::ReleaseResources( )
{
	//m_sun        = nullptr; // This could be done, but only if ReleaseResources is only ever called by the destructor (otherwise it's confusing, the user may not know to call SetSun)
	if(nullptr!=m_integrator)  m_integrator->Release();  m_integrator = nullptr;
///	if(nullptr!=m_rayfactory ) m_rayfactory ->Release();  m_rayfactory  = nullptr;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base::SetSun		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_Base::SetSun( std::unique_ptr<SKTRAN_Sun_Base>& sun )
{
	bool ok = true;

	ok = ok && nullptr!=sun;	// Should have manager class so we can check if threadsafe
	if(ok)
	{
        m_sun = std::move(sun);
    }
    else {
        nxLog::Record(NXLOG_WARNING, "SKTRAN_SolarTransmission_Base::SetSun, Received null sun.");
    }

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base::FillTable		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_Base::FillTable( )
{
	bool ok = true;

	ok = ok && (nullptr!=m_sun);
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_SolarTransmission_Base::FillTable, Must have sun and ray factory defined before filling table!");
	ok = ok && FillTable_ClassSpecific( );

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base::ConfigureOptical		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_Base::ConfigureOptical(  std::weak_ptr< const SKTRAN_RayFactory_Base>	rayfactory, SKTRAN_OpticalPropertiesIntegrator_Base* integrator)
{
	bool ok = true;

//	if (   rayfactory != nullptr )  rayfactory->AddRef();
//	if ( m_rayfactory != nullptr )m_rayfactory->Release();
	if (   integrator != nullptr)   integrator->AddRef();
	if ( m_integrator != nullptr) m_integrator->Release();

	m_rayfactory = rayfactory;
	m_integrator = integrator;

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base::EstimatePolarizationVector		 2014- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_SolarTransmission_Base::EstimateNormalizedPolarizationVector( const HELIODETIC_POINT& pt, const HELIODETIC_UNITVECTOR& look, SKTRAN_Stokes_NC& vec, const HELIODETIC_BASIS& basis) const
{
	bool ok = true;
	double cosEta_scatToObs, sinEta_scatToObs;
	double cosToSource = CosAngleToSource( look, &pt );	// Get angle of this scatter

	ok = ok && m_integrator->GetOpticalProps()->GetResultOfUnpolarizedScatterCM2( pt, cosToSource, vec );		// Particle tracing goes through same angles as path tracing
	
	// Rotate the basis for the electric field components used in this vector through angle eta (Tony's thesis page 44), rotates polarization plane CW about look dir
	if( 1e-6 < 1.0 - cosToSource*cosToSource ){
		// Find angles to rotate scattering plane into observer plane
		double oneOverSinScatt = 1.0 / sqrt( 1.0 - cosToSource*cosToSource);
		cosEta_scatToObs =  CosAngleToSource( basis.y, &pt) * oneOverSinScatt; // cosEta =  dot(refDir,y) = dot(+sunDir,y)
		sinEta_scatToObs = -CosAngleToSource( basis.z, &pt) * oneOverSinScatt; // sinEta = -dot(refDir,z) =  -dot(+sunDir,z)
	

	} else{
		// Ignore rotation to avoid divide by zero
		cosEta_scatToObs = 1.0;
		sinEta_scatToObs = 0.0;
	}

	vec.RotatePolarPlaneThru( cosEta_scatToObs, sinEta_scatToObs ); // Result of pmatrix acting on source vector is in scatter frame -- rotate polarization into my frame
    vec.Normalize();

	return ok;
}

bool SKTRAN_SolarTransmission_Base::EstimateNormalizedPolarizationVector(const double& wavelength, const HELIODETIC_POINT& pt, const HELIODETIC_UNITVECTOR& look, SKTRAN_Stokes_NC& vec, const HELIODETIC_BASIS& basis) const
{
	bool ok = true;
	double cosEta_scatToObs, sinEta_scatToObs;
	double cosToSource = CosAngleToSource(look, &pt);	// Get angle of this scatter

	ok = ok && m_integrator->GetOpticalProps()->GetResultOfUnpolarizedScatterCM2(wavelength, pt, cosToSource, vec);		// Particle tracing goes through same angles as path tracing

	// Rotate the basis for the electric field components used in this vector through angle eta (Tony's thesis page 44), rotates polarization plane CW about look dir
	if (1e-6 < 1.0 - cosToSource * cosToSource) {
		// Find angles to rotate scattering plane into observer plane
		double oneOverSinScatt = 1.0 / sqrt(1.0 - cosToSource * cosToSource);
		cosEta_scatToObs = CosAngleToSource(basis.y, &pt) * oneOverSinScatt; // cosEta =  dot(refDir,y) = dot(+sunDir,y)
		sinEta_scatToObs = -CosAngleToSource(basis.z, &pt) * oneOverSinScatt; // sinEta = -dot(refDir,z) =  -dot(+sunDir,z)


	}
	else {
		// Ignore rotation to avoid divide by zero
		cosEta_scatToObs = 1.0;
		sinEta_scatToObs = 0.0;
	}

	vec.RotatePolarPlaneThru(cosEta_scatToObs, sinEta_scatToObs); // Result of pmatrix acting on source vector is in scatter frame -- rotate polarization into my frame
	vec.Normalize();

	return ok;
}

bool SKTRAN_SolarTransmission_Base::SourceTermAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source ) const
{
	bool ok = true;
	double				source_scalar;
	
//	basis.ProduceBasis(pt, look);

	ok = ok && EstimateNormalizedPolarizationVector( qobj.GetPoint(), qobj.GetLookAway(), *source, qobj.GetBasis() );
	ok = ok && SourceTermAtPoint(qobj, &source_scalar);
	*source *= source_scalar;

	return ok;
}

bool SKTRAN_SolarTransmission_Base::SourceTermAtPoint(const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source) const
{
	bool ok = true;
	double				source_scalar;

	//	basis.ProduceBasis(pt, look);

	ok = ok && EstimateNormalizedPolarizationVector(wavelength, qobj.GetPoint(), qobj.GetLookAway(), *source, qobj.GetBasis());
	ok = ok && SourceTermAtPoint(wavelength, qobj, &source_scalar);
	*source *= source_scalar;

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base::GroundSourceAtPoint		 2016- 12- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_Base::GroundSourceAtPoint( const SKTRAN_SourceTermQueryObject_Base& outboundray, double* source ) const
{
	bool ok = true;
	double brdf;
	double transmission;
	double mu_in;
	double	mu_out;
	double cosdphi;
	const HELIODETIC_POINT&        point = outboundray.GetPoint();
	HELIODETIC_UNITVECTOR		   lookv = outboundray.GetLookAway();				// Get the look direction away from the source of the ray which is towards the ground point
	HELIODETIC_UNITVECTOR		   sunvector;

	lookv.Negate();																	// Get the look direction away from the ground point	
	Sun()->UpdateSun();
	Sun()->SunUnitVector( &sunvector);
	ok = ok && TransmissionAtPoint( point, transmission);

	// Check to see if we are in the upper hemisphere
	if (point.Vector().Z() > 0)
	{
		SKTRAN_TableOpticalProperties_Base::GroundBRDFAngles(point, sunvector, lookv, &mu_in, &mu_out, &cosdphi);

		ok = ok && Integrator()->GetOpticalProps()->GetBRDF(point, mu_in, mu_out, cosdphi, &brdf);

		*source = transmission*mu_in*brdf;
	}
	else
	{
		// Solar ray does not hit the ground point
		*source = 0.0;
	}

	return ok;
}

bool SKTRAN_SolarTransmission_Base::GroundSourceAtPoint(const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& outboundray, double* source) const
{
	bool ok = true;
	double brdf;
	double transmission;
	double mu_in;
	double	mu_out;
	double cosdphi;
	const HELIODETIC_POINT&        point = outboundray.GetPoint();
	HELIODETIC_UNITVECTOR		   lookv = outboundray.GetLookAway();				// Get the look direction away from the source of the ray which is towards the ground point
	HELIODETIC_UNITVECTOR		   sunvector;

	lookv.Negate();																	// Get the look direction away from the ground point	
	Sun()->UpdateSun();
	Sun()->SunUnitVector(&sunvector);
	ok = ok && TransmissionAtPoint(wavelength, point, transmission);

	// Check to see if we are in the upper hemisphere
	if (point.Vector().Z() > 0)
	{
		SKTRAN_TableOpticalProperties_Base::GroundBRDFAngles(point, sunvector, lookv, &mu_in, &mu_out, &cosdphi);

		ok = ok && Integrator()->GetOpticalProps()->GetBRDF(wavelength, point, mu_in, mu_out, cosdphi, &brdf);

		*source = transmission * mu_in*brdf;
	}
	else
	{
		// Solar ray does not hit the ground point
		*source = 0.0;
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_Base::GroundSourceAtPoint		 2016- 12- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_Base::GroundSourceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source ) const
{
	bool ok = true;
	double source_scalar;

	ok = ok && GroundSourceAtPoint( qobj, &source_scalar );
	source->SetTo(0.0); // Depolarize at ground
	source->Assign_I ( source_scalar );

	return ok;
}

bool SKTRAN_SolarTransmission_Base::GroundSourceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source) const
{
	bool ok = true;
	double source_scalar;

	ok = ok && GroundSourceAtPoint(wavelength, qobj, &source_scalar);
	source->SetTo(0.0); // Depolarize at ground
	source->Assign_I(source_scalar);

	return ok;
}


bool SKTRAN_SolarTransmission_Base::MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const
{
	bool ok = true;
	double radiance_scalar;

	ok = ok && EstimateNormalizedPolarizationVector( qobj.GetPoint(), qobj.GetLookAway(), radiance, qobj.GetBasis());
	ok = ok && MonteCarlo_SingleScatteredRadianceAtPoint(qobj, radiance_scalar);
	radiance *= radiance_scalar;

	return ok;
}

bool SKTRAN_SolarTransmission_Base::MonteCarlo_SingleScatteredRadianceAtPoint(const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance) const
{
	bool ok = true;
	double radiance_scalar;

	ok = ok && EstimateNormalizedPolarizationVector(wavelength, qobj.GetPoint(), qobj.GetLookAway(), radiance, qobj.GetBasis());
	ok = ok && MonteCarlo_SingleScatteredRadianceAtPoint(wavelength, qobj, radiance_scalar);
	radiance *= radiance_scalar;

	return ok;
}

bool SKTRAN_SolarTransmission_Base::MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const
{
	bool ok = true;
	double radiance_scalar;

	ok = ok && MonteCarlo_GroundScatteredRadianceAtPoint(qobj, radiance_scalar);
	radiance.SetTo(0.0); // Depolarize at ground
	radiance.Assign_I ( radiance_scalar );

	return ok;
}

bool SKTRAN_SolarTransmission_Base::MonteCarlo_GroundScatteredRadianceAtPoint( const double& wavelength, const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance) const
{
	// the interface is here for consistency but no wavelength-dependent ground effects have been implemented yet
	return MonteCarlo_GroundScatteredRadianceAtPoint(qobj, radiance);
}
