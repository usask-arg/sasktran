#include "../sktran_common.h"
#include <omp.h>


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::SKTRAN_SolarTransmission_2D		2013-07-23*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SolarTransmission_2D::SKTRAN_SolarTransmission_2D()
{
//	m_integrator = NULL;
//	m_raytracer  = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::~SKTRAN_SolarTransmission_2D		2013-07-23*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SolarTransmission_2D::~SKTRAN_SolarTransmission_2D()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::ReleaseResources		2013-07-23*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_SolarTransmission_2D::ReleaseResources()
{
//	if( NULL != m_integrator ) m_integrator->Release();
//	if( NULL != m_raytracer ) m_raytracer->Release();

//	m_raytracer = NULL;
//	m_integrator = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::MakeThreadSafeFor		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

 bool SKTRAN_SolarTransmission_2D::MakeThreadSafeFor( size_t numthreads )
 { 
	if( 0 < numthreads)
	{
		m_numThreadsOnFill = numthreads; 
    }

	return (0 < m_numThreadsOnFill); 
 }

 

/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::TransmissionAtPoint		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

 bool SKTRAN_SolarTransmission_2D::TransmissionAtPoint( const HELIODETIC_POINT&		point,
													    double&						transmission ) const
{
	bool ok = true;

	ok = ok && TransmissionAtVector( point.Vector(), transmission );
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::TransmissionAtVector		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_2D::TransmissionAtVector( const HELIODETIC_VECTOR&				point,
														double&									transmission ) const
{
	bool ok = true;

	double				altweights[2];
	double				szaweights[2];
	size_t				altindex[2];
	size_t				szaindex[2];

	double				alt;
	double				cossza;

	cossza = point.UnitVector().Z();
//	alt    = RayFactory()->Coords()->RadiusToAltitude( point.Magnitude() );
	alt    = m_coords->RadiusToAltitude( point.Magnitude() );

	size_t	numalt, numsza;
	ok = ok && AltWeightsForProfile( alt, altweights, altindex, numalt );
	ok = ok && CosSzaWeights( cossza, szaweights, szaindex, numsza );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_SolarTransmission_2D::TransmissionAtVector, error looking up transmission at vector");
		transmission = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{

		transmission = 0;
		for( size_t szaidx = 0; szaidx < numsza; szaidx++ )
		{

			for( size_t altidx = 0; altidx < numalt; altidx++ )
			{

				transmission += m_transmission.At( altindex[altidx], szaindex[szaidx] ) * altweights[altidx] * szaweights[szaidx];
			}
		}
	}

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::SourceTermAtPoint		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_2D::SourceTermAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj,
													 double*								source ) const
{
	bool ok = true;

	double transmission;
	double kscatt;

	ok = ok && TransmissionAtPoint( qobj.GetPoint(), transmission ); 
	ok = ok && Integrator()->GetOpticalProps()->GetScatteringCoefficientCM2( qobj.GetPoint(), CosAngleToSource( qobj.GetLookAway(), &qobj.GetPoint() ), &kscatt );
	kscatt *= 100;

	*source = transmission * kscatt;

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::GroundSourceAtPoint		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

//bool SKTRAN_SolarTransmission_2D::GroundSourceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, 
//													   double*                              source ) const
//{
//	bool ok = true;
//	
//	double radiance;
//	double albedo;
//
//	ok = ok && MonteCarlo_GroundScatteredRadianceAtPoint( qobj, radiance);
//	ok = ok && Integrator()->GetOpticalProps()->GetBRDF( qobj.GetPoint(), &albedo );
//
//	*source = radiance * albedo;
//
//	return ok;
//}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::MonteCarlo_SingleScatteredRadianceAtPoint		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_2D::MonteCarlo_SingleScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, 
																  double&                              radiance ) const
{
	bool ok = true;
	
	double transmission;
	double scattCoeff_solar;
	double kscatt;				
	ok = ok && TransmissionAtPoint( qobj.GetPoint(), transmission );
	ok = ok && Integrator()->GetOpticalProps()->GetScatteringCoefficientCM2( qobj.GetPoint(), CosAngleToSource(qobj.GetLookAway(), &qobj.GetPoint()), &scattCoeff_solar);	// New look is heliodetic --> cosScatAngle = dot(-zhat,newLook) = -newLook.z
	kscatt = Integrator()->GetOpticalProps()->ScatteringExtinctionPerCM(qobj.GetPoint());
	scattCoeff_solar = 0.0==kscatt ? 0.0 : scattCoeff_solar / kscatt;			// GetScatteringCoefficientCM2 returns P_normalized*kscatt

	radiance = transmission * scattCoeff_solar;

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::MonteCarlo_GroundScatteredRadianceAtPoint		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_2D::MonteCarlo_GroundScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, 
																  double&                              radiance ) const
{
	bool ok = true;

	double transmission;
	double scattCoeff_solar;

	ok = ok && TransmissionAtPoint( qobj.GetPoint(), transmission );
	scattCoeff_solar = qobj.GetPoint().CosSZA()/nxmath::Pi;					// Lambertian scatter --> intensity of light at lambertian scatter point ~ cos(angleBetweenSurfaceNormalAndSolarDirection)  ~=  scatterVector_{heliodetic}.dot(solarVector_{heliodetic})  =  scatterVector_{heliodetic}.dot([0 0 1]) = scatterVector_{heliodetic}.Z(); factor \pi comes from normalization over rotation about \hat{r}
	if(scattCoeff_solar < 0.0) scattCoeff_solar = 0.0;

	radiance = transmission * scattCoeff_solar;

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::AltWeightsForProfile		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_2D::AltWeightsForProfile( double alt, double* altweights, size_t* altindex, size_t& numindex ) const
{
	bool ok = true;

	size_t lowindex;
	size_t highindex;
	double upval;
	double lowval;
	m_heightgrid.IndexOfPointEqualOrAbove( alt, &highindex );
	lowindex = highindex - 1;
	if( 0==highindex )
	{
		// truncate to bottom of grid
		numindex = 1;
		altweights[0] = 1;
		altindex[0] = 0;
	}
	else if ( highindex >= m_heightgrid.NumShells() )
	{
		// truncate to top of grid
		numindex = 1;
		altweights[0] = 1;
		altindex[0] = m_heightgrid.NumShells()-1;
	}
	else
	{
		// inside the grid
		upval = m_heightgrid.At( highindex );
		lowval = m_heightgrid.At( lowindex );
		if( 1e-3 > (upval-alt) ){
			// Exactly on a grid point
			numindex       = 1;
			altweights [0] = 1.0;
			altindex   [0] = highindex;
		} else{
			numindex      = 2;
			altweights[1] = (alt - lowval) / (upval - lowval );
			altweights[0] = (upval - alt ) / (upval - lowval );
			altindex[0]   = lowindex;
			altindex[1]   = highindex;
		}
	}

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::CosSzaWeights		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_2D::CosSzaWeights( double cossza, double* szaweights, size_t* szaindex, size_t& numindex ) const
{
	bool ok = true;

	size_t lowindex;
	size_t highindex;
	double upval;
	double lowval;
	
	m_cosszagrid.IndexOfPointEqualOrAbove( cossza, &highindex );
	lowindex = highindex - 1;
	if( 0==highindex )
	{
		// outside of grid, force result to be 0
		numindex = 1;
		szaweights[0] = 0;
		szaindex[0] = 0;
	}
	else
	{
		// inside the grid
		numindex = 2;
		upval = m_cosszagrid.CosineSZA( highindex );
		lowval = m_cosszagrid.CosineSZA( lowindex );
		szaweights[1] = (cossza - lowval) / (upval - lowval );
		szaweights[0] = (upval - cossza ) / (upval - lowval );
		szaindex[0] = lowindex;
		szaindex[1] = highindex;
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::SetGeometry		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_2D::SetGeometry( std::shared_ptr< const SKTRAN_CoordinateTransform_V2>&	coords,
											   const SKTRAN_GridDefRayTracingShells_V21&				heights, 
											   const SKTRAN_GridDefCosSZA_V21&							cosszagrid )
{
	bool ok = true;

//	nxLog::Record( NXLOG_INFO,"SKTRAN_SolarTransmission_2D::SetGeometry, **** TODO **** class SKTRAN_SolarTransmission_2D and SKTRAN_SolarTransmission_3D should use only one transmission object rather than two copies, one for each class");
	m_coords  = coords;
	m_heightgrid.DeepCopy( heights );
	m_cosszagrid.DeepCopy( cosszagrid );

	// nx2dArray is column major
	m_transmission.SetSize( m_heightgrid.NumShells(), m_cosszagrid.NumAngles() );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::FillTable_ClassSpecific		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_2D::FillTable_ClassSpecific ( )
{
	bool ok = true;

	HELIODETIC_UNITVECTOR	locunit;
	HELIODETIC_VECTOR		loc;

	#pragma omp parallel for schedule(dynamic,1) private(locunit, loc)
	for( int szaidx = 0; szaidx < (int)m_cosszagrid.NumAngles(); szaidx++ )
	{
		double cossza = m_cosszagrid.CosineSZA( szaidx );
		double sinsza = sqrt( 1 - cossza*cossza );
        double transmission; 
		locunit.SetCoords(sinsza, 0, cossza );
		for( size_t altidx = 0; altidx < m_heightgrid.NumShells(); altidx++ )
		{
			loc.SetCoords( locunit, m_coords->AltitudeToRadius( m_heightgrid.At(altidx) ));
			ok = ok && CreateRayAndCalcTransmission( loc, transmission );
			m_transmission.At(altidx, szaidx) = transmission;
			//m_transmission.At(altidx,szaidx) = 1.0;
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SolarTransmission_2D::CreateRayAndCalcTransmission		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SolarTransmission_2D::CreateRayAndCalcTransmission( const HELIODETIC_VECTOR& loc, double& transmission ) const
{
	bool ok;

	HELIODETIC_UNITVECTOR						sun;
	std::unique_ptr<SKTRAN_RayOptical_Base>		ray;

	sun.SetCoords(0,0,1);
	ok =        RayFactory()->CreateRayObject( &ray);
	ok = ok &&  ray->MoveObserver( loc, sun );
	ok = ok &&  ray->TraceRay_NewMethod();
	if(ok){
		if( ray->Storage()->GroundIsHit() )
		{
			transmission = 0;
		}
		else
		{
			Integrator()->CalculateRayScalarTransmission_withMinContainer( ray.get(), NULL, true, false );
			transmission = exp(-1.0*ray->TotalOpticalDepth());
		}
	} else{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_SolarTransmission_2D::CreateRayAndCalcTransmission, Couldn't calculate transmission.");
	}

	return ok;
}
