#include "../sasktranv21_internals.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_V21::SKTRAN_SpecsInternal_V21		2010-5-18*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_V21::SKTRAN_SpecsInternal_V21	()
{
	m_raytracingspecs   = NULL;
	m_diffusespecs      = NULL;
	m_groundpointspecs  = NULL;
	m_quadraturespecs   = NULL;
	m_opticaltablespecs = NULL;
	m_coordinatesystem  = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_V21::~SKTRAN_SpecsInternal_V21		2010-5-18*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_V21::~SKTRAN_SpecsInternal_V21	()
{
	ReleaseResources();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_V21::ReleaseResources		2010-5-18*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_SpecsInternal_V21::ReleaseResources()
{
	if ( m_raytracingspecs   != NULL ) m_raytracingspecs  ->Release();
	if ( m_diffusespecs      != NULL ) m_diffusespecs     ->Release();
	if ( m_groundpointspecs  != NULL ) m_groundpointspecs ->Release();
	if ( m_quadraturespecs   != NULL ) m_quadraturespecs  ->Release();
	if ( m_opticaltablespecs != NULL ) m_opticaltablespecs->Release();
	

	m_raytracingspecs   = NULL;
	m_diffusespecs      = NULL;
	m_groundpointspecs  = NULL;
	m_quadraturespecs   = NULL;
	m_opticaltablespecs = NULL;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_V21::Initialize		2010-5-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsInternal_V21::Initialize( const SKTRANSO_SpecificationsUser* userspecifications )
{
	bool	ok;
	double	groundht;
	double	toaht;

	ReleaseResources();
	groundht = userspecifications->RayTracingSpecs ()->GroundAltitude();
	toaht    = userspecifications->RayTracingSpecs ()->TOAAltitude();

	ok =       userspecifications->CreateCoordinateSystem( &m_coordinatesystem, groundht, toaht );
	ok = ok && userspecifications->RayTracingSpecs ()->CreateInternalSpecs( &m_raytracingspecs   );
	ok = ok && userspecifications->DiffuseSpecs    ()->CreateInternalSpecs( &m_diffusespecs, m_coordinatesystem     );
	ok = ok && userspecifications->GroundPointSpecs()->CreateInternalSpecs( &m_groundpointspecs  );
	ok = ok && userspecifications->QuadratureSpecs ()->CreateInternalSpecs( &m_quadraturespecs   );
	ok = ok && userspecifications->OpticalGridSpecs()->CreateInternalSpecs( &m_opticaltablespecs );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsInternal_V21::Initialize, Error fetching internal specifications from user specifications. Thats a problem");
		ReleaseResources();
	}
	else
	{
		// -- Configure the ray factory for transmission only calculations within the model

		std::unique_ptr< SKTRAN_RayFactory<SKTRAN_RayOptical_Straight,
			                               SKTRAN_RayTracer_Shells,
										   SKTRAN_RayStorage_Straight> >	rayfactorytransmission ( new SKTRAN_RayFactory <SKTRAN_RayOptical_Straight,
											                                                                                SKTRAN_RayTracer_Shells,
											                                                                                SKTRAN_RayStorage_Straight>(m_coordinatesystem) );

		rayfactorytransmission->RayTracer()->Initialize( m_raytracingspecs->RayTracingShells() );
		m_rayfactory_transmissiononly = std::move( rayfactorytransmission);
//		m_rayfactory_transmissiononly->SetThisFactory( m_rayfactory_transmissiononly);

		// Configure the ray factory for trasnmission and scatter calculations within the model

		std::unique_ptr< SKTRAN_RayFactory<SKTRAN_RayOptical_Straight,
			                               SKTRAN_RayTracer_Shells,
		                                   SKTRAN_RayStorage_Straight> > rayfactory_scattered ( new SKTRAN_RayFactory<SKTRAN_RayOptical_Straight,
											                                                                          SKTRAN_RayTracer_Shells, 
											                                                                          SKTRAN_RayStorage_Straight>(m_coordinatesystem)  );

//		rayfactorytransmission->RayTracer()->Initialize( m_raytracingspecs->RayTracingShells() );
		rayfactory_scattered->RayTracer()->Initialize( m_raytracingspecs->RayTracingShells() );
		m_rayfactory_scattered = std::move( rayfactory_scattered);
//		m_rayfactory_scattered->SetThisFactory( m_rayfactory_scattered );

		std::unique_ptr< SKTRAN_RayFactory<SKTRAN_RayOptical_Straight,
			                               SKTRAN_RayTracer_Shells,
		                                   SKTRAN_RayStorage_Straight> >	rayfactory_los       ( new SKTRAN_RayFactory<SKTRAN_RayOptical_Straight, 
											                                                                             SKTRAN_RayTracer_Shells,
											                                                                             SKTRAN_RayStorage_Straight>(m_coordinatesystem) );

//		rayfactorytransmission->RayTracer()->Initialize( m_raytracingspecs->RayTracingShells() );
		rayfactory_los->RayTracer()->Initialize( m_raytracingspecs->RayTracingShells() );
		m_rayfactory_lineofsights = std::move( rayfactory_los);
//		m_rayfactory_lineofsights->SetThisFactory(m_rayfactory_lineofsights);
	}

	return ok;

}
