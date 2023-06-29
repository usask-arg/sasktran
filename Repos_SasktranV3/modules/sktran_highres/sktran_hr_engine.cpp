#include "include/sktran_hr_internals.h"
#include <omp.h>
#include <boost/stacktrace.hpp>
#include <iostream>

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Engine::SKTRAN_HR_Engine		2013-06-12*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Engine::SKTRAN_HR_Engine()
{
	m_coordinates = nullptr;
	m_calcwf = false;
	m_sun.SetInvalid();
	m_wf_extinction = nullptr;

	m_enablelosscattering = true;

	m_numcossza = 1801;
	m_prefillsolartransmission = false;
	m_usesolartableforsinglescattering = false;
}

GEODETIC_INSTANT SKTRAN_HR_Engine::ReferencePoint() const
{
	GEODETIC_INSTANT point;

	point.latitude  = m_coordinates->ReferencePtLatitude();
	point.longitude = m_coordinates->ReferencePtLongitude();
	point.heightm   = 0.0;
	point.mjd = m_coordinates->ReferencePointMJD();

	return point;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Engine::~SKTRAN_HR_Engine		2013-06-12*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Engine::~SKTRAN_HR_Engine()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Engine::ConfigureModel		2013-06-12*/
/** Initializes the model parameters.  Here we create the internal objects:
*   - Lines of Sight Table
 *  - Coordinate Transform
 *  - Optical Properties Table
 *  - Integrator
 *  - Raytracer
 *  - Diffuse Table
 *  - Solar Transmission Table
 *  - Thread Manager
 *  These objects are created based on user input from modelspecifications.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Engine::ConfigureModel( SKTRAN_SpecsUser_Base& modelspecifications, const SKTRAN_LineOfSightArray_V21& linesofsight, size_t numthreads )
{
#if SKTRAN_HR_VERBOSE_TIMING
	boost::timer::cpu_timer t;
	t.start();
#endif
	bool ok = true;

    const SKTRAN_Specifications_Base* modelspecs = dynamic_cast<const SKTRAN_Specifications_Base*>( &modelspecifications);
    ok = ok && nullptr!=modelspecs;
	const SKTRAN_HR_Specs_User* hruserspecs = dynamic_cast<const SKTRAN_HR_Specs_User*>(modelspecs);

	ok = ok && m_internalspecs.CreateCoordinates(&m_coordinates, m_sun, linesofsight, hruserspecs->RayTracingSpecsConst().GetGroundShiftAlt(), hruserspecs->RayTracingSpecsConst().GetTOAHeight(), hruserspecs->RayTracingSpecsConst().GetNadirReferencePointOnGround());
	m_internalspecs.Configure( *modelspecs, linesofsight );
	m_calcwf = m_internalspecs.CalcWf();
	m_linesofsighttable.SetLinesOfSight( linesofsight, *m_coordinates );
	ok = ok && m_internalspecs.OpticalPropertiesSpecs().CreateOpticalTable( m_opticalpropertiestable, *m_coordinates, m_internalspecs.RayTracerSpecs().TOAHeight() );
	m_opticalpropertiestable->SetCoords( m_coordinates );
	if( m_calcwf )
	{
		ok = ok && m_internalspecs.WeightingFunctionSpecs().MakePerturbationList( m_raymanager, m_coordinates, m_linesofsighttable, m_internalspecs.OpticalPropertiesSpecs() );
	}
	ok = ok && m_internalspecs.IntegratorSpecs().CreateIntegrator( *m_opticalpropertiestable, m_optintegrator, m_srcintegrator );
	ok = ok && m_internalspecs.RayTracerSpecs().CreateDiffuseRayFactory    ( m_diffuserayfactory,      m_coordinates );
	ok = ok && m_internalspecs.RayTracerSpecs().CreateLineOfSightRayFactory( m_linesofsightrayfactory, m_coordinates );
	ok = ok && m_internalspecs.RayTracerSpecs().CreateSolarRayFactory      ( m_solarrayfactory,        m_coordinates );
	ok = ok && m_linesofsighttable.CreateRays( m_linesofsightrayfactory.get() );
	ok = ok && m_internalspecs.CreateDiffuseTable( m_diffusetable );
	ok = ok && CreateSolarTable();
	ok = ok && CreateDiffuseSolarTable();
    ok = ok && CreateEmissionsTables();
	ok = ok && m_threadmanager.Initialize( m_coordinates, *m_optintegrator, *m_srcintegrator, m_diffuserayfactory, *m_diffusesolartable, *m_opticalpropertiestable );
	ok = ok && CreateDiffuseSources();
	ok = ok && m_diffusetable->Initialize( m_coordinates, *m_optintegrator, *m_srcintegrator, m_diffuserayfactory, m_diffusesources, *m_opticalpropertiestable );

	ok = ok && m_threadmanager.SetNumThreads(numthreads);

#if SKTRAN_HR_VERBOSE_TIMING
	t.stop();
	{
		boost::timer::cpu_times times = t.elapsed();
		printf("WALL Time Elapsed for Model Geometry Initialization %f \n", double(times.wall)/1E9 );
	};
#endif
	return ok;
}

bool SKTRAN_HR_Engine::CreateDiffuseSources()
{
    if( m_internalspecs.IntegratorSpecs().GetUseSolarTransmission() ) m_diffusesources.push_back( m_diffusesolartable   .get() );
    if( m_internalspecs.IntegratorSpecs().GetUseEmissions()         ) m_diffusesources.push_back( m_emissiontable.get() );
    //m_diffusesources.push_back( m_solartable.get() );

	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Engine::CalculateRadiance		2013-06-12*/
/**  Calculates the radiance for the given parameters.  The losradiance
 *   vector is resized during this call.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Engine::CalculateRadiance(	std::vector<SKTRAN_StokesScalar>*		losradiance,
											double									wavelen,
											size_t									numordersofscatter,
											SKTRAN_AtmosphericOpticalState_V21*		opticalstate,
											std::vector<skRTStokesVector>*			losvector,
											bool									updateclimatology,
											SKTRAN_DiagnosticInterface*				diag  )
{


#if SKTRAN_HR_VERBOSE_TIMING
	boost::timer::cpu_timer t;
#endif
	bool ok = true;
	GEODETIC_INSTANT point;

#if !defined(NXDEBUG)
 try
 {
#endif

	SKTRAN_HR_Diffuse_Source	diffsource;
	diffsource.SetDiffuseTable( *m_diffusetable );

	std::vector<SKTRAN_Source_Term*> sources;

    // Choose sources for the LOS rays 
    sources.resize(0);
    if( m_internalspecs.IntegratorSpecs().GetUseSolarTransmission() )
	{
		if(m_usesolartableforsinglescattering)
		{
			sources.push_back(m_diffusesolartable.get());
		}
		else
		{
			sources.push_back( m_solartable.get() );
		}
	} 
    if( m_internalspecs.IntegratorSpecs().GetUseEmissions() )         sources.push_back( m_emissiontable.get() );
	if( numordersofscatter > 1 ) sources.push_back( &diffsource );

	point.latitude  = m_coordinates->ReferencePtLatitude();
	point.longitude = m_coordinates->ReferencePtLongitude();
	point.heightm = 0.0;
	point.mjd = m_coordinates->ReferencePointMJD();

	// SetTimeAndLocation does not actually trigger a cache update, we have to call
	// ExtinctionPerCm() to force the cache update to be at the reference point
	opticalstate->SetTimeAndLocation( point, updateclimatology );
	double temp = opticalstate->ExtinctionPercm();

#if SKTRAN_HR_VERBOSE_TIMING
	t.start();
#endif
	
	ok = ok && m_opticalpropertiestable->ConfigureOptical( wavelen, *opticalstate );
    ok = ok && m_emissiontable->         ConfigureOptical( wavelen, opticalstate->EmissionObjectVar(), m_opticalpropertiestable.get() );

	ok = ok && m_linesofsightrayfactory->ConfigureOptical(opticalstate, wavelen, point);
	ok = ok && m_diffuserayfactory->ConfigureOptical(opticalstate, wavelen, point);
	ok = ok && m_solarrayfactory->ConfigureOptical(opticalstate, wavelen, point);

	m_diffusesolartable->FillTable();

#if SKTRAN_HR_VERBOSE_TIMING
	t.stop();
	{
		boost::timer::cpu_times times = t.elapsed();
		printf("WALL Time Elapsed for Optical Properties Table Initialization %f \n", double(times.wall)/1E9 );
	};
#endif

	losradiance->resize( m_linesofsighttable.NumRays() );
	if(nullptr!=losvector) losvector->resize( m_linesofsighttable.NumRays() );

	if( numordersofscatter > 1 || (m_calcwf && m_internalspecs.WeightingFunctionSpecs().WFMode() == SKTRAN_HR_wf_Mode_2d ) )
	{
#if SKTRAN_HR_VERBOSE_TIMING
		t.start();
#endif
#if SKTRAN_HR_VERBOSE_TIMING
		t.stop();
		{
			boost::timer::cpu_times times = t.elapsed();
			printf("WALL Time Elapsed for Solar Transmission Table Creation %f \n", double(times.wall)/1E9 );
		}
#endif
	}

#if SKTRAN_HR_VERBOSE_TIMING
	t.start();
#endif

    
	ok = ok && m_threadmanager.ComputeFieldCPU( m_diffusetable.get(), numordersofscatter, wavelen );


#if SKTRAN_HR_VERBOSE_TIMING
	t.stop();
	{
		boost::timer::cpu_times times = t.elapsed();
		printf("WALL Time Elapsed for Diffuse Radiance Calculation %f \n", double(times.wall)/1E9 );
	}
	t.start();
#endif

	//TODO: move this into thread manager
	#pragma omp parallel for schedule(dynamic, 1) reduction(&&:ok)
	for( int i = 0; i < (int)m_linesofsighttable.NumRays(); i++ )
	{
		ok = ok && m_linesofsighttable.RayEntryAt(i)->TraceRay_NewMethod();
		ok = ok && m_optintegrator->CalculateRayScalarTransmission_withMinContainer( m_linesofsighttable.RayAt(i), &losradiance->at(i), false, false);			
		if(nullptr!=losvector){
			SKTRAN_Stokes_NC stokes;
			ok = ok && m_srcintegrator->IntegrateSourceTerm( m_linesofsighttable.RayAt(i), stokes, sources );
			losvector->at(i).SetTo( stokes.I(), stokes.Q(), stokes.U(), stokes.V() );
			losradiance->at(i) = losvector->at(i).At(1);
		} else{
			ok = ok && m_srcintegrator->IntegrateSourceTerm( m_linesofsighttable.RayAt(i), losradiance->at(i), sources );
		}

		if( SKTRAN_HR_ADD_SOLAR_DISK_TO_LOS )
		{
			double cosangle = m_linesofsighttable.RayAt(i)->LookVector().Z();
			double sinangle = sqrt( 1.0 - cosangle*cosangle);
			double opticaldepth = m_linesofsighttable.RayAt(i)->TotalOpticalDepth();
			double sunradius = 0.00929826 /2;
			if ( m_linesofsighttable.RayAt(i)->Storage()->GroundIsHit() == false )
			{
				losradiance->at(i) += exp(-1.0*opticaldepth) * std::max(0.0, sunradius-sinangle)/sunradius; 
			}
		}


		if( SKTRAN_HR_DUMP_LOS_DIFFUSE && i==0 )
			DumpRay( m_linesofsighttable.RayAt(i), *m_opticalpropertiestable ,*sources[0] );
		if(!ok) nxLog::Record( NXLOG_ERROR, "SKTRAN_HR_Engine::CalculateRadiance, Could not integrate LOS %i.", i );
	}

#if SKTRAN_HR_VERBOSE_TIMING
	t.stop();
	if( SKTRAN_HR_VERBOSE_TIMING )
	{
		boost::timer::cpu_times times = t.elapsed();
		printf("WALL Time Elapsed for Line of Sight Integration %f \n", double(times.wall)/1E9 );
	}
#endif

	if( m_calcwf )
	{	
		std::vector<SKTRAN_Source_Term*> wfsources;
		if( numordersofscatter > 1 )
		{
			wfsources.resize(2);
			wfsources[0] = m_diffusesolartable.get();
			wfsources[1] = &diffsource;
		}
		else if ( m_internalspecs.WeightingFunctionSpecs().WFMode() == SKTRAN_HR_wf_Mode_2d )
		{
			wfsources.resize(1);
			wfsources[0] = m_diffusesolartable.get();
		}
		else if ( m_usesolartableforsinglescattering )
		{
			wfsources.resize(1);
			wfsources[0] = m_diffusesolartable.get();
		}
		else
		{
			wfsources.resize(1);
			wfsources[0] = m_solartable.get();
		}
		if( m_internalspecs.WeightingFunctionSpecs().WFSpecies().size() == 0 )
		{
			nxLog::Record(NXLOG_WARNING, "ERROR: Calculating weighting functions without setting the weighting function species");
			m_diffusetable->CleanDiffuseIndexes();
			return ok;
		}

        if(m_internalspecs.MaxPolarizationOrder() > 0) {
            CalculateWeightingFunctionsPolarized( wavelen, wfsources, *opticalstate );
        } else {
            CalculateWeightingFunctions( wavelen, wfsources, *opticalstate );
        }
	}

	m_diffusetable->CleanDiffuseIndexes();

#if !defined(NXDEBUG)
 }
 catch (const std::exception &exc)
 {
	 printf("**** SKTRAN_HR_Engine::CalculateRadiance has thrown an unhandled execption");
	 std::cout << boost::stacktrace::stacktrace();
 }
#endif
	return ok;
}
/* Deprecated 
bool SKTRAN_HR_Engine::CalculateSecondOrderRadiance( std::vector<SKTRAN_StokesScalar>* losradiance, double wavelen, size_t numordersofscatter, SKTRAN_AtmosphericOpticalState_V21* opticalstate, bool updateclimatology, SKTRAN_DiagnosticInterface* diag  )
{
	bool ok = true;
	GEODETIC_INSTANT point;
	SKTRAN_HR_Diffuse_Second_Order_Source* diffsource = new SKTRAN_HR_Diffuse_Second_Order_Source;

	m_internalspecs.DiffuseSpecs().MakeSecondOrderSource( &diffsource );

	std::vector<SKTRAN_Source_Term*> sources;
	sources.resize(2);
	sources[0] = m_solartable;
	sources[1] = diffsource;

	point.latitude = m_coordinates->ReferencePtLatitude();
	point.longitude = m_coordinates->ReferencePtLongitude();
	point.heightm = 0.0;
	point.mjd  = m_linesofsighttable.MeanMJD();

	opticalstate->SetTimeAndLocation( point, updateclimatology );

	m_opticalpropertiestable->ConfigureOptical( wavelen, *opticalstate );
	losradiance->resize( m_linesofsighttable.NumRays() );

	diffsource->SetOptical( *m_diffuseraytracer, *m_integrator, *m_opticalpropertiestable, *m_solartable );

	SKTRAN_RayTracer_Shells_Curved* curvedtracer = NULL;
	curvedtracer = dynamic_cast<SKTRAN_RayTracer_Shells_Curved*>(m_linesofsightraytracer);
	if( NULL != curvedtracer )
	{
		curvedtracer->CalcRefractiveProfile( *opticalstate, wavelen, opticalstate->GetTimeAndLocation() );
	}

	for( size_t i = 0; i < m_linesofsighttable.NumRays(); i++ )
	{
		m_linesofsightraytracer->CreateRay( *m_linesofsighttable.RayAt(i));
		m_integrator->CalculateRayTransmission( m_linesofsighttable.RayAt(i), &losradiance->at(i), false, false);	
		m_integrator->IntegrateSourceTerm( m_linesofsighttable.RayAt(i), losradiance->at(i), sources );
	}

	diffsource->Release();
	return ok;

}
*/

bool SKTRAN_HR_Engine::CalculateWeightingFunctions( double wavelen, const std::vector< SKTRAN_Source_Term* > sources, SKTRAN_AtmosphericOpticalState_V21& opticalstate )
{
	boost::timer::cpu_timer t;
	t.start();

	SKTRAN_HR_WF_Store& pertlist = m_internalspecs.WeightingFunctionSpecs().PertList();

	auto wfspecies = m_internalspecs.WeightingFunctionSpecs().WFSpecies();

    m_wf.clear();
    m_wf.resize(m_linesofsighttable.NumRays());

	for (int idray = 0; idray < m_linesofsighttable.NumRays(); idray++)
	{
		m_wf[idray].clear();
		m_wf[idray].resize(wfspecies.size());
	}

	SKTRAN_HR_WF_Integrator wfintegrator;
	wfintegrator.Initialize( m_opticalpropertiestable.get() );
	skClimatology* neutral;
	opticalstate.GetAtmosphericStateModel( &neutral );

	std::vector<std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>> wfinfo;
	for (int idspecies = 0; idspecies < wfspecies.size(); idspecies++)
	{
		skOpticalProperties* optprop;
		skClimatology* clim;
		opticalstate.GetSpeciesOpticalProperties(wfspecies[idspecies], &optprop);
		opticalstate.GetSpeciesClimatology(wfspecies[idspecies], &clim);

		if (m_internalspecs.WeightingFunctionSpecs().WFSpeciesMode()[idspecies] == SKTRAN_HR_wf_Species_Mode::SKTRAN_HR_wf_Species_numberdensity)
		{
			wfinfo.push_back(std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>(new SKTRAN_HR_WF_SpeciesInformationStandard(*neutral, *optprop, pertlist, *m_coordinates, wavelen, m_internalspecs.OpticalPropertiesSpecs())));
		}
		else if (m_internalspecs.WeightingFunctionSpecs().WFSpeciesMode()[idspecies] == SKTRAN_HR_wf_Species_Mode::SKTRAN_HR_wf_Species_LogNormal_ModeRadius)
		{
			wfinfo.push_back(std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>(new SKTRAN_HR_WF_SpeciesInformationAerosolLogNormal(*neutral, *optprop, pertlist, *m_coordinates, wavelen, opticalstate, clim, wfspecies[idspecies], m_internalspecs.WeightingFunctionSpecs().AerosolSizePercentChange() / 100, 0, m_internalspecs.OpticalPropertiesSpecs())));
		}
		else if (m_internalspecs.WeightingFunctionSpecs().WFSpeciesMode()[idspecies] == SKTRAN_HR_wf_Species_Mode::SKTRAN_HR_wf_Species_LogNormal_ModeWidth)
		{
			wfinfo.push_back(std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>(new SKTRAN_HR_WF_SpeciesInformationAerosolLogNormal(*neutral, *optprop, pertlist, *m_coordinates, wavelen, opticalstate, clim, wfspecies[idspecies], 0.0, m_internalspecs.WeightingFunctionSpecs().AerosolSizePercentChange() / 100, m_internalspecs.OpticalPropertiesSpecs())));
		}
	}

	wfintegrator.CachePerturbationInformation(*m_coordinates, pertlist);

	SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
		               SKTRAN_RayTracer_Straight_Generic,
					   SKTRAN_RayStorage_Straight_HR>*		solarrayfactory  =  new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
																									   SKTRAN_RayTracer_Straight_Generic,
																									   SKTRAN_RayStorage_Straight_HR> (m_coordinates);

	solarrayfactory->RayTracer()->SetEarthRadius( m_coordinates->AltitudeToRadius(0.0) );
	solarrayfactory->RayTracer()->SetUpperAtmoRadius( m_coordinates->AltitudeToRadius(0.0)+m_internalspecs.RayTracerSpecs().TOAHeight() );
	m_internalspecs.WeightingFunctionSpecs().AddWfGeometryToRayTracer( *solarrayfactory->RayTracer(), m_coordinates );
	
	wfintegrator.SetSolarRayFactory( solarrayfactory );
	#pragma omp parallel for schedule(dynamic, 1)
	for( int i = 0; i < (int) m_linesofsighttable.NumRays(); i++ )
	{
		if (m_internalspecs.WeightingFunctionSpecs().WFPrecision() == SKTRAN_HR_wf_precision_all)
		{
			wfintegrator.CalculateWeightingFunctions(*m_linesofsighttable.RayAt(i), pertlist, sources, wfinfo, m_wf[i]);
		}
		else
		{
			// Approx wf
			wfintegrator.CalculateWeightingFunctions(*m_linesofsighttable.RayAt(i), pertlist, sources, wfinfo, m_wf[i], true);
		}
	}

	for (size_t pertidx = 0; pertidx < pertlist.StoreSize(); pertidx++)
	{
		for (size_t speciesidx = 0; speciesidx < wfspecies.size(); speciesidx++)
		{

			double extxs;
			extxs = wfinfo[speciesidx]->CrossSection()[pertidx];
			for (size_t rayidx = 0; rayidx < m_wf.size(); rayidx++)
			{
				m_wf[rayidx][speciesidx].WeightAt(pertidx) = m_wf[rayidx][speciesidx].WeightAt(pertidx) * extxs * 100; // convert /cm to /m

				if (m_wf[rayidx][speciesidx].WeightAt(pertidx) != m_wf[rayidx][speciesidx].WeightAt(pertidx))
				{
					// nan value, set it to 0
					m_wf[rayidx][speciesidx].WeightAt(pertidx) = 0.0;
				}
			}
		}
	}

	t.stop();
	if( SKTRAN_HR_VERBOSE_TIMING )
	{
		boost::timer::cpu_times times = t.elapsed();
		printf("WALL Time Elapsed for WF Calculation %f \n", double(times.wall)/1E9 );
	}
	
	return true;
}


bool SKTRAN_HR_Engine::CalculateWeightingFunctionsPolarized(double wavelen,
                                                            const std::vector<SKTRAN_Source_Term *> sources,
                                                            SKTRAN_AtmosphericOpticalState_V21 &opticalstate) {
    boost::timer::cpu_timer t;
    t.start();

    SKTRAN_HR_WF_Store& pertlist = m_internalspecs.WeightingFunctionSpecs().PertList();

    auto wfspecies = m_internalspecs.WeightingFunctionSpecs().WFSpecies();

    m_wf_pol.clear();
    m_wf_pol.resize(m_linesofsighttable.NumRays());

    m_wf.clear();
    m_wf.resize(m_linesofsighttable.NumRays());

    for (int idray = 0; idray < m_linesofsighttable.NumRays(); idray++)
    {
        m_wf_pol[idray].clear();
        m_wf_pol[idray].resize(wfspecies.size());

        m_wf[idray].clear();
        m_wf[idray].resize(wfspecies.size());
    }

    SKTRAN_HR_WF_Integrator wfintegrator;
    wfintegrator.Initialize( m_opticalpropertiestable.get() );
    skClimatology* neutral;
    opticalstate.GetAtmosphericStateModel( &neutral );

    std::vector<std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>> wfinfo;
    for (int idspecies = 0; idspecies < wfspecies.size(); idspecies++)
    {
        skOpticalProperties* optprop;
        skClimatology* clim;
        opticalstate.GetSpeciesOpticalProperties(wfspecies[idspecies], &optprop);
        opticalstate.GetSpeciesClimatology(wfspecies[idspecies], &clim);

        if (m_internalspecs.WeightingFunctionSpecs().WFSpeciesMode()[idspecies] == SKTRAN_HR_wf_Species_Mode::SKTRAN_HR_wf_Species_numberdensity)
        {
            wfinfo.push_back(std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>(new SKTRAN_HR_WF_SpeciesInformationStandard(*neutral, *optprop, pertlist, *m_coordinates, wavelen, m_internalspecs.OpticalPropertiesSpecs())));
        }
        else if (m_internalspecs.WeightingFunctionSpecs().WFSpeciesMode()[idspecies] == SKTRAN_HR_wf_Species_Mode::SKTRAN_HR_wf_Species_LogNormal_ModeRadius)
        {
            wfinfo.push_back(std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>(new SKTRAN_HR_WF_SpeciesInformationAerosolLogNormal(*neutral, *optprop, pertlist, *m_coordinates, wavelen, opticalstate, clim, wfspecies[idspecies], m_internalspecs.WeightingFunctionSpecs().AerosolSizePercentChange() / 100, 0, m_internalspecs.OpticalPropertiesSpecs())));
        }
        else if (m_internalspecs.WeightingFunctionSpecs().WFSpeciesMode()[idspecies] == SKTRAN_HR_wf_Species_Mode::SKTRAN_HR_wf_Species_LogNormal_ModeWidth)
        {
            wfinfo.push_back(std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>(new SKTRAN_HR_WF_SpeciesInformationAerosolLogNormal(*neutral, *optprop, pertlist, *m_coordinates, wavelen, opticalstate, clim, wfspecies[idspecies], 0.0, m_internalspecs.WeightingFunctionSpecs().AerosolSizePercentChange() / 100, m_internalspecs.OpticalPropertiesSpecs())));
        }
    }

    wfintegrator.CachePerturbationInformation(*m_coordinates, pertlist);

    SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
            SKTRAN_RayTracer_Straight_Generic,
            SKTRAN_RayStorage_Straight_HR>*		solarrayfactory  =  new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
            SKTRAN_RayTracer_Straight_Generic,
            SKTRAN_RayStorage_Straight_HR> (m_coordinates);

    solarrayfactory->RayTracer()->SetEarthRadius( m_coordinates->AltitudeToRadius(0.0) );
    solarrayfactory->RayTracer()->SetUpperAtmoRadius( m_coordinates->AltitudeToRadius(0.0)+m_internalspecs.RayTracerSpecs().TOAHeight() );
    m_internalspecs.WeightingFunctionSpecs().AddWfGeometryToRayTracer( *solarrayfactory->RayTracer(), m_coordinates );

    wfintegrator.SetSolarRayFactory( solarrayfactory );
#pragma omp parallel for schedule(dynamic, 1)
    for( int i = 0; i < (int) m_linesofsighttable.NumRays(); i++ )
    {
        if (m_internalspecs.WeightingFunctionSpecs().WFPrecision() == SKTRAN_HR_wf_precision_all)
        {
            wfintegrator.CalculateWeightingFunctionsPolarized(*m_linesofsighttable.RayAt(i), pertlist, sources, wfinfo, m_wf_pol[i]);
        }
        else
        {
            // Approx wf
            wfintegrator.CalculateWeightingFunctionsPolarized(*m_linesofsighttable.RayAt(i), pertlist, sources, wfinfo, m_wf_pol[i], true);
        }
    }

    for (size_t speciesidx = 0; speciesidx < wfspecies.size(); speciesidx++)
    {

        for (size_t rayidx = 0; rayidx < m_wf_pol.size(); rayidx++) {
            m_wf[rayidx][speciesidx].Allocate(pertlist.StoreSize());
        }
    }


    for (size_t pertidx = 0; pertidx < pertlist.StoreSize(); pertidx++)
    {
        for (size_t speciesidx = 0; speciesidx < wfspecies.size(); speciesidx++)
        {

            double extxs;
            extxs = wfinfo[speciesidx]->CrossSection()[pertidx];
            for (size_t rayidx = 0; rayidx < m_wf_pol.size(); rayidx++)
            {
                m_wf_pol[rayidx][speciesidx].WeightAt(pertidx) = m_wf_pol[rayidx][speciesidx].WeightAt(pertidx) * extxs * 100; // convert /cm to /m

                if (m_wf_pol[rayidx][speciesidx].WeightAt(pertidx).I() != m_wf_pol[rayidx][speciesidx].WeightAt(pertidx).I())
                {
                    // nan value, set it to 0
                    m_wf_pol[rayidx][speciesidx].WeightAt(pertidx).SetTo(0.0);
                }
                m_wf[rayidx][speciesidx].WeightAt(pertidx) = m_wf_pol[rayidx][speciesidx].WeightAt(pertidx).I();
            }
        }
    }

    t.stop();
    if( SKTRAN_HR_VERBOSE_TIMING )
    {
        boost::timer::cpu_times times = t.elapsed();
        printf("WALL Time Elapsed for WF Calculation %f \n", double(times.wall)/1E9 );
    }

    return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Engine::CreateSolarTable		2013-06-27*/
/**  Creates the internal solar transmission source used for all calculations
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Engine::CreateSolarTable()
{
	bool ok = true;
	std::unique_ptr<SKTRAN_Sun_Base> sun ( new SKTRAN_Sun_Point );
	//SKTRAN_SolarTransmission_NoTable_reuseRays* solartable = new SKTRAN_SolarTransmission_NoTable_reuseRays;

	if (m_enablelosscattering)
	{
		std::unique_ptr<SKTRAN_SolarTransmission_NoTable> solartable(new SKTRAN_SolarTransmission_NoTable);
		solartable->ConfigureOptical(m_solarrayfactory, m_optintegrator.get());
		m_solartable = std::move(solartable);
	}
	else
	{
		std::unique_ptr<SKTRAN_SolarTransmission_NoTable_NoLOSSource> solartable(new SKTRAN_SolarTransmission_NoTable_NoLOSSource);
		solartable->ConfigureOptical(m_solarrayfactory, m_optintegrator.get());
		m_solartable = std::move(solartable);
	}

	m_solartable->AddRef();


	m_solartable->SetSun( sun );

	return ok;
}


bool SKTRAN_HR_Engine::CreateEmissionsTables()
{
    bool ok = true;

	std::vector<double> heights = m_internalspecs.RayTracerSpecs().SolarShellHeights();
	SKTRAN_GridDefRayTracingShells_V21 altgrid;
	altgrid.ConfigureHeights( &heights[0], heights.size() );
	if (!m_internalspecs.RayTracerSpecs().UseManualSolarShells())
	{
		altgrid.SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
	}
	// and some cossza
	SKTRAN_GridDefCosSZA_V21	cosszagrid;
	double refcossza = m_coordinates->ReferencePoint(0).CosSZA();
	cosszagrid.AllocateGridArray( 1 );
    cosszagrid.AtVar(0) = refcossza;
	cosszagrid.SetGridSearchMode( SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM );
	SKTRAN_GridDefSLON_V21		slongrid; 
    slongrid.AllocateGridArray( 1 );
    slongrid.AtVar(0) = 0.0;

    if( m_internalspecs.IntegratorSpecs().GetUseEmissions() ){
        m_emissiontable = std::unique_ptr<SKTRAN_EmissionTable_Base> (new SKTRAN_EmissionTable_1D );
    } else{
        m_emissiontable = std::unique_ptr<SKTRAN_EmissionTable_Base> (new SKTRAN_EmissionTable_DoNothing );
    }

    //std::unique_ptr<SKTRAN_SolarTransmission_3D> solartable ( new SKTRAN_SolarTransmission_3D );
	//ok = ok && solartable->InitializeGeometry( m_coordinates, altgrid, cosszagrid, slongrid );
	//ok = ok && solartable->ConfigureOptical( m_solarrayfactory, m_optintegrator.get() );
	//ok = ok && solartable->SetSun( sun );
	//ok = ok && solartable->FillTable();
	//m_diffusesolartable = std::move( solartable );

	//m_diffusesolartable->AddRef();

    ok = ok && nullptr!=m_emissiontable;
    if(ok){
        m_emissiontable->AddRef(); 
        ok = ok && m_emissiontable->InitializeGeometry( m_coordinates, altgrid, cosszagrid, slongrid ); 
    }
    
    if(!ok) nxLog::Record( NXLOG_WARNING, "SKTRAN_HR_Engine::CreateEmissionsTables, Could not create emission tables. ");
    
    return ok;
}


//bool SKTRAN_HR_Engine::CreateEmissionsTables()
//{
//    bool ok = true;
//
//    std::unique_ptr< SKTRAN_EmissionTable_NoTable > emissiontable( new SKTRAN_EmissionTable_NoTable );
//    emissiontable->AddRef( );
//    std::unique_ptr< SKTRAN_EmissionTable_NoTable > diffuseemissiontable( new SKTRAN_EmissionTable_NoTable );
//    diffuseemissiontable->AddRef( );
//
//    emissiontable->SetTableFile( std::string("C:/ARGsoftware/a-band-spectra/tabulatedEmission_nxFormat.txt" ) );
//    emissiontable->InitializeGeometry( m_coordinates.get() );
//    //emissiontable->SetGeometry( altgrid, refpt );
//    diffuseemissiontable->SetTableFile(std::string("C:/ARGsoftware/a-band-spectra/tabulatedEmission_nxFormat.txt"));
//    diffuseemissiontable->InitializeGeometry( m_coordinates.get() );
//    //diffuseemissiontable->SetGeometry( altgrid, refpt );
//
//    emissiontable        ->Initialize( );
//    diffuseemissiontable ->Initialize( );
//
//    m_emissiontable        = std::move( emissiontable );
//    m_diffuseemissiontable = std::move( diffuseemissiontable );
//
//    return ok;
//}


bool SKTRAN_HR_Engine::CreateDiffuseSolarTable()
{
	double maxcossza, mincossza;
	size_t numcossza;

	bool ok = true;
    std::unique_ptr<SKTRAN_Sun_Base> sun ( new SKTRAN_Sun_Point );
	std::unique_ptr<SKTRAN_SolarTransmission_3D> solartable ( new SKTRAN_SolarTransmission_3D(m_prefillsolartransmission, m_internalspecs.RayTracerSpecs().SolarRefractionEnabled()) );

	std::vector<double> heights = m_internalspecs.RayTracerSpecs().SolarShellHeights();
	SKTRAN_GridDefRayTracingShells_V21 altgrid;
	altgrid.ConfigureHeights( &heights[0], heights.size() );
	if (!m_internalspecs.RayTracerSpecs().UseManualSolarShells())
	{
		altgrid.SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
	}
	// and some cossza
	SKTRAN_GridDefCosSZA_V21	cosszagrid;
	double refcossza = m_coordinates->ReferencePoint(0).CosSZA();

	double refsza, maxlossza, minlossza;
	m_internalspecs.RayManager().GetSZA(&refsza, &minlossza, &maxlossza);

	if(m_prefillsolartransmission && m_internalspecs.ScatterOrder() == 1 && !m_internalspecs.RayTracerSpecs().LineOfSightRefractionEnabled() && !m_internalspecs.RayTracerSpecs().SolarRefractionEnabled())
	{
		// Have to do some annoying calculations to get TOA SZA that corresponds to SZA at the earth
		double groundradius = m_coordinates->AltitudeToRadius(m_coordinates->GroundAltitude());
		double toaradius = m_coordinates->AltitudeToRadius(m_coordinates->TOAAltitude());
		double h = (toaradius - groundradius); // altitude of atmosphere

		// quadratic equation for distance from ground to TOA vertically at sza=minlossza, a=1
		double b = (2.0*groundradius*nxmath::cosd(minlossza));
		double c = -2.0*groundradius*h - h*h;

		double x = (-b + std::sqrt(b*b - 4*c)) / 2.0;
		// deviation in sza 
		double beta = nxmath::asind(x * nxmath::sind(minlossza) / (toaradius));

		double mintablesza = std::max(0.0, minlossza - 1 - beta);
		double maxtablesza = std::min(180.0, maxlossza + 1);

		double prefillmaxcossza = nxmath::cosd(mintablesza);
		double prefillmincossza = nxmath::cosd(maxtablesza);
		double deltacossza = 2.0 / double(m_numcossza);
		int prefillnumcossza = std::max(int((prefillmaxcossza - prefillmincossza) / deltacossza), 10);
		// Scale the numcossza by sza
		double scalefactor;
		if(refsza < 88)
		{
			scalefactor = 1.0 / (5.0 * nxmath::cosd(refsza));
		}
		else
		{
			scalefactor = 1.0 / (5.0 * nxmath::cosd(88));
		}
		prefillnumcossza = int(scalefactor * double(prefillnumcossza));
		// std::printf("%f -- %f -- %f -- %f -- %f -- %f -- %d\n", minlossza, maxlossza, maxcossza, mincossza, beta, scalefactor, prefillnumcossza);

		SKTRAN_GridDefCosSZA_V21	prefillcosszagrid;
		prefillcosszagrid.AllocateGridArray( prefillnumcossza );
		for( size_t szaidx = 0; szaidx < prefillnumcossza; szaidx++ )
		{
			prefillcosszagrid.AtVar( szaidx ) = prefillmincossza + ( (prefillmaxcossza - prefillmincossza) * double(szaidx) / (double(prefillnumcossza) -1 ));
		}
		prefillcosszagrid.SetGridSearchMode( SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM );
		solartable->SetPrefillGrid(prefillcosszagrid);

		maxcossza = 1;
		mincossza = -1;
		numcossza = m_numcossza;
	}
	else
	{
		maxcossza = 1;
		mincossza = -1;
		numcossza = m_numcossza;
	}

	cosszagrid.AllocateGridArray( numcossza );
	for( size_t szaidx = 0; szaidx < numcossza; szaidx++ )
	{
		cosszagrid.AtVar( szaidx ) = mincossza + ( (maxcossza - mincossza) * double(szaidx) / (double(numcossza) -1 ));
	}
	cosszagrid.SetGridSearchMode( SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM );
	

	SKTRAN_GridDefSLON_V21		slongrid = m_internalspecs.OpticalPropertiesSpecs().MakeSLONGrid( *m_coordinates );
	ok = ok && solartable->InitializeGeometry( m_coordinates, altgrid, cosszagrid, slongrid );
	ok = ok && solartable->ConfigureOptical( m_solarrayfactory, m_optintegrator.get() );
	ok = ok && solartable->SetSun( sun );
	m_diffusesolartable = std::move( solartable );

	m_diffusesolartable->AddRef();


	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Engine::ReleaseResources		2013-06-12*/
/**  Have to double release 
 *
 **/
/*---------------------------------------------------------------------------*/

void SKTRAN_HR_Engine::ReleaseResources()
{
	if( nullptr != m_opticalpropertiestable )
	{
		m_opticalpropertiestable->Release();
		m_opticalpropertiestable.release();
	}
//	if( nullptr != m_linesofsightrayfactory )
//	{
//		m_linesofsightrayfactory->Release();
//		m_linesofsightrayfactory.release();
//	}
//	if( nullptr != m_diffuserayfactory ) 
//	{
//		m_diffuserayfactory->Release();
//		m_diffuserayfactory.release();
//	}
//	if( nullptr != m_solarrayfactory )
//	{
//		m_solarrayfactory->Release();
//		m_solarrayfactory.release();
//	}
	if( nullptr != m_solartable )
	{
		m_solartable->Release();
		m_solartable.release();
	}
    if( nullptr != m_diffusesolartable )
    {
        m_diffusesolartable->Release();
        m_diffusesolartable.release();
    }

    if( nullptr != m_emissiontable )
    {
        m_emissiontable->Release();
        m_emissiontable.release();
    }

	if( nullptr != m_optintegrator )
	{
		m_optintegrator->Release();
		m_optintegrator.release();
	}
	if( nullptr != m_srcintegrator )
	{
		m_srcintegrator->Release();
		m_srcintegrator.release();
	}
}	

