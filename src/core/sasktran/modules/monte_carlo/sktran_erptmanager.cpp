#include "include/sktran_montecarlo_internals.h"

//HELIODETIC_VECTOR m_losScatterPoint;
//HELIODETIC_VECTOR m_offLosScatterPoint;
//double            p_theta1;
//double            p_theta2;
//int               m_maxOrder;

SKTRAN_ERPTCache::SKTRAN_ERPTCache()
{
    ClearCache();
}

void SKTRAN_ERPTCache::CachePhotonState( const SKTRAN_MCPhoton_Base& photon, int order )
{
    m_currentOrder = order;

    m_totalRadiance = photon.photonRadiance().GetScalar();

    switch(m_currentOrder){
	case 1:
        m_firstScatteredRadiance  = photon.photonRadiance().GetRecentContribSca();
        p_xi1 = photon.m_distanceProb;
        break;
	case 2:
        m_offLosScatterPoint      = photon.m_scatterVector;
        m_secondScatteredRadiance = photon.photonRadiance().GetRecentContribSca();
        p_xi2                     = photon.m_distanceProb;
        break;
	case 3:
        m_higherRayDirection = photon.photonOptical()->LookVector();
        m_o3ScatterWeight    = 0.0<photon.ScatterWeight() ? photon.Albedo() /photon.ScatterWeight(): 0.0;
	default:
        m_higherScatteredRadiance = m_o3ScatterWeight * (m_totalRadiance - m_firstScatteredRadiance - m_secondScatteredRadiance);
        break;
	}
    
	return;
}

void SKTRAN_ERPTCache::ClearCache( )
{
    m_currentOrder            = 0;
    p_xi1                     = 1.0;
    p_xi2                     = 1.0;
    p_theta1                  = 1.0;
    p_theta2                  = 1.0;
    m_o3ScatterWeight         = 1.0;
	m_firstScatteredRadiance  = 0.0;
    m_secondScatteredRadiance = 0.0;
    m_higherScatteredRadiance = 0.0;
	m_totalRadiance           = 0.0;

    return;
}

double SKTRAN_ERPTCache::EnergyToDistribute( int expectedNumMutations ) const
{
    return m_higherScatteredRadiance / ( (double)expectedNumMutations ) ;
}

SKTRAN_ERPTManager::SKTRAN_ERPTManager( ){
    m_raytracer = nullptr;
    m_opi       = nullptr;
    m_optprops  = nullptr;
}

SKTRAN_ERPTManager::~SKTRAN_ERPTManager( ){
    ReleaseResources( );
}

void SKTRAN_ERPTManager::ReleaseResources( ){
//    if(nullptr!=m_raytracer) m_raytracer->Release(); m_raytracer = nullptr;
    if(nullptr!=m_opi)       m_opi      ->Release(); m_opi       = nullptr;
    if(nullptr!=m_optprops)  m_optprops ->Release(); m_optprops  = nullptr;
}

bool SKTRAN_ERPTManager::SetJoinTools( std::shared_ptr<const SKTRAN_RayTracer_Base> rt, SKTRAN_OpticalPropertiesIntegrator_Base* opi, SKTRAN_TableOpticalProperties_Base* optprops ){

	bool ok = true;

    ok = ok && nullptr!=rt && nullptr!=opi && nullptr!=optprops;

    if( ok ){
 //       rt      ->AddRef();
        opi     ->AddRef();
        optprops->AddRef();
        ReleaseResources ( );
        m_raytracer = rt;
        m_opi       = opi;
        m_optprops  = optprops;
	} else{
        nxLog::Record(NXLOG_ERROR, "SKTRAN_ERPTManager::SetJoinTools, Could not set join tools, received at least one nullptr.");
	}

    return ok;
}

bool SKTRAN_ERPTManager::TransitionProbabilityRatio( const SKTRAN_ERPTCache& originalRay,  const SKTRAN_ERPTCache& mutationRay ) const
{
    bool ok = true;

    return ok;
}

bool SKTRAN_ERPTManager::AcceptMutation(       SKTRAN_ERPTCache& origToMutate, const SKTRAN_ERPTCache& partialProp ) const
{
    bool ok = true;

    return ok;
}

