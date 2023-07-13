#include "sktran_montecarlo_internals.h"

class SKTRAN_ERPTManager;

class SKTRAN_ERPTCache
{
    private:
        HELIODETIC_VECTOR     m_offLosScatterPoint;
        HELIODETIC_UNITVECTOR m_higherRayDirection;
        double            m_firstScatteredRadiance;
        double            m_secondScatteredRadiance;
        double            m_higherScatteredRadiance;
        double            m_totalRadiance;
        double            m_o3ScatterWeight;
        double            p_xi1;
        double            p_xi2;
        double            p_theta1;
        double            p_theta2;
        int               m_currentOrder;

    public:
             SKTRAN_ERPTCache();
        void CachePhotonState( const SKTRAN_MCPhoton_Base& photon, int order );
        void ClearCache( );
        double EnergyToDistribute( int expectedNumMutations ) const;

        friend class SKTRAN_ERPTManager;
};


class SKTRAN_ERPTManager
{
    private:
        std::shared_ptr<const SKTRAN_RayTracer_Base>    m_raytracer;
        SKTRAN_OpticalPropertiesIntegrator_Base*		m_opi;
        SKTRAN_TableOpticalProperties_Base*				m_optprops;
    private:
        void ReleaseResources ( );

    public:
             SKTRAN_ERPTManager         ( );
            ~SKTRAN_ERPTManager         ( );
        bool TransitionProbabilityRatio ( const SKTRAN_ERPTCache& originalRay,  const SKTRAN_ERPTCache& mutationRay ) const;
        bool AcceptMutation             (       SKTRAN_ERPTCache& origToMutate, const SKTRAN_ERPTCache& partialProp ) const;
        bool SetJoinTools               ( std::shared_ptr<const SKTRAN_RayTracer_Base> rt, SKTRAN_OpticalPropertiesIntegrator_Base* opi, SKTRAN_TableOpticalProperties_Base* optprops );
};

