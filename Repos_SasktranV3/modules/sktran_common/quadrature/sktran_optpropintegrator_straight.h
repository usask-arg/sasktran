//#include "sktran_common_internals.h"

class SKTRAN_SolarTransmission_Base;


/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Base		2014-2-7*/
/** @ingroup odintegrate
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_OpticalPropertiesIntegrator_Straight : public SKTRAN_OpticalPropertiesIntegrator_Base
{
	protected:
		virtual double									OpticalDepthOfCell							( const SKTRAN_RayOptical_Base* ray, size_t cellidx ) const;
		virtual double                                  OpticalDepthOfCell_advanceCachePoints		( const SKTRAN_RayOptical_Base* ray,
																										size_t cellidx, 
																										HELIODETIC_POINT& cacheSP, 
																										HELIODETIC_POINT& cacheEP,
																										double* kstart = nullptr,
																										double* kend = nullptr
																										) const;
		virtual double                                  OpticalDepthOfCell_withMinCache             ( const SKTRAN_RayOptical_Base* ray, size_t cellidx ) const;
		virtual bool									GetOpticalDepthFromParams                   ( double r0, double r1, double t0, double t1, double rt, double sigmak, double sigmaf, double& opticaldepth ) const;

		virtual double									OpticalDepthOfSegment						( size_t cellidx, const SKTRAN_Distance& startintercept, const SKTRAN_Distance& endintercept, const SKTRAN_RayOptical_Base* ray ) const;
		virtual double                                  OpticalDepthOfSegment_advanceCachePoints    ( size_t cellidx, 
																										const SKTRAN_Distance& startintercept, 
																										const SKTRAN_Distance& endintercept, 
																										const SKTRAN_RayOptical_Base* ray, 
																										HELIODETIC_POINT& cacheSP, 
																										HELIODETIC_POINT& cacheEP,
																										double* kstart = nullptr,
																										double* kend = nullptr
																										) const;
		virtual double                                  OpticalDepthOfSegment_withMinCache          ( size_t cellidx, const SKTRAN_RayOptical_Base* ray ) const;

		virtual bool									CalculateDistanceToTargetOpticalDepth		( double r0, double r1, double t0, double t1, double rt, double sigma0, double sigma1, double resolution, double targetTau, double& scatterDistance ) const;

		// all calls to optical property table (m_opticalprops) should go through wrappers
		virtual double									TotalExtinctionPerCM						( const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& point ) const;
		virtual bool									GetEffectiveExtinctionPerCMWithHeight1		( const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& startpoint, HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1 ) const;
		virtual bool									GetEffectiveExtinctionPerCMWithHeight1		( const SKTRAN_RayOptical_Base* ray, size_t startPtIndex, double* sigma0, double* sigma1 ) const;
	public:		
		virtual bool									CalculateRayScalarTransmission					( SKTRAN_RayOptical_Base*			ray,
																									  double*							transmisison,
																									  bool								totaltransmissiononly,
																									  bool								usecachedtransmission
																									) const override;
		virtual bool                                    CalculateRayScalarTransmission_withMinContainer   ( SKTRAN_RayOptical_Base* ray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission  ) const override;

		virtual bool									SingleScatterRadiance						( SKTRAN_RayOptical_Base* ray, double& radiance, const SKTRAN_Source_Term& solartable ) const override;
		virtual bool									FindNextScatterPosition						( double rand, double resolution, SKTRAN_TableOpticalProperties_Base* cachedOpticalProperties, SKTRAN_RayOptical_Base const* originalRay, HELIODETIC_VECTOR& scatterVector, double& scatterProb, double& randomOpticalDepth, double userImposedDistance ) const override;
		virtual bool									CalculatePartialOpticalDepth				( SKTRAN_RayOptical_Base const* ray, size_t cellindex, double fraction, double& opticaldepth ) const override;


		virtual bool									CalculateRayScalarTransmissionVector				( SKTRAN_RayOptical_Base* ray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission, std::vector<double>* sigmak = NULL, std::vector<double>* sigmaf = NULL ) const override;
};

class SKTRAN_OpticalPropertiesIntegrator_ConstantLayers : public SKTRAN_OpticalPropertiesIntegrator_Straight
{
	protected:
		virtual bool									GetOpticalDepthFromParams						( double r0, double r1, double t0, double t1, double rt, double sigmak, double sigmaf, double& opticaldepth) const override;
		virtual bool									CalculateDistanceToTargetOpticalDepth			( double r0, double r1, double t0, double t1, double rt, double sigma0, double sigma1, double resolution, double targetTau, double& scatterDistance ) const override;

};
