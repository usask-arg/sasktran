//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_WF_Integrator		2014-10-30*/
/** Calculates the weighting functions for a set of perturbations and lines
 *  of sight.  Algorithm is described in a forthcoming publication
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_WF_Integrator : private SKTRAN_OpticalPropertiesIntegrator_Adaptive
{
	private:
		std::shared_ptr< const SKTRAN_RayFactory_Base>				m_rayfactory;
		nx2dArray<double>											m_phasefunction;
		SKTRAN_GridDefScatterAngle_V21								m_phaseanglegrid;
		bool														m_isscatterer;

		std::vector<double>											m_pertalts;
		std::vector<double>											m_pertwidths;


	private:
	
		bool											AddedOpticalDepth					( const SKTRAN_HR_WF_Store&	perturbations, 
																							  std::vector<double>&													addedopticaldepth,
																							  const SKTRAN_RayOptical_Base&											ray
																							) const;


		double											PhaseFunction						( const SKTRAN_HR_WF_SpeciesInformationBase& info, size_t pertindex, double cosangle ) const;
        SKTRAN_ScatMat_MIMSNC                           PhaseMatrix                         ( const SKTRAN_HR_WF_SpeciesInformationBase& info, size_t pertindex, double cosangle ) const;

	public:
														SKTRAN_HR_WF_Integrator				()	{ m_isscatterer = false; };
		virtual										   ~SKTRAN_HR_WF_Integrator				()	{};
		
		bool											Initialize							( const SKTRAN_TableOpticalProperties_Base* opttable ) { return SetOpticalProps( opttable ); }

		bool											CreatePhaseFunctionTableParticleSize(skClimatology&														neutral,
																							 SKTRAN_AtmosphericOpticalState_V21&								opticalstate,
																							 skOpticalProperties&												optprop,
																							 double																wavel,
																							 SKTRAN_HR_WF_Store&											    perturbations,
																							 const SKTRAN_CoordinateTransform_V2&								coords,
																							 double																moderadiuschange,
																							 double																modewidthchange
																							);

		bool											CreateSolarRayFactory				( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
			                                                                                  const SKTRAN_HR_WF_Store& perturbation, 
			                                                                                  double toaHeight ); 

		bool											SetSolarRayFactory					( SKTRAN_RayFactory_Base*	rayfactory );

		bool											CalculateWeightingFunctions			( const SKTRAN_RayOptical_Base& ray,
																							  const SKTRAN_HR_WF_Store&	perturbations,
																							  const std::vector< SKTRAN_Source_Term* >&								sources,
																							  const std::vector<std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>>&               wfinfo,
																							  std::vector<SKTRAN_HR_WF_Ray>&										wfs,
																							  bool directeffectsonly = false
																							) const;

    bool											    CalculateWeightingFunctionsPolarized( const SKTRAN_RayOptical_Base& ray,
                                                                                             const SKTRAN_HR_WF_Store&	perturbations,
                                                                                             const std::vector< SKTRAN_Source_Term* >&								sources,
                                                                                             const std::vector<std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>>&               wfinfo,
                                                                                             std::vector<SKTRAN_HR_WF_Ray_Polarized>&										wfs,
                                                                                             bool directeffectsonly = false
    ) const;

		bool											CachePerturbationInformation(const SKTRAN_CoordinateTransform_V2&								coords,
			const SKTRAN_HR_WF_Store&	perturbations);
};