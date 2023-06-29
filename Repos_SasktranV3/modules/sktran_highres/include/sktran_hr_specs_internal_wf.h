


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_wf		2014-10-30*/
/** Factory class responsible for creating objects required for weighting
 *  function calculations within the HR model.
 **/
/*---------------------------------------------------------------------------*/


class SKTRAN_HR_Specs_Internal_wf
{
	private:
		bool															m_doextinction;
		bool															m_doscatextinction;
		double															m_manualwfresolution;
		double															m_maxwfheight;
		double															m_wfinterpwidth;
		SKTRAN_HR_wf_aerosol_Mode										m_wfaerosolmode;
		std::vector<nxVector>											m_normals;
		std::vector<double>												m_alts;
		SKTRAN_HR_WF_Store												m_pertlist;
		SKTRAN_HR_wf_Mode												m_wfmode;
		SKTRAN_HR_wf_precision											m_wfprecision;
		std::vector<CLIMATOLOGY_HANDLE>									m_wfspecies;
		std::vector<SKTRAN_HR_wf_Species_Mode>							m_wfspeciesmode;
		std::vector<double>												m_manualwfheights;
		std::vector<double>												m_manualwfwidth;
		std::vector<double>												m_manualwfwidthleft;
		std::vector<double>												m_manualwfwidthright;

		double															m_aerosolsizepercentchange;

	private:
		bool				MakeOneDimUniformPerturbations		();

		bool				MakeOneDimLOSPerturbation			(  SKTRAN_HR_LinesOfSightTable&								linesofsight,
																   std::shared_ptr< const SKTRAN_CoordinateTransform_V2> &	coords );

		bool				MakeTwoDimPerturbation				( std::shared_ptr< const SKTRAN_CoordinateTransform_V2> &	coords,
																  const SKTRAN_HR_Specs_Internal_OpticalPropertiesTable&	opttablespecs );													

	public:
		bool				Configure							( const SKTRAN_HR_Specs_User& specs );

		bool				MakePerturbationList				( const SKTRAN_RayTracingRegionManager&						raymanager,
																  std::shared_ptr< const SKTRAN_CoordinateTransform_V2> &	coords,
																  SKTRAN_HR_LinesOfSightTable&								linesofsight,
																  const SKTRAN_HR_Specs_Internal_OpticalPropertiesTable&	opttablespecs );

		bool				DoWfCalculation						() const { return m_wfmode != SKTRAN_HR_wf_Mode_None; }

		bool				AddWfGeometryToRayTracer			( SKTRAN_RayTracer_Straight_Generic&						raytracer,
																  std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&		coords ) const;

		SKTRAN_HR_WF_Store&											PertList	()  { return m_pertlist; }
		const SKTRAN_HR_WF_Store&											PertList() const { return m_pertlist; }

		const std::vector<CLIMATOLOGY_HANDLE>&											WFSpecies	() const { return m_wfspecies; }
		const std::vector<SKTRAN_HR_wf_Species_Mode>&									WFSpeciesMode() const { return m_wfspeciesmode; }

		SKTRAN_HR_wf_precision	WFPrecision() const { return m_wfprecision; }
		SKTRAN_HR_wf_Mode WFMode() const { return m_wfmode; }
		SKTRAN_HR_wf_aerosol_Mode AerosolWFMode() const { return m_wfaerosolmode; }
		double AerosolSizePercentChange() const { return m_aerosolsizepercentchange; }
};