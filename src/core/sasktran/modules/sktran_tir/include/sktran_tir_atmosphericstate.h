
/*-----------------------------------------------------------------------------
 *				class SKTRAN_TIR_AtmosphericOpticalState		2018-06-11	*/
/** \ingroup atmosState
 * A modified version of SKTRAN_AtmosphericOpticalState_V21 for use in the
 * Thermal InfraRed engine, this class is used to by the radiative transfer
 * model to calculate absorption at any point and time in the atmosphere.
 * 
 * The main difference between this class and SKTRAN_AtmosphericOpticalState_V2
 * is that scattering is no longer computed as it is negligible in the thermal
 * IR spectral regime. Additionally, the ability to save the absorption cross
 * sections of individual species and return them to the user has been added
 * as they are needed for the calculation of analytical weighting functions
 * by the TIR engine.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TIR_AtmosphericOpticalState : public SKTRAN_AtmosphericOpticalState_V21
{
	private:
		skClimatology_UserDefinedTable*		m_perturbedatmosphericstate;	//!< The atmospheric state with a perturbed temperature for computing temperature derivatives

	private:
		bool						CalculateMultiWaveCrossSections( SKTRAN_AtmosphericOpticalStateEntry_V21*   entry,
																	 const std::vector<double>&					wavenumber,
																	 skClimatology*								neutralatmosphere,
																	 const GEODETIC_INSTANT&				    placeandtime,
																	 std::vector<double>*					    absxs,
																	 std::vector<double>*					    extxs,
																	 std::vector<double>*					    scattxs
								                                   );


	public:
									SKTRAN_TIR_AtmosphericOpticalState();
		virtual					   ~SKTRAN_TIR_AtmosphericOpticalState();
		bool						CalculateCrossSections ( const std::vector<double>&                               wavelength, 
																   std::vector<double>*                               kabs,
																   std::map<CLIMATOLOGY_HANDLE, std::vector<double>>* speciesxs,
																   std::vector<double>*								  dkabs=nullptr);	//!< Calculate absorption and species cross sections for a range of wavenumbers, using multithreading to improve performance

		bool						GetSpeciesNumberDensity(const CLIMATOLOGY_HANDLE& speciesinlist, double* numberdensity);
		bool						UpdateCache();
		bool						ContainsSpecies(const CLIMATOLOGY_HANDLE& species);
		bool						SetNumThreads(size_t numthreads);
		bool						ConfigurePerturbedAtmosphere(double perturbationfactor, const std::vector<double>& altitudes);
		bool						GetPerturbedAtmosphericStateModel(skClimatology** statemodel);
};
