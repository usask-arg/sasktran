/*-----------------------------------------------------------------------------
 *					class SKTRAN_Engine_Base		2013-9-27*/
/**	@ingroup baseclass
 *	Base class for the sasktran engines.
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_Engine_Base
{
	public:
													SKTRAN_Engine_Base() {}
		virtual									   ~SKTRAN_Engine_Base() {}

		virtual bool								ConfigureModel							(       SKTRAN_SpecsUser_Base&		  modelspecifications,
																							  const SKTRAN_LineOfSightArray_V21&  linesofsight,
																							  size_t                              numthreads ) = 0;

		virtual bool								CalculateRadiance						( std::vector<SKTRAN_StokesScalar>*		losradiance,
																							  double								wavelen,
																							  size_t								numordersofscatter,
																							  SKTRAN_AtmosphericOpticalState_V21*	opticalstate,
																							  std::vector<skRTStokesVector>*        losvector=nullptr, 
																							  bool									updateclimatology = false,
																							  SKTRAN_DiagnosticInterface*			diag = NULL);

		virtual bool								CalculateMultiWavelengthRadiance		( std::vector< std::vector<SKTRAN_StokesScalar> > *		losradiance,
																							  const std::vector<double>&							wavelen,
																							  size_t												numordersofscatter,
																							  SKTRAN_AtmosphericOpticalState_V21*					opticalstate,
																							  std::vector< std::vector<skRTStokesVector> > *        losvector=nullptr, 
																							  bool													updateclimatology = false,
																							  SKTRAN_DiagnosticInterface*							diag = NULL);


};


/** \fn virtual bool SKTRAN_Engine_Base::ConfigureModel	(SKTRAN_SpecsUser_Base& modelspecifications, const SKTRAN_LineOfSightArray_V21&  linesofsight, size_t numthreads ) = 0;
*
* Initializes the model for subsequent use. This function must be called before any calls to CalculateRadiance.
*	
*  \param modelspecifications
*		Model dependent specifications. Each model has its own specifications and you must make sure
*		you pass in specifications which are appropriate for the model you are using.
*
*  \param linesofsight
*		The array of lines of sight for which radiance calculations are required.  The lines of sight are
*		passed into the model so it can perform any optimizations it feels are beneficial.
*
*  \param numthreads
*		The user can suggest how many threads are used by the processing. A value of 0 lets the software choose its own value.
*		The software may limit the number of threads actually created if the user asks for a very large number. 
*
*  \returns
*		True if everything is ok, False otherwise. You should not call CalculateRadiance unless this function succeeds
*/


/**	\fn  virtual bool SKTRAN_Engine_Base::CalculateRadiance( std::vector<SKTRAN_StokesScalar>* losradiance, double wavelen, size_t numordersofscatter, SKTRAN_AtmosphericOpticalState_V21*	opticalstate, bool updateclimatology = false, SKTRAN_DiagnosticInterface* diag = NULL) = 0;
*	Calculates the radiance at the wavelength #wavelen for the "N" lines of sight passed into the model during
*	the last call  to #ConfigureModel. 
*
*  \param losradiance
*		Returns the radiance as a vector array. The array contains the radiance calculated at
*		wavelength #wavelen for the "N" lines of sight defined in the last call to ConfigureModel. The
*		radiances are usually normalized to unit solar irradiance incident upon the atmosphere at each wavelength.
*
*  \param wavelen
*		The wavelength in nanometers of the radiance calculations.
*
*  \param numordersofscatter 
*		The number of scattering orders to compute, note a value of 1 indicates single scatter only. A value of zero
*		is reserved so the model can decide the appropriate number.
*
*  \param opticalstate
*		The object containing the collection of skOpticalProperties and skClimatologies that will be used to define
*		the optical properties of the atmosphere.
*
*  \param updateclimatology
*		An optional argument, which can be used as an efficiency when calculating radiances for the same lines of sight at multiple wavelengths.
*		If true then the climatologies in variable #opticalstate will be requested to update their caches otherwise no explicit request is made.
*
*  \param diag
*		An optional argument, it is possible to pass a diagnostic object so the model can perform callbacks during the calculatio. This is an advanced feature
*		and is normally only used by developers for debugging purposes.
*
*	\returns
*		True if everything is ok, false otherwise. Returned variables are not trustworthy unless tthis function returns true.
*/
