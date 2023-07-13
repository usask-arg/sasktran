


/*-----------------------------------------------------------------------------
 *				skOpticalProperties_MartHybridProfile		2013-1-25*/
/**	\ingroup aerosolskopticalprop
 *	Class for merging the phase matrices of aerosol particles and ice crystals
 *	in the UTLS region.  Scattering properties are computed for both according
 *	to the specified size distributions, with particle number densities set up
 *	in their own climatologies.
 *
 *	A class that makes a hybrid aerosol-ice crystal species.  This species is
 *	to be used in the MART O3 retrieval to retrieve the extinction due to aerosol
 *	and ice crystals.
 *	This class uses an optical properties height profile that is a mixture of
 *	ice crystal and aerosol scattering/absorption characteristics.  A (temporary
 *	input) parameter to this class is the tropopause height. To use the idea of a
 *	UTLS region, in which characteristics of both troposphere and stratosphere are
 *	used, a layer
 *	To retrieve the extinction, the scattering cross section is set to unity.
 */
/*---------------------------------------------------------------------------*/

class skOpticalProperties_MartHybridProfile : public skOpticalProperties
{
	private:
		double										m_tropopauseheightkm;			//!< height of the tropopause (in km) at reference point
		double										m_utlsthicknesskm;				//!< thickness of the UTLS layer (in km) at reference point
		double										m_sulphatefraction;				//!< Fraction of sulphate applied to cross-sections
		double										m_cirrusfraction;				//!< Fraction of cirrus applied to cross-sections;			
		skOpticalProperties_AerosolProfileH2SO4*	m_sulphateaerosol;				//!< calculates scattering properties of sulphate aerosols. User can modify mode radius and width
		skOpticalProperties_AerosolProfileIce*		m_icecirrus;					//!< calculates scattering properties of ice cirrus. User can modify effective radius and gamma rate.

	private:
		void										init										();
		void										ReleaseObjects								();
		bool										UpdateComponentFractions					( double h_km );
		skClimatology_UserDefinedTable*				CreateParameterClimatology					( double* defaultparams, size_t numpoints, CLIMATOLOGY_HANDLE* species );
		bool										ConfigureDefaultParameterProfiles			();
													skOpticalProperties_MartHybridProfile		( const skOpticalProperties_MartHybridProfile& other );		// Dont allow copy constructor
		skOpticalProperties_MartHybridProfile&		operator =									( const skOpticalProperties_MartHybridProfile& other );		// Dont allow assignment operator
//		bool										DeepCopy									( const skOpticalProperties_MartHybridProfile& other );

	public:
													skOpticalProperties_MartHybridProfile		();
		virtual									   ~skOpticalProperties_MartHybridProfile		();
		void										SetTropopauseHeight							( double tropopauseHeightKm );
		void										SetUTLSThickness							( double UTLSThicknessKm );
		skOpticalProperties_AerosolProfileH2SO4*	SulphateAerosol								() { return m_sulphateaerosol;}
		skOpticalProperties_AerosolProfileIce*		IceCirrus									() { return m_icecirrus;}				//!< calculates scattering properties of ice cirrus. User can modify effective radius and gamma rate.
		bool										SetAerosolModeRadiusAndWidthProfile			( const double* altmeters, const double* moderadius_microns, const double* modewidth, size_t numalt);
		bool										SetIceModeRadiusAndWidthProfile				( const double* altmeters, const double* moderadius_microns, const double* modewidth, size_t numalt);
		bool										ExtinctionPerParticle						( const GEODETIC_INSTANT& geopt, double wavelen_nm, double* extinctionperparticle );

	public:
		virtual bool								InternalClimatology_UpdateCache					( const GEODETIC_INSTANT& pt)  override;
		virtual bool								IsScatterer									() const  override;
		virtual bool								IsAbsorber									() const  override;
		virtual bool								SetAtmosphericState							( skClimatology* neutralatmosphere)  override;
		virtual bool								SetLocation									( const GEODETIC_INSTANT& pt, bool* crosssectionschanged )  override;
		virtual bool								CalculateCrossSections						( double wavenumber, double* absxs, double* extxs, double* scattxs )  override;
		virtual bool								CalculatePhaseMatrix						( double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix  )  override;

	public:
		static bool									ConvertV507AerosolNumberDensityToExtinction	( std::vector<double>* profile );
		static bool									ConvertV507AerosolExtinctionToNumberDensity	( std::vector<double>* profile );

//		static bool 								ConvertV507AerosolExtinctionToNumberDensity( const std::vector<double>&		gridAlts,
//																								 const std::vector<double>&		extinctionperkm,
//																								 std::vector<double>*			numberdensity,
//																								 double							wavelen_nm,
//																								 double							ignorevalue);
};
