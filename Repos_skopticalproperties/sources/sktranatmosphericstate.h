
class skOpticalProperties;
class skRTPhaseMatrix;

/*-----------------------------------------------------------------------------
 *					class SKTRAN_AtmosphericOpticalStateEntry_V2		2008-2-26*/
/** \internal **/
/*---------------------------------------------------------------------------*/

class SKTRAN_AtmosphericOpticalStateEntry_V21
{
	private:
		CLIMATOLOGY_HANDLE				m_species;							//!< The species identification used to get number density from climatology
		skClimatology*					m_climatology;						//!< The climatology of the species. Provides number density per cm3
		skOpticalProperties*			m_particleprops;					//!< The optical properties of "one" particle
		double							m_numberdensity;
		double							m_absxs;							// !< absorption cross-section from ParticleOpticalProps
		double							m_extxs;
		double							m_scattxs;

	private:
		void							ReleaseResources			();

//		bool							DeepCopy								( const SKTRAN_AtmosphericOpticalStateEntry_V21& other ); 
		SKTRAN_AtmosphericOpticalStateEntry_V21& operator =						( const SKTRAN_AtmosphericOpticalStateEntry_V21& other );		// Dotn allow assignment operator

	public:
		bool							operator ==					( const SKTRAN_AtmosphericOpticalStateEntry_V21& other ) const {return (other.m_species == m_species) == TRUE;}

	public:
										SKTRAN_AtmosphericOpticalStateEntry_V21	();
										SKTRAN_AtmosphericOpticalStateEntry_V21	( const CLIMATOLOGY_HANDLE& m_species);
										SKTRAN_AtmosphericOpticalStateEntry_V21 ( const SKTRAN_AtmosphericOpticalStateEntry_V21& other );		// Only allow copy constructor when copying blank objects
									   ~SKTRAN_AtmosphericOpticalStateEntry_V21	();
		double							AbsorptionCrossSection					() const { return m_absxs;}
		double							ExtinctionCrossSection					() const { return m_extxs;}
		double							ScatteringCrossSection					() const { return m_scattxs;}

		skOpticalProperties*			ParticleOpticalProps				() { return m_particleprops;}
		bool							Configure							( CLIMATOLOGY_HANDLE species, skClimatology* numberdensityclimatology, skOpticalProperties* particleopticalprops);
		bool							CalculateCrossSections				( double wavenumber, skClimatology* neutralatmosphere, const GEODETIC_INSTANT& placeandtime );
		bool							CalculateMultiWaveExtinctionsPerCM	( const std::vector<double>&	wavenumber, skClimatology* neutralatmosphere, const GEODETIC_INSTANT& placeandtime, std::vector<double>* absxs, std::vector<double>* extxs,  std::vector<double>* scattxs );
		bool							UpdateNumberDensityPerCM3			( const GEODETIC_INSTANT& placeandtime, bool updatecache );
		double 							CurrentNumberDensityPerCM3			( ) { return m_numberdensity;}
		bool							UpdateClimatology					( skClimatology* numberdensityclimatology );
		skClimatology*					GetClimatology						( ) { return m_climatology; }
		CLIMATOLOGY_HANDLE				GetSpecies							() { return m_species; }
};

/*---------------------------------------------------------------------------
 *					class SKTRAN_AtmosphericOpticalState_V2					2003-12-5*/
/** \ingroup atmosState
 *	A class used by end users to specify the optical properties of an atmosphere. This
 *	class is normally used by radiative transfer models to calculate the total
 *	extinction, scattering, absorption  and scattering phase matrices at any point and time
 *	in the atmosphere.
 *
 *	The class allows the user to pick and add an unlimited set of components that define
 *	the complete optical properties of the atmosphere. Each entry consists of a unique label, a
 *	"number density" climatology and the optical properties of one "particle"  The unique label is
 *	implemented as a CLIMATOLOGY_HANDLE ( see www.guidgen.com to make your own). We provide some standard
 *	names for convenience.
 *
 *	The "number density" is implemented as classes derived from skClimatology. The
 *	climatology classes explicitly implement the concept of caching to avoid excessive model updating
 *	during radiative transfer calculations. The skClimatology must support the species identified by the unique label.
 *	Only one entry of the unique label can be in the set of components stored in this class.
 *
 *	The optical properties are implemented as classes derived from skOpticalProperties. The skOpticalProperties classes
 *	are allowed to have an optional internal climatology to provide exotic parameters like mode radius and mode width.
 *	The standard atmospheric state parameters of pressure and temperature are provided through a separate, explicit
 *	atmospheric state climatology.
 *	
 *  We show an example below where user creates 4 entries, modifies the entries to meet his needs and isnerts them into
 *	the SKTRAN_AtmosphericOpticalState_V21 object. Note that all of the skClimatologies and skOpticalProperties must
 *	be created on the heap (with "new") as the object lifetimes are managed through the nxUnknown interface.
 *
 	\code
 	skOpticalProperties_NO2_OSIRISRes*			m_optno2;			
	skOpticalProperties_O3_OSIRISRes*			m_opto3;	
	skOpticalProperties_MartHybridProfile*		m_optaer;		
	skOpticalProperties_RayleighDryAir*			m_optair;

	skClimatology*								m_clmstrataer;
	skClimatology_Pratmo*						m_clmno2;
	skClimatology_LabowOzoneVMR*				m_clmo3;
	skClimatology*								m_clmair;


	m_clmo3       = new skClimatology_LabowOzoneVMR;
	m_clmno2      = new skClimatology_Pratmo;
	m_clmstrataer = CreateOSIRISAerosolAprioriClimatology_V507();
	ok            = (::CreateOsirisEcmwfClimatologyInstance( &m_clmair ) == S_OK);	

	m_optno2 = new skOpticalProperties_NO2_OSIRISRes;				
	m_opto3  = new skOpticalProperties_O3_OSIRISRes;	
	m_optair = new skOpticalProperties_RayleighDryAir;
	m_optaer = new skOpticalProperties_MartHybridProfile;			
	

	// settings for MARTHybridAerosol optical properties
	m_optaer->SulphateAerosol()->SetLogNormalProfileClimatology( g_sizeparamalts_data, g_moderadius_sulphate_data, g_modewidth_sulphate_data, 2 ); 
	m_optaer->IceCirrus()->SetGammaProfileClimatology          ( g_sizeparamalts_data, g_moderadius_ice_data,      g_modewidth_ice_data,      2 );	
	m_optaer->SetTropopauseHeight( -9999 );						//only sulphate aerosol is used
	m_optaer->SetUTLSThickness( 1.0 );							//this is in km!
	
	species->erase();
	ok = ok && species->AddSpecies( SKCLIMATOLOGY_AEROSOL_CM3,          m_clmstrataer, m_optaer);
	ok = ok && species->AddSpecies( SKCLIMATOLOGY_NO2_CM3,              m_clmno2,      m_optno2);
	ok = ok && species->AddSpecies( SKCLIMATOLOGY_O3_CM3,               m_clmo3,       m_opto3);
	ok = ok && species->AddSpecies( SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, m_clmair,      m_optair);
	ok = ok && species->SetAtmosphericStateModel( m_clmair );  											// set ecmwf as the climatology for the temperature profile 

	\endcode
 */
/*------------------------------------------------------------------------- */

class SKTRAN_AtmosphericOpticalState_V21
{
	protected:
		SKTRAN_AtmosphericEmission									m_atmosphericemission;					// Embed the emission table within the optical property table
		skClimatology*           									m_defaultatmosphericstate;				//!< The default atmospheric state used for the pressure and temperature of cross-sections
		skClimatology*												m_atmosphericstate;						//!< The atmospheric state used for pressure and temperature for cross-sections
		skBRDF*														m_albedo;		
		SKTRAN_BRDF_Lambertian*										m_constantalbedo;
		GEODETIC_INSTANT											m_placeandtime;
		bool														m_updateclimatologycache;
		bool														m_isdirty;
		double														m_wavenumber;
		double														m_kabs;
		double														m_kext;
		double														m_kscat;
		double                                                      m_kdelta;  // sum{ (Eddington delta forward scatterer cross section) * (number density) }

		    std::list<SKTRAN_AtmosphericOpticalStateEntry_V21>					m_species;
	typedef std::list<SKTRAN_AtmosphericOpticalStateEntry_V21>::iterator		iterator;
	typedef std::list<SKTRAN_AtmosphericOpticalStateEntry_V21>::const_iterator	const_iterator;

	protected:
		void							ReleaseResources					();
		void							SetDirty							()		{ m_isdirty = true; }
		void							SetPendingCacheUpdate				()		{ m_updateclimatologycache = true;}
		bool							CheckDirtyAndUpdate					();
		bool							CalculateCrossSections				();
		void							CheckCosineRange					(double * cosscatteringangle);
		bool							TimeAndPlaceIsValid					() const {return (m_placeandtime.mjd > 0.0) && (m_placeandtime.latitude >= -90.1);} 
		bool							CheckClimatologyCacheIsValid		(bool warnaboutbadtime);


	public:
										SKTRAN_AtmosphericOpticalState_V21	();
		virtual						   ~SKTRAN_AtmosphericOpticalState_V21	();

	// ----- Methods generally used by SASKTRAN internals
	public:
		SKTRAN_AtmosphericEmission*		EmissionObjectVar			() { return &m_atmosphericemission;} 		
		bool							VectorPhaseMatrix			( double cosscatteringangle, skRTPhaseMatrix* P);
		bool							ScalarPhaseMatrix           ( std::pair<double, size_t> cosscatterandindex, double& p11 );
		bool							ScalarScatteringCoefficient	( double cosscatteringangle, double*          P);
		double							ExtinctionPercm				()		{ CheckDirtyAndUpdate(); return  m_kext;}				//!< Get the total extinction per cm from all of the species at the current model location */
		double							AbsorptionPercm				()		{ CheckDirtyAndUpdate(); return  m_kabs;}				//!< Get the total absorption per cm from all of the species at the current model location */
		double							ScatteringPercm				()		{ CheckDirtyAndUpdate(); return  m_kscat;}				//!< Get the total non-delta forward scattering per cm from all of the species at the current model location */
		double                          DeltaForwardPercm           ()      { CheckDirtyAndUpdate(); return  m_kdelta;}             //!< Get the total delta forward scattering per cm from all the species at the current model location */
		bool							GetAlbedoObject				(skBRDF** albedo);
		bool							SetWavelength				( double wavelen_nm );
		GEODETIC_INSTANT				GetTimeAndLocation			()		{ return m_placeandtime; }
//		bool							DeepCopy					( const SKTRAN_AtmosphericOpticalState_V21& other );

	// ----- Methods generally called by users
	public:
		void							erase						( )						{ReleaseResources();}
		bool							SetAlbedoObject				(skBRDF*	albedo);
		bool							SetAlbedo					(double      albedo);
		bool							SetAtmosphericStateModel	( skClimatology* atmosphericstate);
		bool							SetTimeAndLocation			( const GEODETIC_INSTANT& point, bool updateclimatologycache);	//!< Update the cross sections and internal phase matrix lookup table
		bool							AddEmission					( const CLIMATOLOGY_HANDLE&  species, skEmission* emissionobject);	
		bool							AddSpecies					( const CLIMATOLOGY_HANDLE&  species, skClimatology* numberdensityclimatology, skOpticalProperties* particleopticalprops);
		bool							RemoveSpecies				( const CLIMATOLOGY_HANDLE& species);
		bool							UpdateSpeciesClimatology	( const CLIMATOLOGY_HANDLE& speciesinlist, skClimatology* numberdensityclimatology ); 
		bool							GetSpeciesClimatology		( const CLIMATOLOGY_HANDLE& speciesinlist, skClimatology**       numberdensityclimatology );
		bool							GetSpeciesOpticalProperties	( const CLIMATOLOGY_HANDLE& speciesinlist, skOpticalProperties** opticalprops);
		bool							GetAtmosphericStateModel	( skClimatology** statemodel );
		bool							CalculateMultiWaveCrossSections (const std::vector<double>& wavenumber, std::vector<double>* absxs, std::vector<double>* extxs, std::vector<double>* scattxs );
		bool							CalculateMultiWaveCrossSectionsAndPhaseMatrix(const std::vector<double>& wavenumber,
																					   std::vector<double>* kabs,
																					   std::vector<double>* kext,
																					   std::vector<double>* kscat,
																					   const std::vector<double>& cosangles,
																					   nx2dArray<skRTPhaseMatrix>* P );
		bool							PhaseGridHint(const std::vector<double>& cosscatterangles);
};
