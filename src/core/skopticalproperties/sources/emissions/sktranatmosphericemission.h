class skEmission;
class skSolarSpectrum;

/*-----------------------------------------------------------------------------
 *					class SKTRAN_AtmosphericEmissionEntry		2008-2-26*/
/** \internal **/
/*---------------------------------------------------------------------------*/

class SKTRAN_AtmosphericEmissionEntry
{
	private:
		CLIMATOLOGY_HANDLE					m_species;
		skEmission*							m_emission;					//!< The optical properties of "one" particle
		double								m_radiance;

	private:
		void								ReleaseResources				();
		SKTRAN_AtmosphericEmissionEntry&	operator =						( const SKTRAN_AtmosphericEmissionEntry& other );		// Dotn allow assignment operator

	public:
		bool								operator ==						( const SKTRAN_AtmosphericEmissionEntry& other ) const {return (other.m_species == m_species) == true;}
		skEmission*							EmissionObject					() { return m_emission;}

	public:
											SKTRAN_AtmosphericEmissionEntry			();
											SKTRAN_AtmosphericEmissionEntry			( const CLIMATOLOGY_HANDLE& species);
											SKTRAN_AtmosphericEmissionEntry			( const SKTRAN_AtmosphericEmissionEntry& other );		// Only allow copy constructor when copying blank objects
										   ~SKTRAN_AtmosphericEmissionEntry			();
		double								IsotropicRadiance						() const { return m_radiance;}
		bool								Configure								( CLIMATOLOGY_HANDLE species, skEmission* emission);
		bool								CalculateEmission						( double wavenumber, /* skClimatology* neutralatmosphere, */ const GEODETIC_INSTANT& placeandtime, bool isground );
		bool								CalculateMultiWaveEmission				( const std::vector<double>&	wavenumber, /* skClimatology* neutralatmosphere, */ const GEODETIC_INSTANT& placeandtime, bool isground, std::vector<double>* emission);
		bool								UpdateInternalClimatologies				( const GEODETIC_INSTANT& placeandtime, bool updatecache );
		CLIMATOLOGY_HANDLE					GetSpecies								() { return m_species; }
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
 *	the SKTRAN_AtmosphericEmission object. Note that all of the skClimatologies and skOpticalProperties must
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

class SKTRAN_AtmosphericEmission
{
	private:
//		skClimatology_MSIS90*										m_defaultatmosphericstate;				//!< The default atmospheric state used for the pressure and temperature of cross-sections
		skSolarSpectrum*											m_solarspectrum;
//		skClimatology*												m_atmosphericstate;						//!< The atmospheric state used for pressure and temperature for cross-sections
		GEODETIC_INSTANT											m_placeandtime;							//!< The current place and time.
		bool														m_isground;								//!< Flags thatthe current place and time represents a ground point
		double														m_wavenumber;
		bool														m_updateclimatologycache;
		bool														m_isdirty;
		double														m_radiance;

		    std::list<SKTRAN_AtmosphericEmissionEntry>					m_species;
	typedef std::list<SKTRAN_AtmosphericEmissionEntry>::iterator			iterator;
	typedef std::list<SKTRAN_AtmosphericEmissionEntry>::const_iterator	const_iterator;

	private:
		void							SetDirty							()		{ m_isdirty = true; }
		void							SetPendingCacheUpdate				()		{ m_updateclimatologycache = true;}
		bool							CheckDirtyAndUpdate					();
		bool							CalculateEmissions					();
		bool							TimeAndPlaceIsValid					() const {return (m_placeandtime.mjd > 0.0) && (m_placeandtime.latitude >= -90.1);} 
		bool							CheckClimatologyCacheIsValid		(bool warnaboutbadtime);


	public:
										SKTRAN_AtmosphericEmission	();
		virtual						   ~SKTRAN_AtmosphericEmission	();
		void							ReleaseResources			();

	// ----- Methods generally used by SASKTRAN internals
	public:			
		double							IsotropicRadiance			()				{ CheckDirtyAndUpdate(); return m_radiance;}				//!< Get the total extinction per cm from all of the species at the current model location */
		bool							SetWavelength				( double wavelen_nm );
		GEODETIC_INSTANT				GetTimeAndLocation			() const		{ return m_placeandtime; }

	// ----- Methods generally called by users
	public:
		void							erase						()						{ReleaseResources();}
		bool							SetSolarSpectrum			( skSolarSpectrum* solarspectrum);
		bool							SetTimeAndLocation			( const GEODETIC_INSTANT& point, bool isgroundpoint, bool updateclimatologycache);	//!< Update the cross sections and internal phase matrix lookup table
		bool							AddEmission					( const CLIMATOLOGY_HANDLE&  species, skEmission* emission);
		bool							RemoveEmission				( const CLIMATOLOGY_HANDLE& species);
		bool							GetSpeciesEmissionObject	( const CLIMATOLOGY_HANDLE& speciesinlist, skEmission**			emission);
		bool							CalculateMultiWaveEmissions (const std::vector<double>& wavenumber, std::vector<double>* radiance);
};
