
#if defined(_MSC_VER)
#define BD_TIPS_2017 BD_TIPS_2017
#else
#define BD_TIPS_2017 bd_tips_2017_
#endif

extern "C" void BD_TIPS_2017(int* MOL, double* temp, int* ISO, double* gi, double* QT);


/*-----------------------------------------------------------------------------
 *					skHitranPartitionTableEntry		2013-3-12*/
/** \ingroup hitranoptpropinternals
  *	A small class to hold all of the HITRAN molecule info which is stored
  * Gobal_Data/molparam.txt and Global_Data/parsum.dat. This info is required when calculating
  *	line strength and line shape. One instance of this class is created
  *	for each molecule/isotope defined in molparam.txt. 
  *
  *	The entry stores an internal partiton table in members m_T and m_Q.
  *	The internal partition data is read in from Global_Data/parsum.dat and
  *	is available for most of the molecule isotopologues defined in molparam.txt
  *	but not for all of the more esoteric isotopes.
  *
  *	Note 1) that the isotope CO2 827 in HITRAN 2008 parsum.dat is defined
  *	as CO2 728 in molparam.txtx. I edited my version of molparam.txt so
  *	it also read CO2 827.
  *
  *	Note 2) The file parsum.dat had isotope information for seveal isotopes
  *	of O3 that are not present in molparam.txt. I ignored these extra O3
  *	isotopologues.
**/
/*---------------------------------------------------------------------------*/

class skHitranPartitionTableEntry
{
	public:
		size_t					m_moleculeid;		// The Hitran molecule number (typically 1-41) , eg 1 is H2O
		size_t					m_isotopeid;		// The isotope id. this is typically 3 or 4 digits like 121 but can be up to 4 digits. we allow up to 7
		size_t					m_isotopeorder;		// This is the isotope number (e.g. 1,2,3,4,5) used in the spectral line par files. It follows isotope abundance;
		std::string				m_chemicalname;		// The checmical name as written in HITRAN file molparam.txt
		double					m_abundance;		// The natural abundance as given in HITRAN file molparam.txt 
		double					m_gj;				// The gj as given in HITRAN file molparam.txt 
		double					m_molarmass;		// The molar mass as givn in HITRAN file molparam.txt
		bool					m_hapicompliant;		// HAPI compliant HITRan calculates the Internal Partition Sumns
		std::vector<double>		m_T;				// Temperature in kelvins for Partition table data 70K to 3000K, from HITRAN file parsum.dat
		std::vector<double>		m_Q;				// Partition data for this isostope, same size as m_T, from HITRAN file parsum.dat. Not available for all isotopes of all molecules

	public:
		size_t					ToKey				();
		bool					ReserveSpace		(size_t numtemperature);
		bool					AddEntry			( double t, double q);
		double					InternalPartition	( double T) const;
		size_t					MoleculeNumber		() const { return m_moleculeid;}
		size_t				    IsotopeNumber		() const { return m_isotopeid;}

};



/*-----------------------------------------------------------------------------
 *					skHitranMoleculeManager		2013-3-12*/
/**  \ingroup hitranoptpropinternals
 *	This is a class that keeps a list of all the available molecular
 *	isotopologues available in the hitran version stored on this users computer. 
 *	I assume the user has copied the (2008 or later) file structure from the Hitran FTP site.
 *	I use a registry setting (on both Linux and Windows) to save the base directory
 *	of the HITRAN information. We used to make intermediate binary files of HITRAN text files
 *	but I find modern computers are so blazingly fast that that is no longer necessary.
 *	Consequently I just decode the original HITRAN text files on the fly. Makes
 *	maintenance of the code and databases much easier.
 *
 *	This class is meant to be created as one instance shared between several instances
 *	of skSpectralLineCollection_HitranChemical
 *	
**/
/*---------------------------------------------------------------------------*/

class skHitranMoleculeManager : public nxUnknown
{
	private:
				std::map<size_t, skHitranPartitionTableEntry>					m_molecules;
		typedef std::map<size_t, skHitranPartitionTableEntry>::iterator			iterator;
		typedef std::map<size_t, skHitranPartitionTableEntry>::const_iterator	const_iterator;
	public:
		typedef std::map<size_t, skHitranPartitionTableEntry>::value_type	value_type;

	private:
		bool											m_hapicompliant;

	private:
		bool											FindHitranMolparamFile				( nxString* filename);
		bool											FindHitranParsumFile				( nxString* filename);
		bool											FindHitranGlobalFileFile			( const char* filename, nxString* fullfilename );
		static bool										LoadBaseDirectoryNameFromRegistry	( nxString* filename );
		bool											LoadMoleculeDefinitions				( );
		bool											LoadPartitionDefinitions			( );
		bool											FindMoleculeEntry					( const char* chemicalname, size_t isotopeid, skHitranPartitionTableEntry** entry );

	public:
														skHitranMoleculeManager				(bool hapicompliant);
		virtual										   ~skHitranMoleculeManager				();
		bool											FindHitranMoleculeDirectory			( nxString* moleculedir) const;
		static const skHitranMoleculeManager*			CreateManagerInstance				( bool hapicompliant);
		bool											FindMoleculeId						( const char* chemicalname, size_t*	moleculeid) const;
		bool											FindMoleculeEntry					( size_t moleculeid, size_t isotopeid, const skHitranPartitionTableEntry** entry ) const;
		bool											FetchAllIsotopeEntries				( size_t moleculeid, std::vector<const skHitranPartitionTableEntry*> *table ) const;


};


/*---------------------------------------------------------------------------
 *                skSpectralLine_HitranLineStruct                 2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

struct HitranLineStruct
{
	int										m_molNum;				// molecule number										unitless
	int										m_isoNum;				// isotope order number. This indexes (1 based) the isotope order as written in molparam.txt
	double									m_nuTrans;				// transition wavenumber								inverse centimetres
	double									m_intensity;			// intensity (integrated line strength)					centimetres
	double									m_einsteinA;			// Einstein-A coefficient								inverse seconds
	double									m_airBroad;				// air-broadened half-width at half max					inverse (cm * atm)
	double									m_selfBroad;			// self-broadened half-width at half max				inverse (cm * atm)
	double									m_nuLow;				// lower-state energy									inverse centimetres
	double									m_nAir;					// temperature dependence exponent of airBroad			unitless 
	double									m_deltaAir;				// air-pressure induced shift							inverse (cm * atm)
	double									m_upperStatWt;			// upper-state statistical weight						unitless
	double									m_lowerStatWt;			// lower-state statistical weight						unitless
	char									m_qGlobU[17];			// upper-state "global" quanta							unitless
	char									m_qGlobL[17];			// lower-state "global" quanta							unitless
	char									m_qLocU[17];			// upper-state "local" quanta							unitless
	char									m_qLocL[17];			// lower-state "local" quanta							unitless
	char									m_errCodes[8];			// uncertainty indicies									unitless
	char									m_refCodes[14];			// reference indices									unitless
	char									m_lineFlag[3];			// flag for line-mixing									unitless
};


/*-----------------------------------------------------------------------------
 *					nxSpectralLine_HitranLineStruct		2013-3-8*/
/** \ingroup hitranoptpropinternals
 *	Extracts data form the standard 160-char lines in hitran database. Note
 *	this line format only started in 2002 and onwards. It will not work with
 *	earlier versions of Hitran..
**/
/*---------------------------------------------------------------------------*/

class skSpectralLine_HitranLine : public skSpectralLine		
{
	private:
		const skHitranPartitionTableEntry*		m_moleculeinfo;			// Molecule information that comes from skHitranMoleculeManager
		HitranLineStruct						m_entry;

	private:
		void									ClearRecord						( );
		int										IntegerValFromString			( const char* str, int firstIdx,  int subStrLen );
		double									DoubleValFromString				( const char* str, int firstIdx,  int subStrLen );
		int										ExtendedHexValFromChar			( const char* str, int firstIdx );
		void									SubstringFromString				( char* outputstr, const char* str, int firstIdx,  int subStrLen );

	public:
		bool									Parse160CharRecord				( const char* record, size_t numchar );
		size_t									MoleculeNumber					() const	{ return (size_t)m_entry.m_molNum;}
		size_t									IsotopeOrderNumber				() const	{ return (size_t)m_entry.m_isoNum;}
		bool									SetMoleculeInfo					( const skHitranPartitionTableEntry* info) { m_moleculeinfo = info; return true;}
		const HitranLineStruct&					HitranEntry						() const { return m_entry; }


	public:
												skSpectralLine_HitranLine		( );
												skSpectralLine_HitranLine		( const HitranLineStruct& record);
												skSpectralLine_HitranLine		( const skSpectralLine_HitranLine& other);
		virtual								   ~skSpectralLine_HitranLine		( );
		virtual double							Snm								( )	const	{return m_entry.m_intensity;}		// Spectral Line Intensity from level n to m, at reference temperature. Same value as Hitran database [cm-1/(molecule cm-2)].
		virtual double							Nu								( )	const	{return m_entry.m_nuTrans;}			// The spectral line transition frequency [cm-1] from level n to m (in vacuum)  
		virtual double							GammaAir						( )	const	{return m_entry.m_airBroad;}		// Air  Broadened half width HWHM cm-1/atm at reference temperature 
		virtual double							GammaSelf						( )	const	{return m_entry.m_selfBroad;}		// Self Broadened half width HWHM cm-1/atm at reference temperature
		virtual double							EinsteinA						( )	const	{return m_entry.m_einsteinA;}		// Einstein A coefficient 
		virtual double							ELower							( )	const	{return m_entry.m_nuLow;}			// Lower State energy in cm-1
		virtual double							Nair							( )	const	{return m_entry.m_nAir;}			// Coefficient of temperature dependence of air broadened half width.
		virtual double							Deltaair						( )	const	{return m_entry.m_deltaAir;}		// Air Broadened pressure shift of line transition  in cm-1/atm at reference temeprature
		virtual double							Tref							( )	const	{return 296.0;}						// Reference Temperature in K for the spectral line parameters (often 296K).
};

class skSpectralLineCollection_HitranChemical;

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranIsotope		2013-3-8*/
/** \ingroup spectralline
 *	A class that holds all the lines for one isotope of one molecule from the Hitran database.
 **/
/*---------------------------------------------------------------------------*/

class skSpectralLineCollection_HitranIsotope : public skSpectralLineCollection
{
	private:
		size_t										m_molNum;				// HITRAN molecule number  unitless
		size_t										m_isoNum;				// HITRAN isotopologue number unitless
		const skHitranPartitionTableEntry*			m_moleculeinfo;			// Molecule information that comes from skHitranMoleculeManager
		const skSpectralLineCollection_HitranChemical*	m_parentchemical;
		mutable HitranPartitionTableCache			m_partitioncache;


	public:
		static	size_t								g_numinstances;

	public:
													skSpectralLineCollection_HitranIsotope	( const skSpectralLineCollection_HitranChemical* parentchemical, size_t molNum, size_t isoNum, const skHitranPartitionTableEntry* moleculeinfo);
													skSpectralLineCollection_HitranIsotope	( const skSpectralLineCollection_HitranIsotope& other);
		virtual									   ~skSpectralLineCollection_HitranIsotope	();
		const skHitranPartitionTableEntry*			MoleculeInfo							() const	{ return m_moleculeinfo;}			// Molecule information that comes from skHitranMoleculeManager
		const skSpectralLineCollection_HitranChemical*	ParentChemical						() const	{ return m_parentchemical;}
		size_t										MolNum									() const	{ return m_molNum;}
		size_t										IsoNum									() const	{ return m_isoNum;}

	public:
		virtual double								QPartition			( double T ) const override;
		virtual double								MassAMU				( ) const override;
		virtual double								PartialPressure		( const GEODETIC_INSTANT& geopt, double airpressure, double airtemp) const override;
};

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_IsotopeCollectionHitran		2013-3-8*/
/** **/
/*---------------------------------------------------------------------------*/

class skSpectralLineCollection_HitranChemical
{
	private:
				std::map< size_t, skSpectralLineCollection_HitranIsotope>					m_isotopes;			// The array of isotopes associated with this chemical
		typedef std::map< size_t, skSpectralLineCollection_HitranIsotope>::value_type		value_type;
		typedef std::map< size_t, skSpectralLineCollection_HitranIsotope>::iterator			iterator;
		typedef std::map< size_t, skSpectralLineCollection_HitranIsotope>::const_iterator	const_iterator;

	private:
		CLIMATOLOGY_HANDLE							m_numberdensityguid;				// The number density guid for retrieving number density of this chemical species from the self-broadening climatology for partial pressure
		skClimatology*								m_selfbroadeningclimatology;		// An optional climatology used to fetch the number density of the chemical species used in calculating partial pressure for self broadening
		const skHitranMoleculeManager*				m_manager;							// The Hitran molecule database manager, stores global hitran data.
		size_t										m_molnum;							//
		std::string									m_chemicalname;
		size_t										m_isotopeidfilter;					// Used to select just one isotope (default is 0 which is all isotopes);
		std::string									m_upperstate_globalquanta_filter;	// Filter used to select only lines that belong to the given upper state global quanta
		std::string									m_lowerstate_globalquanta_filter;	// Filter used to select only lines that belong to the given lower state global quanta
		double										m_lowerwavenumber;					// The lower wavenumber of the current micro-window. Used for selecting lines from the HITRAN/spectral line database
		double										m_upperwavenumber;					// The upper wavenumber of the current micro-window. Used for selecting lines from the HITRAN/spectral line database
		double										m_windowmargin;						// The number of wavenumbers to add to the edges of each micro-window
		double										m_maxlinestrength;					// The maximum Snm line strength in the window
		double										m_userdefined_maxlinestrength;		// A User defined maximum line strength
		bool										m_use_userdefined_maxlinestrength;
		bool										m_hapicompliant;

	public:
		static size_t								g_numinstances;


	private:
		bool										FindFile							( nxString* filename);
		bool										LoadFile							();
		bool										FindIsotopeId						( size_t isotopeid, const std::vector<const skHitranPartitionTableEntry*>& isotopetable);
		bool										UpdateMaxLineStrength				();
		void										ReleaseResources					();
		bool										InsertSpectralLineEntryNoFilter		( size_t isotopeid, skSpectralLine_HitranLine* spectralline);
		bool										SetChemicalName						( const char* chemicalname );
		bool										SetLineLimitsfromMaxLineStrength	();
		std::string									GlobalStateToString					( const char* globalstate) const;
		bool										MatchesGlobalQuanta					( const HitranLineStruct* spectralline ) const;

	public:
													skSpectralLineCollection_HitranChemical	( const char* chemicalname, double lowerwavenumber, double upperwavenumber, double windowmargin_wavenumbers, bool hapicompliant, int isotopefilter, const char* lowerstate_globalquanta_filter, const char* upperstate_globalquanta_filter);
												   ~skSpectralLineCollection_HitranChemical	();
		bool										SetIsotopeFilter						( size_t filterid);
		bool										SetLineShapeObject						( skSpectralLineShape* lineshapeobject );
		bool										SetUserDefinedMaxLineStrength			( double maxlinestrength);
		bool										SetLineTolerance						( double tolerance );
		bool										UpdateLocation							( const GEODETIC_INSTANT& geopt, skClimatology* atmosphere );
		bool										AbsorptionCrossSection					( double nu, double* absxsec  ) const;
		bool										AddAbsorptionCrossSectionArray			( const std::vector<double>& wavenum, std::vector<double>* absxs);
		bool										SetSelfBroadeningClimatology			( const CLIMATOLOGY_HANDLE& parameterguid, skClimatology* numberdensityclimatology );
		bool										SelfBroadeningClimatology_UpdateCache	( const GEODETIC_INSTANT& geopt );
		double										PartialPressure							( const GEODETIC_INSTANT& geopt, double airpressure, double airtemp) const;
		bool										SetNumThreads							( size_t numthreads);
		const std::map< size_t, skSpectralLineCollection_HitranIsotope>& Isotopes			() const { return m_isotopes;}
		int											MoleculeNumber							() const { return (int)m_molnum;}
		const skHitranMoleculeManager*				MoleculeManager							() const { return m_manager;}					// The Hitran molecule database manager, stores global hitran data.
		bool										InsertAllSpectralLines					( size_t isotopeid, const std::vector<HitranLineStruct>& spectrallines);
		double										MicroWindowMinWavenum					() const { return m_lowerwavenumber; }
		double										MicroWindowMaxWavenum					() const { return m_upperwavenumber; }
};



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical		2013-3-20*/
/** \ingroup hitranoptprop
 *	The skOpticalProperties interface to the HITRAN database. The code treats
 *	the HITRAN database as a collection of about 41 different chemicals. Each chemical
 *	is a collection of various isotopes and each isotope is a collection of
 *	spectral lines.
 *	
 *	This class wraps up all of the detail and provides a flexible interface
 *	that should provide the capability for state-of-the-art spectral line
 *	calculations. The class allows the user to programatically select the
 *	line shape calculation of their choice.
 *
 *	Atmospheric state information (typically pressure and temperature) is
 *	enabled through the standard climatology interfaces in skOpticalProperties
 *	(i.e. SetAtmopshericState) while partial pressure climatologies required
 *	for self broadening calculations are supported through SetNumberDensityClimatology
 *	and InternalClimatology_UpdateCache.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_HitranChemical : public skOpticalProperties,
	                                       public skOpticalProperty_AdditionalStateInfo_PressTemperatureDependent
{
	private:
		skClimatology*									m_atmosphericstateclimatology;	// The climatology used to fetch temperature and pressure 
		skSpectralLineCollection_HitranChemical*		m_hitranchemical;				// The class object that actually stores all of the Hitran isotopes and spectral lines. Not created until needed
		nxString										m_chemicalname;					// The chemical name of this molecule, used when constructing the instance of m_hitranchemical
		Hitran_CrossSection_Cache*						m_xs_optimizer;					// An object to optimize cross-sections calculation for engines such as HR which cannot call CalculateCrossSectionsArray
		skSpectralLineShape*							m_lineshapeobject;				// The line shape object. Typically Voigt-Kuntz
		skClimatology*									m_selfbroadeningclimatology;	// The climatology used to fetch the number density of this molecule. It is used to calculate the partial pressure required for self-broadening
		CLIMATOLOGY_HANDLE								m_chemicalnumberdensityguid;	// The handle used by the self-broadening climatology to fetch the number density (molecules/cm3) of this molecule
		GEODETIC_INSTANT								m_lastpoint;
		bool											m_isdirty;
		std::string										m_lowerstate_global_quanta_filter;
		std::string										m_upperstate_global_quanta_filter;
		int												m_isotopefilter;
		bool											m_manual_microwindow_isset;
		double											m_manualtolerance;
		double											m_manualmaxlinestrength;
		double											m_lowwavenum;
		double											m_hihwavenum;
		double											m_margin_width_wavenum;			// The number of wavenumbers to add to the edges of the micro-window
		size_t											m_numthreads;
		bool											m_hapicompliant;				// uses files generated from the HAPI compliant interfaces

	private:
		void								SetDirty							();
		bool								CheckDirtyAndUpdate					();
		void								init								();
		bool								CalculateCrossSectionsInternal		( double wavenumber, double* absxs, double* extxs, double* scattxs ) const;
		bool								CheckWavenumberIsAscending			( const std::vector<double>&	wavenumber) const;

	private:
											skOpticalProperties_HitranChemical	( const skOpticalProperties_HitranChemical& other );			// Dont allow copy constructor
		skOpticalProperties_HitranChemical& operator =							( const skOpticalProperties_HitranChemical& other );			// Dont allow assignment oeprator

	public:
											skOpticalProperties_HitranChemical	();
											skOpticalProperties_HitranChemical	(const char* chemicalname );
		virtual							   ~skOpticalProperties_HitranChemical	();
		bool								SetChemicalName						( const char*  chemicalname);
		bool								SetWavenumberRange					( double lowwavenumber, double highwavenumber );
		bool								SetMicroWindowMargin				( double margin_wavenum);
		bool								SetSelfBroadeningClimatology		( skClimatology* numberdensityclimatology  );
		bool								SetSelfBroadeningClimatologyHandle	( const CLIMATOLOGY_HANDLE& parameterguid);
		bool								SetLineShapeObject					( skSpectralLineShape* lineshapeobject );
		bool								SetIsotopeFilter					( size_t filterid);
		bool								SetUserDefinedMaxLineStrength		( double maxlinestrength);
		bool								SetLineTolerance					( double tolerance );
		bool								EnableCachedCrossSections			( double* wavenumbers, size_t numwave );
		void								SetHitran2008Compliant				() { m_hapicompliant = false; }
		bool								SetLowerStateGlobalQuantaFilter		( std::string lowerstateglobalquanta);
		bool								SetUpperStateGlobalQuantaFilter		( std::string upperstateglobalquanta);
		skSpectralLineCollection_HitranChemical*	HitranChemical				() { return m_hitranchemical;}
		skClimatology*                      AtmosphericStateClimatology         () { return m_atmosphericstateclimatology; }

	public:
		virtual bool						SetAtmosphericState					( skClimatology* neutralatmosphere) override;
		virtual bool						SetLocation							( const GEODETIC_INSTANT& pt, bool* crosssectionschanged ) override;
		virtual bool						CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs) override;
		virtual bool						CalculateCrossSectionsArray			( const double * wavenumber, int numwavenumber, double *absxs,  double* extxs, double* scattxs) override;
		virtual bool						InternalClimatology_UpdateCache		( const GEODETIC_INSTANT& pt) override;
		virtual bool						IsScatterer							() const override				{ return false;} 
		virtual bool						IsAbsorber							() const override				{ return true;} 
};


/*---------------------------------------------------------------------------
 *          Class skOpticalProperties_HitranChemical2008          2019-05-24 */
/** **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_HitranChemical2008 : public skOpticalProperties_HitranChemical
{
	public:
										skOpticalProperties_HitranChemical2008	();
										skOpticalProperties_HitranChemical2008  (const char* chemicalname, double lowerwavenumber, double upperwavenumber);
	virtual							   ~skOpticalProperties_HitranChemical2008	();


};
