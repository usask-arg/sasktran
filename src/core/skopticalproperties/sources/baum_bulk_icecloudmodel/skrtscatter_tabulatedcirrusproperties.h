
/*-----------------------------------------------------------------------------
 *			skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry		2011-4-18*/
/** \internal 
 *	A class that hold one tabulated phase matrix entry for one wavelength at
 *	one altitude.
 **/
/*---------------------------------------------------------------------------*/
class skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry
{
	private:
		size_t						m_roundedwavelen;				//!< The wavelength in nanometers multipled by 1000 and rounded to nearest integer (used as an exact index)
		SKTRAN_GridDefBase_V2		m_phasematrix;					//!< The scalar phase matrix  array (typically 181 or 361 elements)
		SKTRAN_GridDefBase_V2		m_cosscatterangle;				//!< The cosine of the scattering angle array, must be the same size as the phase matrix array

	private:
		void						ReleaseResources();

	public:
									skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry();
		bool						DeepCopy			( const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry& other );
		bool						GetPhaseMatrix		( double cosangle, double* phasematrix ) const;
		bool						ConfigureEntry		( double wavelennm, const double* phasematrix, const double* cosangle, size_t numpoints);
		size_t						WavelengthIndex		() const { return m_roundedwavelen;}
		static size_t				RoundTheWavelength	( double nm ) { return  size_t(floor(1000*nm + 0.5));}
};

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry		2011-4-16*/
/** \internal
 *	A class that holds an array of phase matrices and cos(scatter angles) for a set of 
 *	wavelengths at any one altitude. The phase matrices for each wavelength at this height
 *	can be read in from a text file. 
 *	\par File format
 *	Format is expressed using line numbers for clarity. In reality arrays can be spread over multiple 
 *	lines. White space and end-of-lines are skipped
 *	Line   1:		numwavelengths		numscatteringangles
 *	Line   2:		An array[numscatteringangles] of cos(scatterangle), can be distributed over several lines.
 *	Line   3:		Wavelen_nm[0]		An Array[numscatteringangle] of phasematrix for wavelength 0  
 *	Line   3:		Wavelen_nm[1]		An Array[numscatteringangle] of phasematrix for wavelength 1  
 *	Line   n+2:		Wavelen_nm[n-1]		An Array[numscatteringangle] of phasematrix for wavelength n-1  
 **/
/*---------------------------------------------------------------------------*/
class skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry
{
	private:
		double														m_heightmeters;					//!< The height assocaited with this set of wavelength entries
		size_t														m_numwavelens;					//!< The number of wavelength entries at this altitude
		size_t														m_maxwavelens;					//!< The maximum number of entries allowed
		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry*		m_wavelenentries;				//!< the array of weavelength entries.

	private:
		void														ReleaseResources();
		bool														AllocateEntries	( size_t maxentries );

	public:
																	skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry	();
																   ~skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry	();
																	skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry	( double heightm ); // reserved for std::comparator functions (like lower_bound)
		bool														operator <										( const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry & other ) const { return m_heightmeters < other.m_heightmeters;}
		double														Height											() const	{ return m_heightmeters;}
		bool														LoadWavelengthEntriesFromHeightFile				( double heightm, const char* filename );
		bool														DeepCopy										( const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry& other );
		bool														FindWavelengthEntry								( double wavelennm, const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry** entry ) const;
};


/*-----------------------------------------------------------------------------
 *				skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB	2011-4-11 */
/**	\ingroup aerosolskopticalprop
 *	A class that interfaces to Bryan Baum's database of cirrus ice crystal 
 *	cross-sections and phase functions, built using capabilities of both 
 *	skOpticalProperties_TabulatedExtinction_HeightWavelength and 
 *	skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength.
 *	Required data are retrieved from the database (whose directory is sepcified 
 *	with a registry key) and are cached locally (in the directory specified in 
 *	m_configdir) using the base class functions.
 *	Before you use this class, set up the directory that will be pointed to with
 *	the registry key.  This directory must contain all of the files containing
 *	cross-sections (De=10..180um) and phase functions (lambda=0.4..2.2um). 
 *	Also see Ice Cloud Research Team web site:
 *		http://www.ssec.wisc.edu/~baum/Cirrus/IceCloudModels.html
 *		and in particular http://www.ssec.wisc.edu/~baum/Cirrus/Solar_Spectral_Models.html
 *		and look at section VIS-SWIR Spectral Models:

 *        University of Wisconsin-Madison Space Science and Engineering Center
 *	Ice Cloud Research Team:  
 *		Bryan Baum (bryan.baum@ssec.wisc.edu), 
 *		Ping Yang (Texas A&M), 
 *		Andrew J. Heymsfield and Carl Schmitt (NCAR)     
 *	The cross-section database contains 18 text files, each one specific to a given 
 *	effective diameter (De), which span a range from 10 to 180 microns in 
 *	10 micron increments. Each file contains entries for each of the 144 
 *	wavelengths at which properties are available. The ice crystal habit 
 *	mixture described above is used for all calculations of microphysical and 
 *	optical properties.  
 *	Cross-sections data provided for each wavelength are the following:
 *		Asymmetry factor (g)					..
 *		Single scattering albedo (omega)		..
 *		Extinction efficiency (Qe)				..
 *		Scattering cross section (sig_scat)		(in units of um^2)
 *		Extinction cross section (sig_ext)		(in units of um^2)
 *	Within each of the (ascii) phase function files, data are formatted as follows:
 *		1st column: scattering angle (from 0 to 180 degrees with 498 angles)
 *		2nd column:  mean value of scattering phase function at	De = 10 microns
 *		3rd column:    standard deviation of phase function at	De = 10 microns
 *		4th column:  mean value of scattering phase function at	De = 20 microns
 *		5th column:    standard deviation of phase function at	De = 20 microns
 *		....	and so forth for the rest of the De values up to 180 microns 
 * NOTE: the files read in by this class have been stripped of header data by a
 *		Python script.  These files are subject to change, and the Python script 
 *		is available in the current directory.  Run it on the supplied data files
 *		to strip away the text. After this, set up a directory for these files, 
 *		and this class will generate a registry key to point to this directory.
 *		The most up-to-date version of the phase functions for solar wavelengths
 *		can be found on the above website.
 */
/*---------------------------------------------------------------------------*/
class skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB : public skOpticalProperties
{
	private:
				  std::list<skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry>					m_heightentries;		// The cache of height entries (each height has phase matrices for an array of wavelengths)
		typedef	  std::list<skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry>::iterator			iterator;				// Iterates over the colelction
		typedef	  std::list<skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry>::const_iterator	const_iterator;
		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry*	m_lowerheight;				// The lower height entry used to interpolate height entries
		double																	m_f1;						// the interpolating factor for the lower height entry
		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry*	m_upperheight;				// the upper height entry
		double																	m_f2;						// the interpolating factor for the upper height entry.
		size_t							m_idxh1;								// height index (m_heights) for lower-height extinct
		size_t							m_idxh2;								// height index (m_heights) for upper-height extinct
		std::vector<double>				m_deltafunctionforwardscatterfraction;	//!< The fraction of total scatter (0-1) dedicated to a delta function forward scatter peak.

	private:
		bool							m_isdirty;
		nxString						m_configdir;			// local dir for keeping working files (i.e. './inputdata')
		size_t							m_numfilecolumns;		// num columns in DB files
		nx1dArray<double>				m_heights;				// table heights in m
		nx1dArray<double>				m_wavelennm;			// table wavelengths in nm
		nx1dArray<double>				m_effectivesize_db;		// all database effective sizes in um
		nx1dArray<double>				m_wavelennm_db;			// all database wavelengths in nm
		nx1dArray<double>				m_scatterangles_db;		// all database scattering angles in deg
		nx2dArray<double>				m_crossection;			// [numheights x numlambdas] cross section in cm^2
		nx2dArray<double>				m_deltafraction;		// fraction of incoming photons scattered into forward peak (per height, per wavelength)
		nx1dArray<double>				m_truncangles_db;		// [1..18: num. sizes] phase function trunc. angles
		bool							m_isabsorber;

	private:
		void							ReleaseResources					( );
		bool							LookupIndicesAndWeights				( const nx1dArray<double>& h, double value, double* w1, size_t* idx1, double* w2, size_t* idx2 ) const;
		nxString						LoadDatabaseDirectoryFromRegistry	( );
		nxString						LoadStorageDirectoryFromRegistry	( );
		bool							WriteHeightWavelengthProfileToFile	( const char *filename,nx2dArray<double> extinction,nx1dArray<double> heightsm,std::vector<double> lambdas ) const;
		bool							WritePhaseFunctionFiles				( const std::vector<double>& lambdasnm,const nx1dArray<double>& reff,const nx1dArray<double>& cosscatangle,const nx3dArray<double>& phasefunct,nx1dArray<nxString>& filenames ) const;
		bool							WriteMasterFile						( const nx1dArray<double>& heightsm,const nx1dArray<nxString>& filenames,const nxString& masterfilename ) const;
		bool							ReadPhaseFunctionFiles				( const nx1dArray<double>& lambdasnm,const nx1dArray<double>& reff,const nx1dArray<double>& cosscatangle,const nx3dArray<double>& phasefunct,nx1dArray<nxString>& filenames );
		bool							ReadMasterFile						( const nx1dArray<double>& heightsm,const nx1dArray<nxString>& filenames,const nxString& masterfilename );
		bool							LoadHeightWavelengthProfileFromFile	( const char* filename );
		bool							ConfigureHeightInterpolation		( double heightm );
		bool							IsValidDatabaseEffectiveSize		( double reff_um)		{ return ((reff_um+1e-6)>m_effectivesize_db.At(0) && (reff_um-1e-6)<m_effectivesize_db.At(m_effectivesize_db.size()-1));	}
		bool							IsValidDatabaseWavelength			( double lambda_nm )	{ return ((lambda_nm+1e-6)>m_wavelennm_db.At(0) && (lambda_nm-1e-6)<m_wavelennm_db.At(m_wavelennm_db.size()-1)); }
		bool							CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs, double* forwardscatterfrac  ) const;
		bool							CalculatePhaseMatrix				( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix) const;

										skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB	( const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB& other ); // Dont allow copy constructor
		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB& operator =						( const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB& other ); // Dont allow assignment operator	

	public:
										skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB	();
		virtual						   ~skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB	()  override;
		bool							ConfigureRayTracingByExtinctionProfile							( nx2dArray<double> *cloudscatextinct,double minraytracedelta,std::vector<double> *raytracealts );
		bool							SetIceRadiusAndCrossSectionProfile								( nx2dArray<double> icesizeprofile,std::vector<double> lambdasnm,nx2dArray<double> *cloudscatextinct,const char *configdir );
		bool							SetIceRadiusAndNumDensityProfileGivenTau						( std::vector<double> lambdasnm,double hCloudTopKm,double cloudthickM,double opticaldelta,double particlesize,double lambdatau,double tau,nx2dArray<double> *iceprofile,nx2dArray<double> *cloudscatextinct,std::vector<double> *raytracealts,double minraytracedelta,double *iwp,const char *configdir );
		bool							SetupCloudPropertiesGaussGeneric								( nx2dArray<double> *icedata,double hCloudTopKm,double cloudthickm,double particlesize,double maxnumdensity,double opticaldelta );
		bool							SetLocalDirectory												( const char *configDir );	
		bool							SetDefaultTruncationAngles										();
		bool							LoadHeightWavelengthProfileFromMasterFile						( const char* masterfilename );
		bool							SetDBArrays														( );
		bool							GetPhaseFunctsFromDBAndSetMembers								( const nx1dArray<double>& heightsm,const nx1dArray<double>& reff,const std::vector<double>& lambdasnm,nxString &masterfilename );
		bool							GetXSectsFromDBAndSetMembers									( const nx1dArray<double>& heightsm,const nx1dArray<double>& reff,const std::vector<double>& lambdasnm,nx2dArray<double> *cloudcrosssection );
		bool							GetXSectsFromDBAndSetMembersGivenTau							( const nx1dArray<double>& heightsm,const nx1dArray<double>& reff,const std::vector<double>& lambdasnm,nx1dArray<double> *numden,nx2dArray<double> *cloudscatextinct,double lambdatau,double tau,std::vector<double> *raytracealts,double minraytracedelta,double *iwp );
		bool							AdjustNumAndExtinctionToFitTau									( nx1dArray<double> heightsM,nx2dArray<double> *extinction,nx1dArray<double> *numden,size_t idxtau,double tau,double opticaldelta );
		void							SetIsScatterer													( bool isscatterer)		{ m_isabsorber = !isscatterer; }
		//bool							SetExtinctionTable												( const nx2dArray<double>& extinction, const std::vector<double>& wavelens, const nx1dArray<double>& heights_meters);
		bool							SetCrossSectionTable											( const nx2dArray<double>& extinction, const std::vector<double>& wavelens, const nx1dArray<double>& heights_meters);
		double							TruncateAndComputeDeltaFraction									( nx1dArray<double> *phasefunct,double thetaC );
		double							ComputeAsymmetryParameter										( nx1dArray<double> *phasefunct,double thetaC );
//		bool							TruncatePhaseFunction											( nx1dArray<double> scatangle,nx3dArray<double> *phasefunct,const nx1dArray<double>& reff,const std::vector<double>& lambdasnm );

	public:
		virtual bool					SetAtmosphericState			( skClimatology* neutralatmosphere, const GEODETIC_INSTANT& pt, bool* crosssectionschanged  );			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool					InternalClimatology_UpdateCache	( const GEODETIC_INSTANT& /*pt*/)  override { return true;}
		virtual bool					CalculateCrossSections		( double wavenumber, double* absxs, double* extxs, double* scattxs)  override;														//!< Calculate cross-sections at the specified wave-number.
		virtual bool					CalculatePhaseMatrix		( double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix)  override;
		virtual bool					IsScatterer					() const  override	{ return !m_isabsorber;}												//!< Returns true if this particles scatters radiation
		virtual bool					IsAbsorber					() const  override	{ return m_isabsorber;}												//!< Returns true if this particles absorbs radiation radiation
};

