

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry		2009-6-16*/
/** \internal
 *	A class that hold one tabulated phase matrix entry for one wavelength at
 *	one altitude.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry
{
	private:
		double						m_roundedwavelen;				//!< The wavelength in nanometers multipled by 1000 and rounded to nearest integer (used as an exact index)
		SKTRAN_GridDefBase_V2		m_phasematrix;					//!< The scalar phase matrix  array (typically 181 or 361 elements)
		SKTRAN_GridDefBase_V2		m_cosscatterangle;				//!< The cosine of the scattering angle array, must be the same size as the phase matrix array

	private:
		void						ReleaseResources();
//		bool						DeepCopy			( const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry& other );

	public:
									skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry();
		bool						GetPhaseMatrix		( double cosangle, double* phasematrix ) const;
		bool						ConfigureEntry		( double wavelennm, const double* phasematrix, const double* cosangle, size_t numpoints);
		double						WavelengthIndex		() const { return m_roundedwavelen;}
		static double				RoundTheWavelength	( double nm ) { return  floor(1000*nm + 0.5);}
};

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry		2009-6-16*/
/** \internal
 *	a class that holds an array of phase matrices and cos(scatter angles) for a set of 
 *	wavelengths at any one altitude. The phase matrices for each wavelength at this height
 *	can be read in from a text file. Multiple instances of this class can be created inside
 *	class #skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength to create height profile
 *	of phase matrices.
 *
 *	\par File format
 *	I have expressed format using line numbers for clarity. In reality arrays can be spread over multiple 
 *	lines. White space and end-of-lines are skipped
 *	Line   1:		numwavelengths		numscatteringangles
 *	Line   2:		An array[numscatteringangles] of cos(scatterangle), can be distributed over several lines.
 *	Line   3:		Wavelen_nm[0]		An Array[numscatteringangle] of phasematrix for wavelength 0  
 *	Line   3:		Wavelen_nm[1]		An Array[numscatteringangle] of phasematrix for wavelength 1  
 *	Line   n+2:		Wavelen_nm[n-1]		An Array[numscatteringangle] of phasematrix for wavelength n-1  
 *
 *
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_TabulatedPhaseMatrix_HeightEntry
{
	private:
		double														m_heightmeters;					//!< The height assocaited with this set of wavelength entries
		size_t														m_numwavelens;					//!< The number of wavelength entries at this altitude
		size_t														m_maxwavelens;					//!< The maximum number of entries allowed
		skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry*	m_wavelenentries;				//!< the array of weavelength entries.


	private:
		void														ReleaseResources();
		bool														AllocateEntries	( size_t maxentries );

	public:
																	skOpticalProperties_TabulatedPhaseMatrix_HeightEntry	();
																   ~skOpticalProperties_TabulatedPhaseMatrix_HeightEntry	();
																	skOpticalProperties_TabulatedPhaseMatrix_HeightEntry	( double heightm ); // reserved for std::comparator functions (like lower_bound)
		bool														operator <										( const skOpticalProperties_TabulatedPhaseMatrix_HeightEntry & other ) const { return m_heightmeters < other.m_heightmeters;}
		double														Height											() const	{ return m_heightmeters;}
		bool														LoadWavelengthEntriesFromHeightFile				( double heightm, const char* filename );
//		bool														DeepCopy										( const skOpticalProperties_TabulatedPhaseMatrix_HeightEntry& other );
		bool														FindWavelengthEntry								( double wavelennm, const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry** entry ) const;
};

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength		2009-6-16*/
/** \ingroup skopticalpropmisc
 *	An Optical Properties class that allows a user to specify the scattering
 *	phase matrix as a function of altitude for a discrete set of wavelengths.
 *	The class interpolates the phase matrix in both scattering angle and altitude.
 *	It does not interpolate in wavelength. The cross section
 *	is provided by another Optical Properties object
 *
 *	\par Wavelength rounding
 *	The code rounds the wavelength to the
 *	nearest 1/1000 th of a nanometer and uses that "integer" value for finding
 *	matching wavelength entries.
 *
 *	\par Specifying the phase matrices
 *	The basic problem is that we are tring to specify the phase matrix over
 *	the scattering angle range (0 to 180 degrees) for a range of wavelengths
 *	and a range of altitudes. The first point to note is that we do NOT specify
 *	the phase matrix on a scattering angle grid  but on a a cosine of scattering angle
 *	grid. This must be in ascending (cosine) order so it starts at -1 (cos 180) and
 *	goes to +1 (cos 0).
 *
 *	There is a primary master file that specifies the altitudes at which phase matrices
 *	are defined. Each line of the file specifies the height at which the phase matrix is
 *	defined and the name of a text file containing the phase matrices and cos(scattering angle) for
 *	each wavlength. 
 *
 *	Line 1:		h_meters[0]			RelativePathToWavelengthFileFor_Height 0
 *	Line 2:		h_meters[1]			RelativePathToWavelengthFileFor_Height 1
 *	Line 3:		h_meters[2]			RelativePathToWavelengthFileFor_Height 2
 *	Line n-1:	h_meters[n-1]		RelativePathToWavelengthFileFor_Height n-1
**/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength : public skOpticalProperties
{
	private:
				  std::list<skOpticalProperties_TabulatedPhaseMatrix_HeightEntry>					m_heightentries;		// The cache of height entries (each height has phase matrices for an array of wavelengths)
		typedef	  std::list<skOpticalProperties_TabulatedPhaseMatrix_HeightEntry>::iterator			iterator;				// Iterates over the colelction
		typedef	  std::list<skOpticalProperties_TabulatedPhaseMatrix_HeightEntry>::const_iterator	const_iterator;

	private:
		skOpticalProperties*									m_crosssectionobject;			// The object used to retrieve cross-sections.
		skOpticalProperties_TabulatedPhaseMatrix_HeightEntry*			m_lowerheight;					// The lower height entry used to interpolate height entries
		double														m_f1;							// the interpolating factor for the lower height entry
		skOpticalProperties_TabulatedPhaseMatrix_HeightEntry*			m_upperheight;					// the uper height entry
		double														m_f2;							// the interpolating factor for the upepr height entry.

	private:
		bool								CheckCrossSectionValid		() const;
		bool								ConfigureHeightInterpolation( double heightm );
		void								ReleaseResources			();
		bool								CalculatePhaseMatrix		( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix ) const;

											skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength( const skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength& other );	// dont allow copy constructor
		skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength& operator =						 ( const skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength& other );	// Dotn allow assignment operator

	public:
											skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength();
		virtual							   ~skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength() override;
		bool								SetCrossSectionObject						( skOpticalProperties* crosssectionobject );
		bool								LoadHeightWavelengthProfileFromMasterFile	( const char* masterfilename );

	public:
		virtual bool						SetAtmosphericState			( skClimatology* neutralatmosphere)  override;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						SetLocation					( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )  override;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						InternalClimatology_UpdateCache	( const GEODETIC_INSTANT& pt )  override;
		virtual bool						CalculateCrossSections		( double wavenumber, double* absxs, double* extxs, double* scattxs) override;														//!< Calculate cross-sections at the specified wave-number.
		virtual bool						CalculatePhaseMatrix        ( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)  override;	//!< Calculate the phase matrix at the specified wavenumber and scattering angle
		virtual bool						IsScatterer					() const  override;
		virtual bool						IsAbsorber					() const  override;
};
