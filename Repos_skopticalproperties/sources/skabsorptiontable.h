
/*-----------------------------------------------------------------------------
 *					class skOpticalProperty_AdditionalStateInfo
Key		2013-11-19*/
/** A class used as part of caching convolved high resolution spectra. The problem to be
 *	addressed is that the convolving code will cache convolved high resolution
 *	cross-sections but the convolved spectra depend upon the same parameters
 *	such as pressure, temperature that the high cross-section
 *	depends upon.  Hence this small "key" class is used to address the problem.
 *
 *	The high resolution skOpticalProperties class provides a virtual method
 *	to set the key to the relevant state parameters. This set of numbers is
 *	then used to generate the key to cache the the convolved spectrum.
 *
 *	NOte that wavelength is 
 *

 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperty_AdditionalStateInfoKey
{
	private:
		double						m_stateparameters[4];				// A static array of state parameters. We could generalise and use std::vector but this is probably sufficient for many years.
		int							m_numstateparams;

	public:
									skOpticalProperty_AdditionalStateInfoKey();
								   ~skOpticalProperty_AdditionalStateInfoKey();
		bool						SetKeyStateParameters( const double params[], size_t numparams);
		void						Clear		();
		bool						operator == (const skOpticalProperty_AdditionalStateInfoKey& other ) const;
		bool						operator <  (const skOpticalProperty_AdditionalStateInfoKey& other ) const;
};


/*-----------------------------------------------------------------------------
 *					class skOpticalProperty_AdditionalStateInfo		2013-6-24*/
/** An interface exposed by some optical property classe that is used to
 *	associate a given atmospheric state with a given set of cross-sections
 *	via a unique index. This was first develeoped for convolution algorithms which
 *	need to cache a set of convolved cross-sections for different atmospheric conditions.
 *
 *	Many of the skOpticalProperty cross-sections are only dependent upon 
 *	temperature (or pressure and temperature) so it relatively easy to assign a unique index
 *	with each atmospheric state using justthe temperature.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperty_AdditionalStateInfo
{
	public:
		virtual bool		KeyedIndexFromAtmosphericState( skClimatology* neutralatmosphere, const GEODETIC_INSTANT& pt, skOpticalProperty_AdditionalStateInfoKey* index  ) = 0;
};


/*-----------------------------------------------------------------------------
 *					skOpticalProperty_AdditionalStateInfo_NotDependent		2013-6-24*/
/** This is a class/interface that can be used to represent skOpticalProperties that
 *	have no dependency upon atmospheric state.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperty_AdditionalStateInfo_NotDependent: public skOpticalProperty_AdditionalStateInfo
{
	public:
		virtual bool		KeyedIndexFromAtmosphericState( skClimatology* /*neutralatmosphere*/, const GEODETIC_INSTANT& /*pt*/, skOpticalProperty_AdditionalStateInfoKey* index ) {double d[] = {0}; index->SetKeyStateParameters( d, 0); return true;}
};

/*-----------------------------------------------------------------------------
 *					class skOpticalProperty_AdditionalStateInfo_TemperatureDependent		2013-6-24*/
/**	This is a class/interface that can be used to represent skOpticalProperties
 *	cross-sections that depend only upon atmospheric temperature.
**/
/*---------------------------------------------------------------------------*/

class skOpticalProperty_AdditionalStateInfo_TemperatureDependent: public skOpticalProperty_AdditionalStateInfo
{
	public:
		virtual bool		KeyedIndexFromAtmosphericState( skClimatology* neutralatmosphere, const GEODETIC_INSTANT& pt, skOpticalProperty_AdditionalStateInfoKey* index ) ;
};

/*-----------------------------------------------------------------------------
 *					class skOpticalProperty_AdditionalStateInfo_TemperatureDependent		2013-6-24*/
/**	This is a class/interface that can be used to represent skOpticalProperties
 *	cross-sections that depend only upon atmospheric temperature and pressure.
**/
/*---------------------------------------------------------------------------*/

class skOpticalProperty_AdditionalStateInfo_PressTemperatureDependent: public skOpticalProperty_AdditionalStateInfo
{
	public:
		virtual bool		KeyedIndexFromAtmosphericState( skClimatology* neutralatmosphere, const GEODETIC_INSTANT& pt, skOpticalProperty_AdditionalStateInfoKey* index ) ;
};
/*-----------------------------------------------------------------------------
 *					skWavelengthToPSF_Table		2012-7-24*/
/** \internal
 *	Pure virtual class used to return the point spread function of an instrument
 *	which was used to measure a cross-section. This is used by several classes
 *	to report the instrument point spread function to convolution code so the
 *	instrument width can be properly accounted for..**/
/*---------------------------------------------------------------------------*/

class skWavelengthToPSF_Table
{
	public:
		virtual			   ~skWavelengthToPSF_Table		(){}	

		/*-----------------------------------------------------------------------------
		 *					GetInstrumentPSF_FWHM		2013-6-24*/
		/** Returns the Full Width Half Max of the point spread function at wavelength "nm"
		 *	used to measure the associated cross-section data. This value is accounted for 
		 *	when convolving the lab-measured cross-section data so we have consistency.
		 **/
		/*---------------------------------------------------------------------------*/

		virtual double		GetInstrumentPSF_FWHM		( double nm ) const = 0;

		/*-----------------------------------------------------------------------------
		 *					GetInstrumentPointSpacing		2013-6-24*/
		/** Returns the intstrument interpoint spacing a the specified wavelength
		 */
		/*---------------------------------------------------------------------------*/

		virtual double		GetInstrumentPointSpacing	( double nm ) const = 0;
};


/*-----------------------------------------------------------------------------
 *					class skWavelengthToPSF_TableConstant		2012-7-24*/
/** \internal
 *	A class used to report an instrument point spread function that is just
 *	a constant value in nanometers.
 **/
/*---------------------------------------------------------------------------*/

class skWavelengthToPSF_TableConstant : public skWavelengthToPSF_Table
{
	private:
		double	m_fwhm;
		double	m_pointspacing;

	public:
							skWavelengthToPSF_TableConstant	()													{ m_fwhm = 0; m_pointspacing= 0.01;}
		bool				DeepCopy						( const skWavelengthToPSF_TableConstant& other )	{ m_fwhm = other.m_fwhm; m_pointspacing = other.m_pointspacing; return true;}
		void				SetInstrumentPSF_FWHM			( double fwhm)										{ m_fwhm = fwhm;}
		void				SetInstrumentPointSpacing		( double spacing_nm)								{ m_pointspacing = spacing_nm;}
		virtual double		GetInstrumentPSF_FWHM			( double /*nm*/  ) const							{ return m_fwhm;}
		virtual double		GetInstrumentPointSpacing		( double /*nm*/ ) const								{ return m_pointspacing;}

};


/*-----------------------------------------------------------------------------
 *							2012-7-24*/
/** \internal
 *	A class used to report an instrument point spread function that is a
 *	constant wavenumber. The FTS used by Voigt to measure O3 is an example.
**/
/*---------------------------------------------------------------------------*/

class skWavelengthToPSF_TableConstantWavenumber : public skWavelengthToPSF_Table
{
	private:
		double	m_deltak;
		double	m_pointspacingcm1;			// Point sapcing in wavenumber space.

	public:
							skWavelengthToPSF_TableConstantWavenumber	()			{ m_deltak = 0; m_pointspacingcm1 = 0;}
		bool				DeepCopy							( const skWavelengthToPSF_TableConstantWavenumber& other )	{ m_deltak = other.m_deltak; m_pointspacingcm1 = other.m_pointspacingcm1; return true;}
		void				SetInstrumentPSF_FWHMWavenumber		( double deltak)    { m_deltak = deltak;}
		void				SetInstrumentPointSpacingWavenumber	( double spacingcm1){ m_pointspacingcm1 = spacingcm1;}

		virtual double		GetInstrumentPSF_FWHM				( double nm  ) const; 
		virtual double		GetInstrumentPointSpacing			( double nm  ) const;
};
/*-----------------------------------------------------------------------------
 *					class skWavelengthToPSF_Table		2012-7-24*/
/** \internal 
 *	Pure virtual class used to return the point spread function of an instrument
 *	used to measure the tables. This is used in class skOpticalProperties_UserDefinedAbsorption
 *	to tell the user the PSF of the instrument used to measure the tables.
 *
 **/
/*---------------------------------------------------------------------------*/

class skWavelengthToPSF_TableArray : public skWavelengthToPSF_Table
{
	private:
	         std::map< double, double >						m_entries;
	typedef  std::map< double, double >::iterator			iterator;
	typedef  std::map< double, double >::const_iterator		const_iterator;
	double		m_pointspacing;


	public:
							skWavelengthToPSF_TableArray	() { m_pointspacing = 0;}
		bool				DeepCopy						( const skWavelengthToPSF_TableArray& other )	{ m_entries = other.m_entries; m_pointspacing = other.m_pointspacing; return true;}
		bool				AddPSFEntry						( double nm, double fwhm);
		void				SetInstrumentPointSpacing		( double spacing_nm)	{m_pointspacing = spacing_nm;}
		virtual double		GetInstrumentPSF_FWHM			( double nm ) const;
		virtual double		GetInstrumentPointSpacing		( double /*nm*/  ) const { return m_pointspacing;}

};

/*---------------------------------------------------------------------------
 *					class sk_AbsorptionTabulatedTableEntry		2003-11-28*/
/**	\internal
 *	Base class for molecules that are purely absorbing.  This aids radiative transfer
 *	calculations as we do not need to consider scattering from molecules that are purely
 *	absorbing.
 **/
 /*-------------------------------------------------------------------------*/

class sk_AbsorptionTabulatedTableEntry
{
	private:
		double				m_t;					// The temperature in Kelvins
		nx1dArray<double>	m_nm;					// The wavelength in nanometers
		nx1dArray<double>	m_xsection;				// The cross-section in cm2.
		double				m_minnm;				// The minimum wavelength in nanometers
		double				m_maxnm;				// The maximum wavelength  in nanometers

	private:
		bool				CheckWavelengths					( );
		void				ClearMinMaxRange					( );


	public:
							sk_AbsorptionTabulatedTableEntry	( );
							sk_AbsorptionTabulatedTableEntry	( double t);
		bool				Configure							( double t, double* nm, intptr_t nmstrides, double *xs, intptr_t xsstrides, size_t npts );
		bool				Configure							( double t, nx1dArray<double>& nm, nx1dArray<double>& xsection);
		bool				Configure							( double t, const std::vector<double>& nm, const std::vector<double>& xsection );
		bool				GetCrossSection						( double nanometer, double* crosssection ) const;
		double				GetTemperature						( ) const												{ return m_t;}
		double				MinimumWavelengthEntry				( ) const												{ return m_minnm;}
		double				MaximumWavelengthEntry				( ) const												{ return m_maxnm;}
		bool				operator <							( const sk_AbsorptionTabulatedTableEntry& other ) const { return m_t <  other.m_t;}
		bool				operator ==							( const sk_AbsorptionTabulatedTableEntry& other ) const { return m_t == other.m_t;}
};



/*-----------------------------------------------------------------------------
 *					struct _sk_xsectarray		2012-7-26*/
/** \internal
 *	A small structure used to specify cross-sections as a function of
 *	wavelength. A lot of the O3 and NO2 cross-sections use this structure to
 *	specify cross-section tables.
 **/
/*---------------------------------------------------------------------------*/

struct _sk_xsectarray
{
	double	nm;
	double	xsect;
};


/*---------------------------------------------------------------------------
 *					class skOpticalProperties_UserDefinedAbsorption						2003-11-28*/
/**	\ingroup skopticalpropmisc
 *	Base class for molecules or particles that are either purely absorbing or purely
 *	scattering and have the cross-sections expressed in tables as a function of
 *	wavelength and temperature.The code was upgraded in July 2012 to account for
 *	the DBM O3 cross-sections which have different wavelength coverage for
 *	different temperatures.
 *
 *	Cross-sections are provided as a table of wavelength versus cross-section for each
 *	temperature. The tables are normally stored in a derived class. Only one table per
 *  temperature should be used. Each wavelength/cross-section table is assumed to be contiguous from the first
 *	to last wavelength and is linearly interpolated. The code will return an error if the user tries
 *	to get a cross-section at a wavelength outside the range of wavelengths of any of
 *	the wavelength/cross-section tables.
 *
 *	As an example the Daumont Brion Malicet (DBM) O3 tables consist of 5 temperature tables (218K,
 *	228K, 243K, 273K and 295K). The 218K table extends from 194.50 nm to 650 nm. The 3 tables 228K, 243K
 *	and 273K tables have a shorter range from 194.5nm to 520 nm. The last table, 295K, extends from
 *	195 nm to 830 nm. When the user supplies a temperature and wavelength the code determines
 *	which of the 5 temperature tables span the selected wavelength and from this selection chooses the temperature
 *	tables just hotter and just colder than the desired temperature. Cross-sections are then
 *	linearly interpolated in wavelength and temperature.
 *
 *	Wavelength interpolation will only interpolate: it will not extrapolate beyond the table. Thus if none of
 *	the temperature tables have any entries that span the desired wavelength then the code will return a fail status.
 *	On the other hand if only one of the temperature tables has entries that span the selected wavelength then it will be used
 *	regardless of the incoming temperature. For example asking for cross-sections at 750 nm from the DBM O3 at any temperature will
 *	always return the same cross-section as only the 295K table spans this wavelength.
 *
 *	Temperature interpolation is linear when multiple temperature tables are available for a given wavelength
 *	and truncates out-of-range temperatures to the nearest valid table.
 *	This is reasonable behaviour for most atmospheric applications as temperature variations in cross-sections are
 *	small and the number of cross-sections measurements with temperature are limited and may not
 *	fully span the range of temepratures encountered in the atmosphere.
 *
 *	Developers should ensure that temperature tables (typically in derived classes)
 *	are added in ascending temperature order. I dont think it actually needs this but I have never tested
 *	anything but this. The cross-sections and corresponding wavelengths (in nanometers)
 *	for each temperature setting must be specified in ascending order (the class will not accept wavelength data that are not in ascending order).
 *
 **/
/*-------------------------------------------------------------------------*/

class skOpticalProperties_UserDefinedAbsorption: public skOpticalProperties
{
	private:
			skClimatology*													m_backgroundatmosphere;
			bool															m_quietwavelengthtruncation;	// If true then dont fail a cross-section request outside the range of wavelengths
			RefractiveIndexDryAirSTP										m_refractiveindexair;		// used to convert wavelength in air to wavelength in vacuum.
			bool															m_wavelengthinvacuum;		// True if the wavelengthtables are in vacuum
			bool															m_isabsorber;
			double															m_temperature;
		    std::list< sk_AbsorptionTabulatedTableEntry >					m_temperature_entries;
	typedef std::list< sk_AbsorptionTabulatedTableEntry >::iterator			iterator;
	typedef std::list< sk_AbsorptionTabulatedTableEntry >::const_iterator	const_iterator;

	private:
		bool										CopyEntries									( const skOpticalProperties_UserDefinedAbsorption& other );
		sk_AbsorptionTabulatedTableEntry*			FetchNewOrExistingEntryAtTemperature		( double t );
		bool										InterpolateCrossSectionInTemperature		( double t, double nm, double* crosssection) const;

	protected:
		bool										AddEntry									( double t, double* nm,     int nmstride, double *xs, int xsstride, int npts );
		bool										AddAscendingWavenumberEntry					( double t, double* usercm, int cmstride, double *xs, int xsstride, int npts );
		bool										DeepCopyWithoutTables						( const skOpticalProperties_UserDefinedAbsorption& other );
		bool										DeepCopy									( const skOpticalProperties_UserDefinedAbsorption& other );

	private:
													skOpticalProperties_UserDefinedAbsorption	( const skOpticalProperties_UserDefinedAbsorption& other );	// Dont allow copy constructor
		skOpticalProperties_UserDefinedAbsorption&	operator =									( const skOpticalProperties_UserDefinedAbsorption& other );	// Dont allow assignment operator

	public:
													skOpticalProperties_UserDefinedAbsorption	()	{ m_quietwavelengthtruncation = false; m_wavelengthinvacuum = false; m_isabsorber = true; m_temperature = -99999;m_backgroundatmosphere = nullptr;}
		virtual									   ~skOpticalProperties_UserDefinedAbsorption	() override { if (m_backgroundatmosphere != nullptr) m_backgroundatmosphere->Release();};
		bool										Set_Temperature								( double kelvin);															// Set the temperature (default is do nothing)
		bool										AddUserEntry								( double kelvin, nx1dArray<double>& nm, nx1dArray<double>& xsection);
		bool										AddUserEntry								( double kelvin, const std::vector<double>& usernm, const std::vector<double>& userxsection);
		double										Temperature									() const				{ return m_temperature;}
		void										SetIsScatterer								( bool isscatterer)		{ m_isabsorber = !isscatterer;}
		void										ClearEntries								()						{ m_temperature_entries.clear();}
		void										SetVacuumWavelengths						( bool isvacuum )		{ m_wavelengthinvacuum = isvacuum;}		// True if the wavelengthtables are in vacuum
		double										AirToVacuumCorrection						( double  nm_airSTP ) const;
		void										SetQuietWavelengthTruncation				( bool bequiet)			{ m_quietwavelengthtruncation = bequiet;}

	public:
		virtual bool								SetAtmosphericState							( skClimatology* neutralatmosphere )  override;
		virtual bool								SetLocation									( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )  override;
		virtual bool								CalculateCrossSections						( double wavenumber, double* absxs, double* extxs, double* scattxs) override;
		virtual bool								IsScatterer									() const  override						{ return !m_isabsorber;}
		virtual bool								IsAbsorber									() const  override						{ return m_isabsorber;}
		virtual bool								InternalClimatology_UpdateCache				( const GEODETIC_INSTANT& /*pt*/)  override { return true;}
};

/*---------------------------------------------------------------------------
 *					class skOpticalProperties_TabulatedExtinction_HeightWavelength						2003-11-28*/
/**	\ingroup skopticalpropmisc
 *	A class for tabulating extinction (or absorption) as a function of altitude
 *	and wavelength. This was originally built for manually defining the extinction
 *	of henyey-Greenstein phase functions.
 **/
/*-------------------------------------------------------------------------*/

class skOpticalProperties_TabulatedExtinction_HeightWavelength: public skOpticalProperties
{
	private:
		bool								m_isabsorber;
		nx1dArray<double>					m_heights;
		nx1dArray<double>					m_wavelennm;
		nx2dArray<double>					m_extinction;
		size_t								m_idxh1;
		size_t								m_idxh2;
		double								m_h1;
		double								m_h2;

	private:
		void								ReleaseResources			();
		bool								LookupIndicesAndWeights		( const nx1dArray<double>& h, double value, double* w1, size_t* idx1, double* w2, size_t* idx2 ) const;
		bool								CalculateCrossSections		( double wavenumber, double* absxs, double* extxs, double* scattxs) const;
											skOpticalProperties_TabulatedExtinction_HeightWavelength( const skOpticalProperties_TabulatedExtinction_HeightWavelength& other );	// Dont allow copy constructor
		skOpticalProperties_TabulatedExtinction_HeightWavelength& operator =						( const skOpticalProperties_TabulatedExtinction_HeightWavelength& other );	// Dont allow assignment operator

	public:
											skOpticalProperties_TabulatedExtinction_HeightWavelength	();
		virtual							   ~skOpticalProperties_TabulatedExtinction_HeightWavelength	() override {};
		void								SetIsScatterer						( bool isscatterer)		{ m_isabsorber = !isscatterer; }
		bool								SetExtinctionTable					( const nx2dArray<double>& extinction, const nx1dArray<double>& wavelens, const nx1dArray<double>& heights_meters);
		bool								LoadHeightWavelengthProfileFromFile	( const char* filename );


	public:
		virtual bool						SetAtmosphericState			( skClimatology* neutralatmosphere)  override;									//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						SetLocation					( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )  override;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						CalculateCrossSections		( double wavenumber, double* absxs, double* extxs, double* scattxs)  override;														//!< Calculate cross-sections at the specified wave-number.
		virtual bool						IsScatterer					() const  override	{ return !m_isabsorber;}												//!< Returns true if this particles scatters radiation
		virtual bool						IsAbsorber					() const   override	{ return m_isabsorber;}												//!< Returns true if this particles absorbs radiation radiation
		virtual bool						InternalClimatology_UpdateCache	( const GEODETIC_INSTANT& /*pt*/)  override { return true;}
};

/*---------------------------------------------------------------------------
 *					class skOpticalProperties_O3_GomeBurrows				  2003-11-28*/
/** \ingroup o3skopticalprop
 *	Table of O3 cross-sections measured with the GOME flight instrument at medium wavelength resolution.
 *
 *	\par Spectral resolution
 *		- 231-307 nm, 0.20 nm (FWHM ?)
 *		- 307-316 nm, 0.20 nm
 *		- 311-405 nm, 0.17 nm
 *		- 405-611 nm, 0.29 nm
 *		- 595-794 nm, 0.33 nm
 *
 *	\par Temperature Range
 *	The cross-sections were measured at 5 temperatures covering normal stratospheric and tropospheric ranges. The paper does not
 *	provide any advice on how to interpolate in temperature. We use the standard (linear interpolation) technique provided by the base class.
 *	The spectra were measured at the following four temperatures.
 *		-# 202 K,
 *		-# 221 K,
 *		-# 241 K,
 *		-# 273 K
 *		-# 293 K
 *
 *	\par References
 *  J. P. Burrows, A. Richter, A. Dehn, B. Deters, S. Himmelmann, S. Voigt, and J. Orphal:
 *  "Atmospheric Remote-Sensing Reference Data from GOME: 2. Temperature-Dependent Absorption Cross Sections of O3 in the 231-794 nm Range",
 *	Journal of Quantitative Spectroscopy and Radiative Transfer 61, 509-517, 1999.
 *
 *
 *	Air measurements.
 *
 */
/*-------------------------------------------------------------------------*/
class skOpticalProperties_O3_GomeBurrows: public skOpticalProperties_UserDefinedAbsorption,
	                                      public skWavelengthToPSF_TableArray,								// Describes the instrument used to measure these cross-sections 
										  public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section

{
	private:
												skOpticalProperties_O3_GomeBurrows  ( const skOpticalProperties_O3_GomeBurrows& other );	// Dont allow copy constructor
		skOpticalProperties_O3_GomeBurrows&		operator =							( const skOpticalProperties_O3_GomeBurrows& other );	// Dont allow assignment operator

	public:
												skOpticalProperties_O3_GomeBurrows();
		virtual								   ~skOpticalProperties_O3_GomeBurrows() override {}
};

/*---------------------------------------------------------------------------
 *					class skOpticalProperties_O3_SciaBogumilV3				  2003-11-28*/
/** \ingroup o3skopticalprop
 *	The version 3.0 cross-sections of O3 measured by Bogumil et al. with the Sciamachyy flight
 *	instrument before launch. The cross-sections are measured from 230 nm to 1070 nm at five
 *	seperate temperatures,
 *		-# 203 K,
 *		-# 223 K,
 *		-# 243 K,
 *		-# 273 K
 *		-# 293 K
 *
 *	All measuremnts are at the instrinsic FWHM spectral resolution of sciamachy,
 *		- 0.32 nm below 311.7 nm
 *		- 0.21 nm between 310.8 nm and 402 nm
 *		- 0.52 nm between 402 nm and 597.8 nm
 *		- 0.47 nm between 597.8 nm and 781.6 nm
 *		- 0.62 nm between 781.6 nm and 1052.5 nm
 *		- 1.45 nm above 1052.5 nm
 *
 *	\par Experimental:
 *		- Light source: Xenon lamp, QTH lamp for 489.5 - 718.9 nm wavelength region
 *		- Cell:         Multiple reflection quartz cell (White optics)
 *		- Optical pathlength 505 cm and 985 cm
 *		- Gas mixture:  O3, O2, N2
 *		- Condition:    Flow measurement
 *		- Total pressure: 50 - 900 mbar
 *		- Spectrometer:   SCIAMACHY PFM Satellite Spectrometer (Eight channel grating spectrometer)
 *
 *	Note: Spectral Regions are cutted and concatenated at 287.7 nm, 302.2 nm, 310.8 nm, 316.3 nm,
 *  334.2 nm, 345.2 nm, 402 nm, 489.5 nm, 597.8 nm, 718.9 nm, 781.6 nm and 1052.5 nm.
 *
 *	\par Version history
 *		- Version 1.0 first published version
 *		- Version 2.0 correction for the memory effect in channel 1 (230 - 310.8 nm) applied
 *		- Version 3.0 baseline correction for channel 1 revised
 *
 *	\par References:
 *	K. Bogumil, J. Orphal, J.P. Burrows, J.M. Flaud: Vibrational progressions in the
 *	visible and near ultra-violet asorption spectrum of ozone. Chem Phys Lett, 349, pages 241-248, 2001
 */
/* ------------------------------------------------------------------------- */

class skOpticalProperties_O3_SciaBogumilV3: public skOpticalProperties_UserDefinedAbsorption,
	                                        public skWavelengthToPSF_TableArray,								// Describes the instrument used to measure these cross-sections 
											public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section


{
	private:
												skOpticalProperties_O3_SciaBogumilV3  ( const skOpticalProperties_O3_SciaBogumilV3& other );	// Dont allow copy constructor
		skOpticalProperties_O3_SciaBogumilV3&	operator =							( const skOpticalProperties_O3_SciaBogumilV3& other );	// Dont allow assignment operator

	public:
												skOpticalProperties_O3_SciaBogumilV3();
		virtual								   ~skOpticalProperties_O3_SciaBogumilV3() override {}
};


/*-----------------------------------------------------------------------------
 *					class  skOpticalProperties_O3_SciaBogumilV4		2012-7-27*/
/** \ingroup o3skopticalprop
 *	Version 4 is a reanalysis of the actual Bogumil data ( version 3). A manuscript was being
 *	prepared for publication (as of July 2012) on the revised cross-sections (private communication
 *	with Mark Weber). The Sciamachy people had some problems in the Huggins band with the
 *	Bogumil data as total ozone retrievals were biased 5%.
 */
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O3_SciaBogumilV4: public skOpticalProperties_UserDefinedAbsorption,
	                                        public skWavelengthToPSF_TableArray,								// Describes the instrument used to measure these cross-sections 
											public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section

{
	private:
												skOpticalProperties_O3_SciaBogumilV4	( const skOpticalProperties_O3_SciaBogumilV4& other );	// Dont allow copy constructor
		skOpticalProperties_O3_SciaBogumilV4&	operator =								( const skOpticalProperties_O3_SciaBogumilV4& other );	// Dont allow assignment operator


	public:
							skOpticalProperties_O3_SciaBogumilV4();
		virtual			   ~skOpticalProperties_O3_SciaBogumilV4() override {}
};

/*---------------------------------------------------------------------------
 *					class skOpticalProperties_O3_BassPaurQuadratic				  2003-11-28*/
/** \ingroup o3skopticalprop
 *	Calculates the O3 cross-section using the Bass-Paur quadratic temperature coefficients.
 *	This is in contrast to matching class #skOpticalProperties_O3_BassPaur which uses tables
 *	of temperature measurements. The data are a world standard but show signs of age and are
 *	of limited wavelength extent (no Chappuis measurements). They are are at high wavelength
 *	sampling and high wavelength resolution.
 *
 *	\par Wavelength Range
 *	The table covers the wavelength range 245.0180 .. 342.7800 nm in 1956 samples,  approx 0.05 nm per step.
 *
 *	\par Temperature Range
 *	The Bass-Paur data set uses a quadratic temperature function 
 *	where the quadratic coeffs are given for each wavelength.
 *	The valid temperature range is 203 K to 298 K. The code truncates temperature dependence
 *	curves to these limits.
 *
 *	\par Spectral Resolution
 *	The spectral resolution is better than 0.025 nm. 
 *
 *	\par Data Source
 *	These data are an exact replication of the data in file bp.par on the IGACO site, http://igaco-o3.fmi.fi/ACSO/files/cross_sections
 *	and was copied between July 15 and July 25 2012.  Note that the entries for 282.470 nm and 282.460 nm were swapped around as they are
 *	in the wrong order in the original files.
 *
 *	\par References
 *  -#	Bass, A. M. and Paur, R. J., The ultraviolet cross-sections of ozone, I. Measurements
 *		Proc. Quadrennial Ozone Symp. Halkidiki, Greece, Reidel, Dordrecht, pp. 606–610, 1984
 *  -#	Bass, A. M. and Paur, R. J., The ultraviolet cross-sections of ozone, Ii. Results and
 *		temperature dependence. Proc. Quadrennial Ozone Symp. Halkidiki, Greece, Reidel,
 *		Dordrecht, pp. 611–616, 1984
 *
 */
/*-------------------------------------------------------------------------*/

class skOpticalProperties_O3_BassPaurQuadratic: public skOpticalProperties,
	                                            public skWavelengthToPSF_TableConstant,								// Describes the instrument used to measure this cross-section
												public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section
{

	public:
		struct sk_dummy_basspaur03xsectentry
		{
			double  nm;
			double	c0;
			double  c1;
			double  c2;
		};

	private:
		skClimatology*													m_backgroundatmosphere;
		bool															m_isdirty;
		double															m_temperature;
		static const sk_dummy_basspaur03xsectentry						m_o3coeffs[1956];


	private:
		bool						CalculateCrossSections					( double wavenum, double* absxs, double* extxs, double* scattxs ) const;
		bool						DeepCopy								( const skOpticalProperties_O3_BassPaurQuadratic& other );
									skOpticalProperties_O3_BassPaurQuadratic( const skOpticalProperties_O3_BassPaurQuadratic& other );	// Dont allow copy constructor
		skOpticalProperties_O3_BassPaurQuadratic&  operator =				( const skOpticalProperties_O3_BassPaurQuadratic& other );	// Dont allow assignment operator

	public:
									skOpticalProperties_O3_BassPaurQuadratic	();
		virtual					   ~skOpticalProperties_O3_BassPaurQuadratic	() override;
		double						Temperature					() const				{ return m_temperature;}
		bool						Set_Temperature				( double kelvin);															// Set the temperature (default is do nothing)
		double						BassPaurCrossSection		( const sk_dummy_basspaur03xsectentry* entry ) const;

	public:
		virtual bool				SetAtmosphericState			( skClimatology* neutralatmosphere)  override;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool				SetLocation					( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )  override;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool				InternalClimatology_UpdateCache	(const GEODETIC_INSTANT& /*pt*/ )  override { return true;}
		virtual bool				CalculateCrossSections		( double wavenumber, double* absxs, double* extxs, double* scattxs )  override;														//!< Calculate cross-sections at the specified wave-number.
		virtual bool				IsScatterer					() const	 override { return false;}												//!< Returns true if this particles scatters radiation
		virtual bool				IsAbsorber					() const	 override { return true;}												//!< Returns true if this particles absorbs radiation radiation
};


/*-----------------------------------------------------------------------------
 *					class skOpticalProperties_O3_BassPaur		2012-7-25*/
/** \ingroup o3skopticalprop
 *	Tables of the O3 cross-section at 6 different temperatures made by Bass-Paur.
 *	A complementary class #skOpticalProperties_O3_BassPaurQuadratic provides the same data
 *	with quadratic temperature interpolation. The data are a world standard but show signs of age and are
 *	of limited wavelength extent (no Chappuis measurements). They are are at high wavelength
 *	sampling and high wavelength resolution.
 *
 *	\par Wavelength Range
 *	The table covers the wavelength range 245.0180 .. 342.7800 nm in 1956 samples,  approx 0.05 nm per step.
 *
 *	\par Temperature Range
 *	The Bass Paur tables provide measurements at 6 temperatures: 203 K, 223 K, 246 K, 273 K, 276 K and 280 K.
 *
 *	\par Spectral Resolution
 *	The spectral resolution is better than 0.025 nm. 
 *
 *	\par Data Source
 *	These data are an exact replication of the data files bp_203clc.dat, bp_223clc.dat, bp_246clc.dat,
 *	bp_273clc.dat, bp_276clc.dat and bp_280clc.dat on the IGACO site, http://igaco-o3.fmi.fi/ACSO/files/cross_sections
 *	and were copied on July 25 2012. Note that the entries for 282.470 nm and 282.460 nm were swapped around as they are
 *	in the wrong order in teh original files.
 *
 *	\par References
 *  -#	Bass, A. M. and Paur, R. J., The ultraviolet cross-sections of ozone, I. Measurements
 *		Proc. Quadrennial Ozone Symp. Halkidiki, Greece, Reidel, Dordrecht, pp. 606–610, 1984
 *  -#	Bass, A. M. and Paur, R. J., The ultraviolet cross-sections of ozone, Ii. Results and
 *		temperature dependence. Proc. Quadrennial Ozone Symp. Halkidiki, Greece, Reidel,
 *		Dordrecht, pp. 611–616, 1984
 */
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O3_BassPaur: public skOpticalProperties_UserDefinedAbsorption,
	                                   public skWavelengthToPSF_TableConstant,								// Describes the instrument used to measure this cross-section
									   public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section
{
	private:
												skOpticalProperties_O3_BassPaur	( const skOpticalProperties_O3_BassPaur& other );	// dont allow copy constructor
		skOpticalProperties_O3_BassPaur&		operator =						(const skOpticalProperties_O3_BassPaur& other );	// Dont allow assignment operator

	public:
												skOpticalProperties_O3_BassPaur();
		virtual								   ~skOpticalProperties_O3_BassPaur(){}

};

/*-----------------------------------------------------------------------------
 *					class skOpticalProperties_O3_DaumontBrionMalicet		2009-11-6*/
/** \ingroup o3skopticalprop
 *	Tabulated high resolution cross-sections of O3 measured by Daumont, Brion and Malicet in the early
 *	1990's. Voigt et al. 2001. The wavelength range is a little variable with temperature but covers the
 *	entire UV to NIR region, from 194.50 nm to 830.00 nm at 0.01 to 0.02 nm resolution.
 *	The cross-section data were collected at 0.01-0.02 nm resolution and each
 *	wavelength/cross-section table varies in size from 22,052 to 63,501 entries. The data consists of
 *	5 tables of wavelength versus cross-section for 5 temperatures.
 *
 *	The class derives from skOpticalProperties_UserDefinedAbsorption which was upgraded in July 2012 to properly
 *	interpolate the DBM wavelength/cross-section tables.
 *
 * \par Temperature Range
 *	measurements are provided at 5 temperatures covering typical stratospheric and tropospheric conditions:
 *		-# 218 K
 *		-# 228 K
 *		-# 243 K
 *		-# 273 K
 *		-# 295 K
 *	
 *	\par Wavelength Range
 *	The wavelength range of each temperature table is slightly different and is given below. Note that most of
 *	the temperature variation occurs in the huggins band between 315 and 360 nm
 *
 *		- 218K -> 194.50nm to 650.01nm
 *		- 228K -> 194.50nm to 520.01nm
 *		- 243K -> 194.50nm to 519.01nm
 *		- 273K -> 299.50nm to 520.01nm
 *		- 295K -> 195.00nm to 830.00nm
 *
 *	We looked into temperature interpolation and while DBM suggest that a quadratic interpolation scheme
 *	(page 269 of paper 3.) they do not indicate an explicit technique. We tested several quadratic fitting routines and
 *	found that a truncated linear fit in temperatuer was visually more appealing than any of the quadratic fits and had none
 *	of the undesirable artifacts (excessive curvature etc) that naturally arises with quadratic curve fitting.
 *
 *
 *	\par Data Source
 *	These data are an exact replication of the data files O3_CRS_BDM_218K.dat, O3_CRS_BDM_228K.dat, O3_CRS_BDM_243K.dat,
 *	O3_CRS_BDM_273K.dat and O3_CRS_BDM_295K.dat on the IGACO site, http://igaco-o3.fmi.fi/ACSO/files/cross_sections.
 *	The files were copied sometime between July 16-July 25 2012.
 *
 *	\par References
 *	Daumont Brion Malicet O3 cross-sections as published in,
 *		-# D. Daumont et al. Ozone UV spectroscopy I: Absorption cross-sections at room temperature
 *		   Journal of Atmospheric Chemistry Vol. 15, 1992.
 *		-# J. Brion et al. High-resolution laboratory absorption cross section of O3. Temperature effect
 *		   Chemical Physics Letters Vol. 213, No. 5,6, 1993.
 *		-# J. Malicet et al. Ozone UV spectroscopy II: Absorption cross-sections and temperature dependence
 *		   Journal of Atmospheric Chemistry Vol. 21, 1995.
 * 		-# J. Brion et al. Absorption spectra measurements for the ozone molecule in the 350-830 nm region
 *		   Journal of Atmospheric Chemistry Vol. 31, 1998.
*/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O3_DaumontBrionMalicet : public skOpticalProperties_UserDefinedAbsorption,
	                                               public skWavelengthToPSF_TableConstant,								// Describes the instrument used to measure this cross-section
												   public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section

{
	private:
												skOpticalProperties_O3_DaumontBrionMalicet	( const skOpticalProperties_O3_DaumontBrionMalicet& other );	// dont allow copy constructor
		skOpticalProperties_O3_DaumontBrionMalicet&		operator =							( const skOpticalProperties_O3_DaumontBrionMalicet& other );	// Dont allow assignment operator



	public:
												skOpticalProperties_O3_DaumontBrionMalicet();
		virtual								   ~skOpticalProperties_O3_DaumontBrionMalicet() override {}

};

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_FTSVoigt		2009-11-6*/
/** \ingroup o3skopticalprop
 *	Tabulated high resolution cross-sections of O3 measured by Voigt et al. 2001. The
 *	resolution is a constant 5 cm-1 (0.027 nm at 230 nm to 0.36 nm at 850 nm)
 *
 *	\par Temperatures
 *	The Voigt paper presents a new exponential interpolation algorithm for
 *	interpolating cross-sections in temperature. This class does not yet use this technique but still
 *	uses truncated linear interpolation. Measurements are provided at the following five temperatures,
 *		-# 203 K
 *		-# 223 K
 *		-# 246 K
 *		-# 280 K
 *		-# 293 K
 *
 *	\par References
 *		-# S. Voigt, J. Orphal, K. Bogumil, J.P. Burrows, The temperature dependence (203–293 K)
 *			of the absorption cross sections of O3 in the 230–850 nm region measured by
 *			Fourier-transform spectroscopy. Journal of Photochemistry and Photobiology A: Chemistry Vol. 143, 2001.
 *		-# J. Orphal et al. A critical review of the absorption cross-sections of O3 and NO2 in the 240-790 nm region
 *			ESA Technical Note MO-TN-ESA-GO-0302, 2002.
*
 *	\par Abstract of Paper 1.
 *	Absolute absorption cross sections of O3 were measured in the 230–850 nm (11765–43478 cm?1) region
 *	at five different temperatures (203–293 K) using a Fourier-transform spectrometer, at a spectral
 *	resolution of 5.0 cm?1 (corresponding to about 0.027 nm at 230 nm and to about 0.36 nm at 850 nm).
 *	The spectral accuracy of the data is better than 0.1 cm?1 — about 0.5 pm at 230 nm and about 7.2 pm
 *	at 850 nm — validated by recording of I2 absorption spectra in the visible using the
 *	same experimental set-up. O3 absorption spectra at different concentrations were recorded at five
 *	different sample temperatures in the range 203–293 K, and at each temperature at two total
 *	pressures (100 and 1000 mbar) using O2/N2 mixtures as buffer gas. Within the limits of
 *	experimental uncertainties, no influence of total pressure on the O3 spectrum was observed in
 *	the entire spectral region, as expected from the short lifetimes of the upper electronic states of
 *	O3. The temperature dependence of the O3 absorption cross sections is particularly strong
 *	in the Huggins bands between 310 and 380 nm, as observed in previous studies. An empirical
 *	formula is used to model the temperature dependence of the O3 absorption cross sections
 *	between 236 and 362 nm, a spectral region that is particularly important for atmospheric remote-sensing
 *	and for photochemical modelling.
 *
 *	\par Header details from distributed Data Files
 *	ESA Study 11340/95/NL/CN UV - Cross-Sections in the UV and Visible
 *  O3 ABSORPTION CROSS-SECTIONS AT 203-293K
 *
 *	Low pressure = 100 mbar total pressure*
 *	High pressure = 1000 mbar total pressure
 *
 *	S. Voigt, J. Orphal, and J. P. Burrows
 *	University of Bremen - Institute of Environmental Physics
 *	P. O. Box 33 04 40
 *	D-28334 Bremen, Germany
 *	Tel. + 49 (0)421 218 3526
 *	e-mail: Susanne.Voigt@iup.physik.uni-bremen.de
 *
 *	Wavelength range:	230 - 850 nm
 *	Wavenumber range:	11752 - 43315 cm-1
 *	Spectral Resolution:	5 cm-1
 *
 *	EXPERIMENTAL
 *	Cell:	Multiple reflection quartz cell (White optics)
 *	optical pathlength 505 cm/120 cm (single path without White optics)
 *
 *	Spectrometer:	BRUKER IFS 120HR Fourier-Transform-Spectrometer
 *
 *	Wavenumber Range	Light Source		Detector
 *	----------------	------------		--------
 *	29500-43500 cm-1	Xe lamp			UV diode
 *	20000-33000 cm-1	Xe lamp			GaP diode
 *	12000-25000 cm-1	QTH lamp		Si diode
 *
 *	Note: Spectral Regions are cutted and concatenated at 22000 cm-1,
 *	31000 cm-1 and 33000 cm-1.
**/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O3_FTSVoigt : public skOpticalProperties_UserDefinedAbsorption,
	                                    public skWavelengthToPSF_TableConstantWavenumber,					// Describes the instrument used to measure this cross-section
										public skOpticalProperty_AdditionalStateInfo_TemperatureDependent	// Describes what state parameters affetct the cross-section
{
	private:
		bool								m_usinglowpressure;

	private:
		bool								UseLow100mBPressureEntries		();
		bool								UseHigh1000mBPressureEntries	();
		bool								DeepCopy						( const skOpticalProperties_O3_FTSVoigt& other );
											skOpticalProperties_O3_FTSVoigt( const skOpticalProperties_O3_FTSVoigt& other );	// Dont allow copy constructor
		skOpticalProperties_O3_FTSVoigt&	operator =				( const skOpticalProperties_O3_FTSVoigt& other );	// Dont allow assignment operator

	public:
											skOpticalProperties_O3_FTSVoigt();
		virtual 						   ~skOpticalProperties_O3_FTSVoigt() override {}

};


/*-----------------------------------------------------------------------------
 *					skWavelengthToPSF_SerdyuchenkoV1		2012-7-27*/
/** \internal
 *	The class that hold the O3 Serdyuchenko point spread function as a function of
 *	wavelength
 **/
/*---------------------------------------------------------------------------*/

class skWavelengthToPSF_SerdyuchenkoV1 : public skWavelengthToPSF_Table
{
	private:
		skWavelengthToPSF_TableConstantWavenumber					m_psfwavenum;

	public:
							skWavelengthToPSF_SerdyuchenkoV1		();	
		virtual			   ~skWavelengthToPSF_SerdyuchenkoV1		(){}
		bool				DeepCopy								( const skWavelengthToPSF_SerdyuchenkoV1& other )	{ return m_psfwavenum.DeepCopy( other.m_psfwavenum);}
		virtual double		GetInstrumentPSF_FWHM					( double nm ) const;
		virtual double		GetInstrumentPointSpacing				( double /*nm*/ ) const { return 0.01;}
};

/*-----------------------------------------------------------------------------
 *					class skOpticalProperties_O3_SerdyuchenkoV1		2012-7-26*/
/** \ingroup o3skopticalprop
 *	Tabulated high resolution cross-sections of O3 measured by Serdyuchenko et al. 2011. These are
 *	a new generation of cross-sections with measurements from 193 K to 293 K in 10 K steps at high resolution.
 *
 *	\par Temperatures
 *	The Serdychenko measurements try to avoid the limitations of interpolation by making twice as many
 *	meautrements. Measurements are provided at the following eleven temperatures,
 *		-# 193 K
 *		-# 203 K
 *		-# 213 K
 *		-# 223 K
 *		-# 233 K
 *		-# 243 K
 *		-# 253 K
 *		-# 263 K
 *		-# 273 K
 *		-# 283 K
 *		-# 293 K
 *
 *	\par References
 *		-# A. Serdyuchenko, V. Gorshelev, M. Weber, J.P. Burrows: New broadband high resolution ozone absorption cross-sections
 *		   in Spectroscopy Europe, http://www.spectroscopyeurope.com/articles/55-articles/3082-new-broadband-high-resolution-ozone-absorption-cross-sections
 *		-# Peer reviewed paper submitted in summer 2012.

 *	\par Abstract 
 *	In this article, we report on the research to improve our knowledge of the ozone absorption cross-sections.
 *	This is required for active and passive remote sensing applications yielding the total column and profiles of ozone.
 *	New laboratory measurements provide data for a wide spectral range in the ultraviolet (UV), visible (vis) and
 *	near infrared (NIR) regions at a spectral resolution of 0.02 nm. An absolute accuracy of about 3% or better
 *	and wavelength accuracy better than 0.005 nm throughout the spectral range have been achieved at 11 
 *	temperatures from 195 K to 293 K.
 *
 *	Comparison of the available ozone cross-sections with our new dataset shows good agreement within the
 *	uncertainty limits. This new cross-section dataset improves the ozone data quality which is required for stratospheric
 *	ozone trend studies and the determination of tropospheric ozone abundance.
 *
 *	Ozone is the most important trace gas in both the stratosphere and the troposphere. The global monitoring of the
 *	ozone concentration, using both satellite-borne and ground-based instruments, plays a key role in the determination
 *	of the long-term trends for the stratospheric ozone layer, which protects the biosphere from harmful UVB
 *	radiation and air quality related studies.
 *
 *	The requirement to measure small changes in stratospheric and tropospheric ozone places strong demands on the
 *	accuracy of the ozone absorption cross-sections used in retrievals of the spectra delivered by remote sensing
 *	spectrometers. Several satellite spectrometers have been used to measure low-resolution cross-sections
 *	pre-flight (SCIAMACHY, GOME and GOME-2 flight models).1,2,3 These datasets have the great advantage of automatically
 *	incorporating the instrumental slit functions. However, use of these datasets is not straightforward for
 *	other instruments. In the report for the Absorption Cross Sections of Ozone (ACSO) committee, Weber et al.4
 *	consider the impact of cross-section choice on total ozone retrieval applied to GOME, SCIAMACHY and GOME-2
 *	and discuss necessary resolution matching, wavelength shifts and scalings. [The ACSO committee was established
 *	by the World Meteorological Organization and the International Association of Meteorology and Atmospheric Sciences
 *	to review and recommend ozone cross-sections for all the commonly used (both ground-based and satellite)
 *	atmospheric ozone monitoring instruments.]
 *
 *	Among high-resolution datasets, the most important are the so-called data of Bass–Paur5,6 and data
 *	of Malicet, Daumont, Brion et al.7 (and references cited therein). Regardless of the high quality of these data,
 *	they have serious limitations, leaving room for improvement. Both datasets are based on experimental data acquired
 *	at only five temperatures, compelling researchers to use interpolation for other temperatures. In addition, these
 *	datasets only cover the limited UV and vis spectral regions.
 *
 *	More details on ozone cross-sections obtained before 2003 can be found in a comprehensive overview by Orphal.8
 *	Relevant data are available, for example, from the online spectral atlas of gaseous molecules of the
 *	Max-Planck-Institute for Chemistry, Mainz9 or from the ACSO home page.10
 *
 *	The new accurate broadband cross-sections determined in this study have inherent advantages over the previous
 *	datasets to the maximum possible extent. The data were obtained for 11 temperatures down to 195 K. For convenient
 *	use in various current and future projects, the new dataset uniquely combines a broad spectral coverage from 220 nm
 *	to 1000 nm with spectral resolution as high as 0.02 nm. This dataset enables accurate convolution with the slit
 *	functions of all currently relevant ground-based and satellite-based remote sensing ­instruments.
 *
 *	\par Header details from distributed Data Files
 *	
 *		- Source: IUP, MolSpec Lab, Serdyuchenko A., Gorshelev V, Weber M.
 *		- Spectrometer:   Echelle Spectrometer ESA 4000 and Bruker HR 120 FTS 
 *		- Double jacket quartz cell, thermo-insulated, pre-cooler, cryogenic cooling 
 *
 *		- Spectral Resolution(HWHM):
 *			- 0.01 nm below 290 nm},
 *			- 1 cm-1 between 290 nm and 350 nm},
 *			- 0.01 nm between 350 nm and 450 nm},
 *			- 1 cm-1 between 450 nm and 1100 nm,
 *			.
 *		- Grid: interpolated on grid 0.01 nm
 *		- Absolute calibration: using pure ozone pressure.
 *		- Relative systematic uncertainty budget:
 *			- pressure:    2%
 *			- temperature: 1%
 *			- absorption length: < 0.1%
 *			.
 *		- Total relative systematic uncertainty <3%
 * 
 *		- Concatenated spectra parameters:
 *			- Spectral regions: Lightsource stability  Optical density limits:
 *			- 213-290 nm	    0.5%                   0.5-2	
 *			- 290-310 nm 	    2%                     0.1-2
 *			- 310-340 nm	    1%                     0.1-2
 *			- 340-450 nm	    1%                     0.05-1
 *			- 450-750 nm	    0.2%                   0.5-2
 *			- 750-1100nm	    0.2%                   0.001-0.1
*/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_O3_SerdyuchenkoV1 : public skOpticalProperties_UserDefinedAbsorption,
	                                          public skWavelengthToPSF_SerdyuchenkoV1,								// Describes the instrument used to measure this cross-section
											  public skOpticalProperty_AdditionalStateInfo_TemperatureDependent		// Describes what state parameters affetct the cross-section

{
	private:
												skOpticalProperties_O3_SerdyuchenkoV1	( const skOpticalProperties_O3_SerdyuchenkoV1& other );	// Dont allow copy constructor
		skOpticalProperties_O3_SerdyuchenkoV1&	operator =								( const skOpticalProperties_O3_SerdyuchenkoV1& other );	// Dont allow assignment operator

	public:
							skOpticalProperties_O3_SerdyuchenkoV1();
		virtual			   ~skOpticalProperties_O3_SerdyuchenkoV1() override {}
};

/*---------------------------------------------------------------------------
 *					class skOpticalProperties_O3				  2003-11-28*/
/** \ingroup o3skopticalprop
 *	These are the cross-sections used in the OSIRIS level 2 analysis for Saskmart V5.07.
 *	This table is based upon the cross-sections of Bogumil Orphal and Burrows.
 *	This may change in the furture if it is decided that betetr O3 cross-sections exist.
 *
 *	Note that this class sets the QuietWavelength Truncation flag, 
 *	m_quietwavelengthtruncation in the base class skOpticalProperties_UserDefinedAbsorption, 
 *	so it will return a zero cross-section without any error for all wavelengths outside its
 *	tabular range.
 */
/*-------------------------------------------------------------------------*/

class skOpticalProperties_O3_OSIRISRes: public skOpticalProperties_UserDefinedAbsorption
{
	private:
												skOpticalProperties_O3_OSIRISRes  ( const skOpticalProperties_O3_OSIRISRes& other );	// Dont allow copy constructor
		skOpticalProperties_O3_OSIRISRes&		operator =							( const skOpticalProperties_O3_OSIRISRes& other );	// Dont allow assignment operator

	public:
												skOpticalProperties_O3_OSIRISRes();
		virtual								   ~skOpticalProperties_O3_OSIRISRes() override {}
};

#include "no2/no2xsections.h"
#include "so2/so2xsections.h"

/*---------------------------------------------------------------------------
 *					class skOpticalProperties_NO2_Burrows98		          2003-11-28*/
/** \ingroup no2skopticalprop
 *	Tabulates the absorption cross section of NO2 molecules at medium
 *	resolution (0.2-0.33 nm) from 230 nm to 795 nm at 4 temperatures. The
 *	cross-sections were measured by Burrows et al. with the GOME instrument before flight.
 *
 *	\par Spectral resolution
 *		- 231-307 nm, 0.20 nm (FWHM ?)
 *		- 307-316 nm, 0.20 nm
 *		- 311-405 nm, 0.17 nm
 *		- 405-611 nm, 0.29 nm
 *		- 595-794 nm, 0.33 nm
 *
 *	\par Temperature Range
 *	The cross-sections were measured at 4 temperatures covering normal stratospheric and tropospheric ranges. The paper does not
 *	provide any advice on how to interpolate in temperature. We use the standard (linear interpolation) technique provided by the base class.
 *	The spectra were measured at the following four temperatures.
 *		-# 221K
 *		-# 241K
 *		-# 273K
 *		-# 293K
 *
 *	\par References
 *  J. P. Burrows, A. Dehn, B. Deters, S. Himmelmann, A. Richter, S. Voigt, and J. Orphal:
 *  "Atmospheric Remote-Sensing Reference Data from GOME: 1. Temperature-Dependent Absorption Cross Sections of NO2 in the 231-794 nm Range",
 *	Journal of Quantitative Spectroscopy and Radiative Transfer 60, 1025-1031, 1998.
 *
 * \par Paper Abstract
 *	Absorption cross-sections of NO2 between 231-794 nm have been measured in the
 *	221-293K temperature range, using the global ozone monitoring experiment (GOME) flightmodel
 *	(FM) satellite spectrometer. The spectra have a resolution of about 0.2 nm below 400 nm
 *	and of about 0.3 nm above 400 nm. These are the first reference spectra of NO2 covering at the
 *	same time the entire UV-visible-NIR spectral range and a broad range of relevant atmospheric
 *	temperatures. The new absorption cross-sections are important as accurate reference data for
 *	atmospheric remote-sensing of NO2 and other minor trace gases.
*/
/*-------------------------------------------------------------------------*/

class skOpticalProperties_NO2_Burrows98: public skOpticalProperties_UserDefinedAbsorption,
	                                     public skWavelengthToPSF_TableArray,									// Describes the instrument used to measure these cross-sections 
										 public skOpticalProperty_AdditionalStateInfo_TemperatureDependent		// Describes what state parameters affetct the cross-section
	                                   
{
	private:
												skOpticalProperties_NO2_Burrows98	( const skOpticalProperties_NO2_Burrows98& other );	// dont allow copy constructor
		skOpticalProperties_NO2_Burrows98&		operator =							(const skOpticalProperties_NO2_Burrows98& other );	// Dont allow assignment operator

	public:
												skOpticalProperties_NO2_Burrows98	();
		virtual								   ~skOpticalProperties_NO2_Burrows98	() {}
};

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_NO2_OSIRISRes		2009-6-16*/
/** \ingroup no2skopticalprop
 *	Calculates the absorption cross section of NO2 molecules from 230 nm to 795 nm
 *	and 221 K to 293 K. The cross-sections have been reduced to the resolution
 *	of OSIRIS.
 */
/*---------------------------------------------------------------------------*/

class skOpticalProperties_NO2_OSIRISRes: public skOpticalProperties_UserDefinedAbsorption
{
	private:
												skOpticalProperties_NO2_OSIRISRes	( const skOpticalProperties_NO2_OSIRISRes& other );	// dont allow copy constructor
		skOpticalProperties_NO2_OSIRISRes&		operator =							(const skOpticalProperties_NO2_OSIRISRes& other );	// Dont allow assignment operator

	public:
												skOpticalProperties_NO2_OSIRISRes();
		virtual								   ~skOpticalProperties_NO2_OSIRISRes() override {}
};


