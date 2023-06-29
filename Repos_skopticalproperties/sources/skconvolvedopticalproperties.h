class skOpticalProperties_ConvolvedDiscreteWavelenEntry;

/*-----------------------------------------------------------------------------
 *					class skconvolvedabsorbtionfuncptr		2012-4-30*/
/** \internal
 *	This is a helper class taht **/
/*---------------------------------------------------------------------------*/

class skconvolvedabsorbtionfuncptr
{
	private:
		skOpticalProperties_ConvolvedDiscreteWavelenEntry*		m_entry;
		double							m_w0;			// The central wavelength of the Guassian distribution in nanometers
		double							m_sd;			// The standard deviation of the Guassian distribution in nanometers
		int								m_xsid;			// cross-section id, 1 = absorption, 2 =scattering, 3 =  extinction, 4 = 1.0( to check normalization)
//		size_t							m_threadindex;

	public:
										skconvolvedabsorbtionfuncptr( );
		static int						XSID_Absorption				() { return 1;}
		static int						XSID_Scatter				() { return 2;}
		static int						XSID_Extinction				() { return 3;}
		static int						XSID_Normalize				() { return 4;}
		bool							Configure					( skOpticalProperties_ConvolvedDiscreteWavelenEntry* entry, double w0, double sd);
		bool							SetVariable					( int xsid )			{m_xsid = xsid; return true;}
		double							operator ()					( double wavelennm);
};

/*-----------------------------------------------------------------------------
 *					class skConvolvedWavelengthEntry				2012-4-30*/
/** \internal
 *	This is a class that generates a new cross-section based upon the
 *	convolution of existing cross-sections. This can be used to generate instrument
 *	specific cross-sections.
 *
 *	This class only manages the cross-section convolution for one wavelength/wavenumber.
 *	At the curruent time the class does not convolve the scattering phase matrix as we are 
 *	we are mostly interested in the effect upon absorbing species.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_ConvolvedDiscreteWavelenEntry
{
	private:
		bool						m_isdirty;					//!< True if the cross section is out of date
		double						m_wavenumber;				//!< The wavenumber (cm-1) of the cross-section
		double						m_psf_fwhm_nm;				//!< The point spread function, FWHM expressed in nanometers
		double						m_highres_stepsize_nm;		//!< Nominal resolution of the high resolution spectrum in nanometers.
		double						m_convolvedabsorption;		//!< Local cache of absorption cross-section, only valid if not dirty
		double						m_convolvedscattering;		//!< Lcoal cache of scattering cross-section, only valid if not dirty
		double						m_convolvedextinction;		//!< Local cache of extinction cross-section, only valid if not dirty
		skOpticalProperties*		m_highresopticalprops;		//!< The optical properties used to generate the convolution
		std::mutex				m_mutex;

	private:
		void						init();
		bool						CheckDirtyAndUpdate				( );
		bool						GenerateConvolvedCrossSections	( );

	public:
									skOpticalProperties_ConvolvedDiscreteWavelenEntry		();
									skOpticalProperties_ConvolvedDiscreteWavelenEntry		( const skOpticalProperties_ConvolvedDiscreteWavelenEntry& other );
									skOpticalProperties_ConvolvedDiscreteWavelenEntry		( double wavenumber );
								   ~skOpticalProperties_ConvolvedDiscreteWavelenEntry		();
		bool						operator <						( const skOpticalProperties_ConvolvedDiscreteWavelenEntry& other ) const		{ return m_wavenumber <  other.m_wavenumber;}
		bool						operator ==						( const skOpticalProperties_ConvolvedDiscreteWavelenEntry& other ) const		{ return m_wavenumber == other.m_wavenumber;}
		void						SetDirty						()														{ m_isdirty = true;}
		skOpticalProperties*		HighResOpticalProperties		()														{ return m_highresopticalprops;}
		bool						SetHighResOpticalProperties		( skOpticalProperties* optprop, double nm_resolution);
		bool						SetWavenumberAndFWHM			( double wavenumber, double fwhm_nm);
		double						ConvolvedAbsorption				( )	{ CheckDirtyAndUpdate( ); return m_convolvedabsorption;}
		double						ConvolvedScattering				( )	{ CheckDirtyAndUpdate( ); return m_convolvedscattering;}
		double						ConvolvedExtinction				( )	{ CheckDirtyAndUpdate( ); return m_convolvedextinction;}
};


/*-----------------------------------------------------------------------------
 *					class skConvolvedWavelengthEntriesTable				2012-4-30*/
/** \internal
 *	This class stores a table of convolved wavelength entry tables. This table
 *	is used by the convolved optical properties class to store the convolved 
 *	cross-section for every wavelength of interest. 
 *	I have puposely configured this table to store unique entries for each
 *	different wavelength as we anticipate that we are normally only using a handful (6-12)
 *	of wavelengths
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable
{
	private:         std::map< double, skOpticalProperties_ConvolvedDiscreteWavelenEntry >						m_entries;
	public:  typedef std::map< double, skOpticalProperties_ConvolvedDiscreteWavelenEntry >::iterator			iterator;
	public:  typedef std::map< double, skOpticalProperties_ConvolvedDiscreteWavelenEntry >::const_iterator		const_iterator;
	public:  typedef std::map< double, skOpticalProperties_ConvolvedDiscreteWavelenEntry >::value_type			value_type;

	private:

	public:
											skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable();
		bool								FindEntry( double wavenumber, skOpticalProperties_ConvolvedDiscreteWavelenEntry** entry  = NULL);
		bool								AddEntry ( double wavenumber, skOpticalProperties* highresptoperties, double nm_resolution, double fwhm_nm );
		bool								SetDirty ();
};


/*-----------------------------------------------------------------------------
 *					class skOpticalProperties_ConvolvedDiscreteWavelenCachedState				2012-4-30 */
/** \ingroup skopticalpropmisc
 *	A class that convolves a high resolution cross-section spectrum to a
 *	requested spectral resolution specified as the FWHM in nanometers. The class caches convolved 
 *	cross-section data as a function of both atmospheric state and wavelength. This allows the class to
 *	work quite efficiently in Sasktran RT calculations where the atmopsheric state is rapidly changed as we
 *	process rays at different locations.
 *
 *	Most (all the ones we have) high resolution cross-sections are only a function of a few atmospheric parameters 
 *	(eg. Temperature and Pressure) and we use class skOpticalProperty_AdditionalStateInfo and its descendants to convert
 *	atmospheric state parameters into a unique key. This class was the next generation of development after
 *	older class #skOpticalProperties_Convolved.
 *
 *	This class compensates for the spectral resolution of the source, high resolution cross-sections
 *	by using a #skWavelengthToPSF_Table object. This object, which is automatically built into most of our
 *	measured high resolution cross-section tables, provides the spectral resolution of the source cross-sections
 *	as a functionof wavelength.
 *
 *	\par Caching of convolved Cross-sections
 *	This class caches cross-sections for each new wavelength as well as atmospheric state. This generally improves performance as it will
 *	only execute the convolution if atmospheric state parameters (eg temperature) change after a call to #SetAtmosphericState. 
 *	The code works well if you have a limited set of discrete wavelengths (eg. less than a few thousand) but there are no hard
 *	limits on the maximum number. Note that our high resolution cross-sections, which are suitable for convolution, have the
 *	appropriate #skWavelengthToPSF_Table and #skOpticalProperty_AdditionalStateInfo interfaces built into their respective
 *	classes through multiple inheritance.
 * 
 *	\par Main Details
 *	The class has two functions that differentiate it from other optical property classes,
 *		-# #SetHighResOpticalProperties.
 *		-# #GetTargetFWHM
 *
 *  \par An example
 *	The Convolved cross-section optical properties are very similar to all other optical properties but here is
 *	small example showing how you set up the code to convolve the SCHIAMACHY BOGUMIL V4 O3 cross section to a fixed FWHM convolution.
 *	Note the class automatically adjusts for the changing resolution of the BOGUMIL cross-sections.
 *
 *	\code
 * void example()
 * {
 *	skOpticalProperties_O3_SciaBogumilV4*						o3sciav4;       // The high resolution cross-section
 *	skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM*		o3conv;	        // The convolving optical properties
 *
 *	o3sciav4 = new skOpticalProperties_O3_SciaBogumilV4;                        // Create the high resolution optical properties
 *	o3conv   = new skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM;       // Create the convolving class
 *	o3sciav4->AddRef();                                                         // manage the lifetime of the high res cross-section
 *	o3conv  ->AddRef();                                                         // manage the lifetime of the convolving class
 *	o3conv  ->SetHighResOpticalProperties( o3sciav4, o3sciav4, o3sciav4 );      // configure the convolver so it uses the Sciamachy high res cross-section
 *	o3conv  ->SetFWHM( 0.9 );                                                   // configure the convolver so it convolves to a spectral resolution of 0.9 nm
 *
 *	MakeCrossectionVersusTemperatureArray( o3conv,   "o3conv.txt");             // Use the convolving optical properties using standard skOpticalProperties intwerface.
 *
 *	o3sciav4->Release();                                                        // Release the high res spectral object
 *	o3conv->Release();                                                          // Release the convolver object.
 * }
 *	\endcode
 *	\par See Also
 *	#skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM	
 *
*/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_ConvolvedDiscreteWavelenCachedState : public skOpticalProperties
{

	private:
				 std::map<skOpticalProperty_AdditionalStateInfoKey, skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable>						m_cachedentriestable;			// The cached convolved cross-sections for all cached atmospheric states (eg pressure temperature).
		typedef  std::map<skOpticalProperty_AdditionalStateInfoKey, skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable>::iterator				iterator;
		typedef  std::map<skOpticalProperty_AdditionalStateInfoKey, skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable>::const_iterator		const_iterator;
		typedef  std::map<skOpticalProperty_AdditionalStateInfoKey, skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable>::value_type			value_type;

	private:
		skClimatology*																m_backgroundatmosphere;
		skOpticalProperties*														m_highresopticalproperties;
		skWavelengthToPSF_Table*													m_highresdetails;				// The high resolution point spread function and spacing
		skOpticalProperty_AdditionalStateInfo*										m_atmosphericstateinfo;			// The object that gerenates indexes for m_cachedentriestable from the atmospheric state
		skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable*					m_entriestable;					// The current convolved cross-sections for the current atmospheric state (eg pressure temperature).

	private:
		bool															CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs ) const;
																		skOpticalProperties_ConvolvedDiscreteWavelenCachedState		( const skOpticalProperties_ConvolvedDiscreteWavelenCachedState& other );	// Dont allow copy copnstructor
		skOpticalProperties_ConvolvedDiscreteWavelenCachedState& 		operator =													( const skOpticalProperties_ConvolvedDiscreteWavelenCachedState& other);	// Dont allow assignment constructor

	public:
																		skOpticalProperties_ConvolvedDiscreteWavelenCachedState		();
		virtual														   ~skOpticalProperties_ConvolvedDiscreteWavelenCachedState		();
		/*-----------------------------------------------------------------------------
		 *					SetHighResOpticalProperties								2011-8-9*/
		/** Function to define the source, high resolution cross-section data. Users will normal call this
		 *	to define the high resolution source cross-section data to convolve. We have written many of the
		 *	high resolution cross-section classes so that one instance of the class can be both the first
		 *	and second parameter; in fact we recommend this technique as we dont do lifetime management on the
		 *	#measurementdetails object.
		 *
		 *	\param [in] highresopticalproperties
		 *	Provides access to the high resolution cross-section data. The lifetime of the object
		 *	is managed by the class.
		 *
		 *	\param [out] measurementdetails
		 *	Pointer to an object that can be used to get the spectral resolution of the source cross-section data
		 *	as a function of wavelength. Note that the lifetime of this object is not managed. This is not that big an issue
		 *	as it is normally the same object as #highresopticalproperties
		 *
		 *	\returns
		 *	True if successful.
		 **/
		/*---------------------------------------------------------------------------*/
		bool								SetHighResOpticalProperties	 ( skOpticalProperties*						highresopticalproperties,
																		   skWavelengthToPSF_Table*					measurementdetails,
																		   skOpticalProperty_AdditionalStateInfo*	atmosphericstateinfo);


	public:

		/*-----------------------------------------------------------------------------
		 *					GetTargetFWHM								2011-8-9*/
		/** A purely virtual function that fetches
		 *	the desired spectral resolution of the target convolution centered at the
		 *	wavelength of interest. End users usually dont call this function as it
		 *	automatically called when doing the convolution calculation. This function must be
		 *	implemented in a derived class.
		 *
		 *	\param [in] wavelen_nm
		 *	The wavelength at the which the spectral resolution is required 
		 *
		 *	\param [out] fwhm_nm
		 *	Returns the spectral resolution at the desired wavelength as a FWHM in nanmeters.
		 *
		 *	\returns
		 *	True if successful.
		 **/
		/*---------------------------------------------------------------------------*/
		virtual bool						GetTargetFWHM						( double wavelen_nm, double* fwhm_nm ) const = 0;

	// -- virtuals derived from base class
	public:									
		virtual bool						SetAtmosphericState					( skClimatology* neutralatmosphere)  override;
		virtual bool						SetLocation							( const GEODETIC_INSTANT& pt, bool* crosssectionschanged )  override;
		virtual bool						InternalClimatology_UpdateCache			( const GEODETIC_INSTANT& pt)  override;
		virtual bool						CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs  ) override;
		virtual bool						CalculatePhaseMatrix				( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix )  override;
		virtual bool						IsScatterer							() const  override;
		virtual bool						IsAbsorber							() const  override;
		virtual double						DeltaFunctionForwardScatterFraction	() const  override;
};


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM		2012-8-2*/
/** \ingroup skopticalpropmisc
 *	An implementation of #skOpticalProperties_ConvolvedDiscreteWavelenCachedState that uses
 *	constant value for the target convolution. See #SetFWHM
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM : public skOpticalProperties_ConvolvedDiscreteWavelenCachedState
{
	private:
		double	m_fwhm;

	private:
					 skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM( const skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM& other ); // Dont allow copy constructor
		skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM& operator = ( const skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM& other );	// Dont allow assignment operator
	public:
									skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM() { m_fwhm = 0.9;}
		virtual					   ~skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM() {}
		/*-----------------------------------------------------------------------------
		 *					SetFWHM								2011-8-9*/
		/** set the desired spectral resolution of the final convolution expressed
		 *	as FWHM in nanometers.
		 *
		 *	\param [in] fwhm
		 *	The desired spectral resolution, FWHM in nanometers.
		 **/
		/*---------------------------------------------------------------------------*/
		void						SetFWHM			(double fwhm )											{ m_fwhm  = fwhm;}
		virtual bool				GetTargetFWHM	( double /*wavelen_nm*/, double* fwhm_nm ) const override 	{ *fwhm_nm = m_fwhm; return true;}
};


