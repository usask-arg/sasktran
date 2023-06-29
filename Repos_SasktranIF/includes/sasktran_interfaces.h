
/**** DO NOT CHANGE VARIABLE NAMES IN THIS FILE AS SWIG RELIES UPON EXACT PATTERN MATCHING OF VARIABLE NAMES ****/

extern "C" bool SKTRAN_IFSetRegistryDirectory(const char* registrydirname);							// An internal function that allows Python SasktranIF to set the registry directory for all SasktranIF component DLL's
extern "C" bool SKTRAN_IFCreateRegistryEntriesForDLL(const char* dllname, const char* paramstr);	// An internal function that allows Python SasktranIF to load a compnent DLL and ask it create registry entries.


/*-----------------------------------------------------------------------------
 *					ISKModuleBase		 2015- 8- 26*/
/** The base class for all the Python/Sasktrn objects exported by this
 *  class.
 **/
/*---------------------------------------------------------------------------*/

class ISKModuleBase
{
	private:
		char*					m_dllname;

	protected:
		virtual bool			SetPropertyScalar		( const char* propertyname, double value ) = 0;
		virtual bool			SetPropertyArray		( const char* propertyname, const double* value, int numpoints ) = 0;
		virtual bool			SetPropertyObject		( const char* propertyname, ISKModuleBase* object ) = 0;
		virtual bool			SetPropertyString		( const char* propertyname, const char* str ) = 0;
		char**					DllNamePtr				() { return &m_dllname;}

	public:
								ISKModuleBase			();
		virtual				   ~ISKModuleBase			();

		virtual nxUnknown*		RawObjectUnknown		() = 0;
		bool					SetProperty				( const char* propertyname, void* valueorobject, int numpoints_or_type  );
		virtual bool			GetProperty				( const char* propertyname, const double** propertyvalue, int* numpoints );
};

/*-----------------------------------------------------------------------------
 *					ISkClimatology		2014-2-8*/
/** This is the class that users use in their program to access the
 *	skClimatology objects via a simplified interface.
 **/
/*---------------------------------------------------------------------------*/

class ISKClimatology : public ISKModuleBase
{
	private:
		ISKClimatology_Stub*	m_climatology;

	protected:
		virtual nxUnknown*		RawObjectUnknown			();
		virtual bool			SetPropertyScalar			( const char* propertyname, double value ) override;
		virtual bool			SetPropertyArray			( const char* propertyname, const double* value, int numpoints ) override;
		virtual bool			SetPropertyObject			( const char* propertyname, ISKModuleBase* object ) override;
		virtual bool			SetPropertyString			( const char* propertyname, const char* str ) override;

	public:
								ISKClimatology				( const char* climatologyname);
		virtual				   ~ISKClimatology				();
		ISKClimatology_Stub*	Stub						() { return m_climatology;}
		bool					Create_New_ClimatologyName	( const char* name );
		bool					IsValidObject				() const { return m_climatology != nullptr;}
		bool					UpdateCache					( const GEODETIC_INSTANT& location );
		bool					GetParameter				( const char * climatology_handle_name,  const GEODETIC_INSTANT& location, double* valueout );
		bool					GetHeightProfile			( const char * climatology_handle_name,  const GEODETIC_INSTANT& location, const double* altitude, double *profile, int numalts );
		bool					SetPropertyUserDefined		( const char * climatology_handle_name,  double* profilevalues, int numpoints);
};


/*-----------------------------------------------------------------------------
 *					ISkOpticalProperty		2014-2-8*/
/** This is the class that users use in their program to access the
 *	skOpticalProperty objects via a simplified interface.
 **/
/*---------------------------------------------------------------------------*/

class ISKOpticalProperty: public ISKModuleBase
{
	private:
		ISKOpticalProperty_Stub*	m_optprop;

	protected:
		virtual nxUnknown*		RawObjectUnknown					();
		virtual bool			SetPropertyScalar					( const char* propertyname, double value ) override;
		virtual bool			SetPropertyArray					( const char* propertyname, const double* value, int numpoints ) override;
		virtual bool			SetPropertyObject					( const char* propertyname, ISKModuleBase* object ) override;
		virtual bool			SetPropertyString					( const char* propertyname, const char* str ) override;

	public:
								ISKOpticalProperty					();
								ISKOpticalProperty					( const char* optpropname);
		virtual				   ~ISKOpticalProperty					();
		ISKOpticalProperty_Stub* Stub								() { return m_optprop;}
		bool					IsValidObject						() const { return m_optprop != nullptr;}
		bool					SetAtmosphericState					( ISKClimatology& atmosphere);
		bool					SetLocation							( const GEODETIC_INSTANT& pt );
		bool					InternalClimatology_UpdateCache		( const GEODETIC_INSTANT& pt);
		bool					CalculateCrossSections				( const double * wavenumber,  double *absxs, double* extxs, double* scattxs, int numortype);
		bool					CalculatePhaseMatrix				( double wavenumber, double cosscatterangle, double phasematrix[16] );
		bool					AddUserDefined						( double temperature,  double* wavelen_nm, int numwave, double* crosssection, int numcross);
		bool                    AddUserDefinedPressure              ( double* pressure, int numpressure, double* temperature, int numtemperature, double* wavelen_nm, int numwavel, double* crosssection, int numcross, double broadnervmr );

};


/*-----------------------------------------------------------------------------
 *					ISKEmission		2014-2-8*/
/** This is the class that users use in their program to access the
 *	skEmission objects via a simplified interface.
 **/
/*---------------------------------------------------------------------------*/

class ISKEmission: public ISKModuleBase	
{
	private:
		ISKEmission_Stub*		m_emissionprop;

	protected:
		virtual nxUnknown*		RawObjectUnknown					();
		virtual bool			SetPropertyScalar					( const char* propertyname, double value ) override;
		virtual bool			SetPropertyArray					( const char* propertyname, const double* value, int numpoints ) override;
		virtual bool			SetPropertyObject					( const char* propertyname, ISKModuleBase* object ) override;
		virtual bool			SetPropertyString					( const char* propertyname, const char* str ) override;



	public:
								ISKEmission							( const char* emissionname);
		virtual				   ~ISKEmission							();
		ISKEmission_Stub*		Stub								() { return m_emissionprop;}
		bool					IsValidObject						() const { return m_emissionprop != nullptr;}
		bool					UpdateLocation						( const GEODETIC_INSTANT& pt, bool isground );
		bool					UpdateCache							( const GEODETIC_INSTANT& pt);
		bool					IsotropicEmission					( const double * wavenumber, double * isotropicradiance, int numorscalar);	
};


/*-----------------------------------------------------------------------------
 *					ISKBrdf										2016-12-12*/
/** This is the class that users use in their program to access the
 *	SKBRDF objects via a simplified interface.
 **/
/*---------------------------------------------------------------------------*/

class ISKBrdf: public ISKModuleBase	
{
	private:
		ISKBrdf_Stub*			m_brdfprop;

	protected:
		virtual nxUnknown*		RawObjectUnknown					();
		virtual bool			SetPropertyScalar					( const char* propertyname, double value ) override;
		virtual bool			SetPropertyArray					( const char* propertyname, const double* value, int numpoints ) override;
		virtual bool			SetPropertyObject					( const char* propertyname, ISKModuleBase* object ) override;
		virtual bool			SetPropertyString					( const char* propertyname, const char* str ) override;



	public:
								ISKBrdf								( const char* brdfname);
		virtual				   ~ISKBrdf								();
		ISKBrdf_Stub*			Stub								() { return m_brdfprop;}
		bool					IsValidObject						() const { return m_brdfprop != nullptr;}
		bool					BRDF								( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* return_brdf);
};

/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum		 2015- 1- 7*/
/** This is a class **/
/*---------------------------------------------------------------------------*/

class ISKSolarSpectrum: public ISKModuleBase
{
	private:
		ISKSolarSpectrum_Stub*	m_solarspectrum;

	protected:
		virtual nxUnknown*		RawObjectUnknown					();
		virtual bool			SetPropertyScalar					( const char* propertyname, double value ) override;
		virtual bool			SetPropertyArray					( const char* propertyname, const double* value, int numpoints ) override;
		virtual bool			SetPropertyObject					( const char* propertyname, ISKModuleBase* object ) override;
		virtual bool			SetPropertyString					( const char* propertyname, const char* str ) override;


	public:
								ISKSolarSpectrum					( const char* solarspectrumname);
		virtual				   ~ISKSolarSpectrum					();
		bool					IsValidObject						() const { return m_solarspectrum != nullptr;}
		bool					Irradiance							( const double* wavelen_nm_vacuum, double* irradiance, int numpoints );
		bool					IrradianceAt1AU						( const double* wavelen_nm_vacuum, double* irradiance, int numpoints );
		bool					SetSolarDistanceFromMjd				( double mjd );		
		bool					SetSolarDistanceFromAU				( double au );		
		bool					MinValidWavelength					( double* minwavelength_nm);
		bool					MaxValidWavelength					( double* maxwavelength_nm);
		bool					NanometerResolutionFWHM				( const double* wavelen_nm_vacuum, double* resolution_nm_fwhm, int numpoints) ;
		bool					SampleSpacing						( const double* wavelen_nm_vacuum, double* sample_spacing, int numpoints    ) ;
};





/*-----------------------------------------------------------------------------
 *					class ISKStokesVectorIF	 2015- 11- 24*/
/** A wrapper class for the ISKSTokesVector object. This class is only
  * used in the MATLAB installations.
 **/
/*---------------------------------------------------------------------------*/

class ISKStokesVectorIF: public ISKModuleBase
{
	public:
		ISKStokesVector			m_stokes;

	protected:
		virtual bool			SetPropertyScalar		( const char* /*propertyname*/, double /*value*/ ) override { return false;}
		virtual bool			SetPropertyArray		( const char* /*propertyname*/, const double* /*value*/, int /*numpoints*/ ) override { return false;}
		virtual bool			SetPropertyObject		( const char* /*propertyname*/, ISKModuleBase* /*object*/ ) override  {return false;}
		virtual bool			SetPropertyString		( const char* /*propertyname*/, const char* /*str*/ ) override {return false;}

	public:

	public:
		virtual				   ~ISKStokesVectorIF(){};
		virtual nxUnknown*		RawObjectUnknown		() { return NULL;}
};
/*-----------------------------------------------------------------------------
 *					ISKEngine										2014-2-8*/
/** \par Purpose
 *	The ISKEngine exposes an interface to one of the Sasktran radiative transfer engines.
 *	The purpose of the engine is to calculate atmospheric radiance for
 *	observers either inside or outside the atmosphere looking in any direction.
 *	The ISKEngine class is configured in its contructor to use the engine of choice. The user will
 *	typically create an instance ISKEngine, define the lines of sight and the
 *	atmospheric optical state and then calculate the radiance seen along the lines of sight.
 *	Other optional configuration steps may be performed.
 *
 *	\par I just want to run the code
 *	The actual sasktran engines that do all of the radiative transfer work are complex pieces
 *	of code and we found that many users who want to just run the code find
 *	the level of complexity too intimidating to be useful. These classes are for
 *	those users who just want to run the code and dont want to mess around building
 *	and compiling vast swaths of code.
 *
 *	\par Plug and Play
 *  The purpose of the ISKEngine class is to develop a
 *	simple "plug and play" interface for people who want to use radiative
 *	transfer codes developed by others and dont want to be bogged down in
 *	the configuration details. We have intentionally hidden the details
 *	from the end user. The engine developers are responsible for providing
 *	good, out-of-the-box performance. The developers have pre-compiled
 *	all of the complex engine code into pre-built dynamic link libraries (DLL)
 *	on Windows. Windows users can simply download and install a new engine
 *	on their system. Unfortunately Linux users must still build the code
 *	from scratch and install a shareable object.
 *
 **/
/*---------------------------------------------------------------------------*/

class ISKEngine: public ISKModuleBase
{
	private:
		ISKEngine_Stub*			m_engine;																//!< Pointer to the class that provides the virtualization to all the engines

	private:
								ISKEngine					( const ISKEngine& /*other*/ ) {throw("ISKEngine Copy Constructor not allowed");}					//!< Users cannot copy this object.
		ISKEngine&				operator =					( const ISKEngine& /*other*/ ) {throw("ISKEngine assignment operator not allowed"); return *this;}	//!< Users cannot assign this obejct

	protected:
		virtual nxUnknown*		RawObjectUnknown			() { return m_engine->RawObjectPointer();}
		virtual bool			SetPropertyScalar			( const char* propertyname, double value ) override;
		virtual bool			SetPropertyArray			( const char* propertyname, const double* value, int numpoints ) override;
		virtual bool			SetPropertyObject			( const char* propertyname, ISKModuleBase* object ) override;
		virtual bool			SetPropertyString			( const char* propertyname, const char* str ) override;


	public:
								ISKEngine					( const char* enginename );
		virtual				   ~ISKEngine					();
		ISKEngine_Stub*			Stub						() { return m_engine;}
		bool					IsValidObject				() const { return m_engine != nullptr;}
		bool					AddLineOfSight				( double mjd,  const nxVector& observer, const nxVector& lookvector, int* losindex  );					//!< [Mandatory]
		bool					AddSpecies					( const char* climatology_handle_name, ISKClimatology& climatology, ISKOpticalProperty& opticalproperty); //!< [Mandatory], note climatology_handle_name is SWIG specific]
		bool					AddEmission					( const char* climatology_handle_name,   ISKEmission&    emission );									  //!< [Optional], note climatology_handle_name is SWIG specific]
		bool					SetAtmosphericState			( ISKClimatology& climatology );									//!< [Optional]
		bool					SetAlbedo					( double albedo );													//!< [Mandatory for most engines]
		bool					SetBRDF						( ISKBrdf* brdf );
		bool					SetPolarizationMode			( int polarizationmode);
		bool					SetWavelengths				( const double* wavelen, int numwavelen );							//!< [Mandatory]
		bool					InitializeModel				();																	//!< [Mandatory]
		bool					CalculateRadiance			( const double**          radiance,  int* numwavelens, int* numlinesofsight);				//!< [Mandatory]
		bool					CalculateStokesVector		( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight);				//!< [Mandatory but not fully Supported by all engines]
		bool					GetWeightingFunctions		( const double**          wf,        int* numwavelens, int* numlinesofsight, int* numwf );	//!< [Not Supported by all engines]
		virtual bool			GetProperty					( const char* propertyname, const double** propertyvalue, int* numpoints ) override;				//!< [Optional]
};



/*-----------------------------------------------------------------------------
 *					class ISKGeodetic								2014- 5- 6*/
/** A helper class for the SasktranIF code. The class provides a useful
 *	software interface for Geodetic coordinate calculations. It is very useful
 *	for specifiying observer positions and line of sight directions in 
 *	local and global coordinates.
 **/
/*---------------------------------------------------------------------------*/

class ISKGeodetic: public ISKModuleBase
{
	private:
		ISKGeodetic_Stub*						m_geoid;															//!< Pointer to the class that provides the virtualization to all the engines

	protected:
		virtual nxUnknown*		RawObjectUnknown				() { return m_geoid;}
		virtual bool			SetPropertyScalar				( const char* propertyname, double value ) override;
		virtual bool			SetPropertyArray				( const char* propertyname, const double* value, int numpoints ) override;
		virtual bool			SetPropertyObject				( const char* propertyname, ISKModuleBase* object ) override;
		virtual bool			SetPropertyString				( const char* propertyname, const char* str ) override;


   	public:	
	     						ISKGeodetic						();
		virtual				   ~ISKGeodetic						();
		bool					IsValidObject					() const { return m_geoid != nullptr;}
      	bool					SetLocationLatLonAlt			( double latitude, double longitude,  double alt );				//!< Set the current location using the specified geodetic coordinates
      	bool					SetLocationXYZ					( const nxVector& geocentric  );									//!< Set the current location from the specified geocentric X,Y,Z vector. (all in meters).
		bool					SetLocationFromTangentPoint		( const nxVector& r, const nxVector& lookv );						//!< Set the current location from the implied tangent point
		bool					SetLocationFromTangentAltitude	( double requiredheight, const nxVector& spacecraftlocation, const nxVector& boresightplane, nxVector* requiredlookvector);	//!< Set the current location from the tangent point at the specified height. Also return the \e look \e vector necessary to do this.
      	nxVector				GetLocalWest					( );			//!< Get the topocentric unit vectors at the current location
      	nxVector				GetLocalSouth					( );			//!< Get the topocentric unit vectors at the current location
      	nxVector				GetLocalUp						( );			//!< Get the topocentric unit vectors at the current location
		nxVector				GetLocationXYZ					( );
		double					GetLongitude					( );
		double					GetLatitude						( );
		double					GetAlt							( );
		bool					GetAltitudeIntercepts			( double H, const nxVector& observerposition, const nxVector& look, nxVector* entrypoint, nxVector* exitpoint);
		nxVector				GetOsculatingSpheroidCenter		( );
		double					GetOsculatingSpheroidRadius		( );
};

class ISKMie : public ISKModuleBase
{
private:
	ISKMie_Stub* m_mie;

protected:
	virtual nxUnknown* RawObjectUnknown() { return m_mie; }
	virtual bool			SetPropertyScalar(const char* propertyname, double value) override;
	virtual bool			SetPropertyArray(const char* propertyname, const double* value, int numpoints) override;
	virtual bool			SetPropertyObject(const char* propertyname, ISKModuleBase* object) override;
	virtual bool			SetPropertyString(const char* propertyname, const char* str) override;

public:
	ISKMie(const char* name);
	virtual ~ISKMie();

	virtual bool Calculate(double lambda, double radius, double refrac_real, double refrac_imag);

	// Output functions
	double								Qext();
	double								Qsca();
	double								Qabs();
	double								Cext();
	double								Csca();
	double								Cabs();

	void								S1(std::complex<double>** s, int* numpoints);
	void								S2(std::complex<double>** s, int* numpoints);
	void								PMom(double** pmom, int* numlegendre);

	void								Sforward(double* real, double* imag);
	void								SBackward(double * real, double* imag);
	void								TForward(int i, double* real, double* imag);
	void								TBackward(int i, double* real, double* imag);
	double								Spike();
};
