#include <map>
#include <functional>


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_BaseEngine		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKEngine_Stub_SO : public ISKEngine_Stub
{
	private:
		SKTRANSO_Engine							m_engine;
		bool									m_modelisinitialized;
		SKTRAN_AtmosphericOpticalState_V21		m_opticalstate;
		SKTRANSO_SpecificationsUser_Legacy		m_specs;
		SKTRAN_LineOfSightArray_V21				m_linesofsight;
		bool									m_isdirty;
//		std::vector<double>						m_getpropertybuffer;
		std::vector<double>						m_wavelen;
		nx2dArray<double>						m_radiance;
		nx2dArray<ISKStokesVector>				m_radiancepolarized;
		size_t									m_numordersofscatter;
		bool									m_updateclimatology;

	private:
		bool									CheckModelNotInitialized				();
		void									MakeScalarSetFunctions					();
		void									MakeVectorSetFunctions					();
		void									MakeGetFunctions						();

	public:
												ISKEngine_Stub_SO						();
		virtual 							   ~ISKEngine_Stub_SO						();
		virtual bool							AddLineOfSight							( double mjd,  const nxVector& observer, const nxVector& lookvector, int* losindex  ) override;
		virtual bool							AddSpecies								( const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty) override ;
		virtual bool							AddEmission								( const EMISSION_HANDLE& species,    ISKEmission_Stub* emission ) override;
		virtual bool							SetAlbedo								( double albedo ) override;
		virtual bool							SetBRDF									( ISKBrdf_Stub* brdf ) override;
		virtual bool							SetPolarizationMode						( int polarizationmode) override;
		virtual bool							SetAtmosphericState						( ISKClimatology_Stub* climatology ) override;
		virtual bool							SetWavelengths							( const double* wavelen, int numwavelen ) override; 
		virtual bool							InitializeModel							() override;
		virtual bool							CalculateRadiance						( const double** radiance, int* numwavelens, int* numlinesofsight) override;
		virtual bool							CalculateStokesVector					( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight) override;
		virtual bool							GetWeightingFunctions					( const double** wf, int* numwavel, int* numlinesofsight, int* numwf ) override;
};

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_BaseEngine		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKEngine_Stub_OCC : public ISKEngine_Stub
{
	private:
		SKOCCULT_OCC_Engine						m_engine;
		SKOCCULT_Specs_User						m_specs;
		SKTRAN_AtmosphericOpticalState_V21		m_opticalstate;
		SKTRAN_LineOfSightArray_V21				m_linesofsight;
		std::vector< std::vector<double> >		m_extinctionbuffer;
		nx2dArray<double>						m_extinction;
		nx2dArray<ISKStokesVector>				m_extinctionpolarized;
		bool									m_isdirty;
		std::vector<double>						m_wavenumber;							 
		bool									m_updateclimatology;
		bool									m_isconfigured;

	private:
		void									MakeScalarSetFunctions();
		void									MakeVectorSetFunctions();

	public:
												ISKEngine_Stub_OCC						();
		virtual 							   ~ISKEngine_Stub_OCC						();
		virtual bool							AddLineOfSight							( double mjd,  const nxVector& observer, const nxVector& lookvector, int* losindex  )override;
		virtual bool							AddSpecies								( const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty)override ;
		virtual bool							AddEmission								( const EMISSION_HANDLE& species,    ISKEmission_Stub* emission ) override;
		virtual bool							SetAlbedo								( double albedo ) override;
		virtual bool							SetBRDF									( ISKBrdf_Stub* brdf ) override;
		virtual bool							SetPolarizationMode						( int polarizationmode) override;
		virtual bool							SetAtmosphericState						( ISKClimatology_Stub* climatology ) override;
		virtual bool							SetWavelengths							( const double* wavenumber, int numwavenumber ) override;
		virtual bool							InitializeModel							() override;
		virtual bool							CalculateRadiance						( const double** radiance, int* numwavelens, int* numlinesofsight)override;
		virtual bool							CalculateStokesVector					( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight) override;
		virtual bool							GetWeightingFunctions					( const double** wf, int* numwavel, int* numlinesofsight, int* numwf ) override;
};

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKEngine_Stub_HR : public ISKEngine_Stub
{
	private:
		SKTRAN_HR_Engine						m_engine;
		SKTRAN_AtmosphericOpticalState_V21		m_opticalstate;
		SKTRAN_HR_Specs_User					m_specs;
		SKTRAN_LineOfSightArray_V21				m_linesofsight;
		size_t									m_numthreads;
		bool									m_isdirty;
		std::vector<double>						m_wavelen;
		nx2dArray<double>						m_radiance;
		nx2dArray<skRTStokesVector>             m_radiancePol;
		nx2dArray<ISKStokesVector>				m_radiancepolarized;
		bool                                    m_storeStokes;
		size_t									m_numordersofscatter;
		bool									m_updateclimatology;
		bool									m_geometryisconfigured;
		nx3dArray<double>					    m_wfbuffer;
		nx2dArray<double>						m_brdfwfbuffer;
		bool									m_issetdiagnostics;
		bool 									m_storeOpticalDepth;
		std::vector<std::vector<double>>        m_cellOpticalDepthBuffer;
		std::vector<std::vector<double>>        m_cellDistanceBuffer;

	private:
		bool									GetBasisHelio					( GEOGRAPHIC_BASIS* basis, size_t losindex);
		bool 									GetBasisGeo 					( GEOGRAPHIC_BASIS* basis, size_t losindex);
		bool									SetInternalPolarizationMode		( int specifier);
		bool									MakeDefaultOpticalState();				// to be removed once climatology interface is done
		bool									MakeScalarSetFunctions();
		bool									MakeVectorSetFunctions();
		bool									MakeObjectSetFunctions();
		bool									MakeVectorGetFunctions();
		bool									MakeScalarGetFunctions();
		bool									MakeStringSetFunctions();
//		bool									GetPropertyScalar( const char* propertyname, double* value );
//		bool									ParseCommandAndIndex( const nxString& input, nxString& command, int& index );
		bool									CheckModelNotInitalized( const char* propertystr) const;

	public:
												ISKEngine_Stub_HR						();
		virtual 							   ~ISKEngine_Stub_HR						()override;
		virtual bool							AddLineOfSight							( double mjd,  const nxVector& observer, const nxVector& lookvector, int* losindex  ) override;
		virtual bool							AddSpecies								( const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty) override;
		virtual bool							AddEmission								( const EMISSION_HANDLE& species, ISKEmission_Stub* emissionobject) override;
		virtual bool							SetAlbedo								( double albedo )override;
		virtual bool							SetBRDF									( ISKBrdf_Stub* brdf ) override;
		virtual bool							SetPolarizationMode						( int polarizationmode) override;
		virtual bool							SetAtmosphericState						( ISKClimatology_Stub* climatology )override;
		virtual bool							SetWavelengths							( const double* wavelen, int numwavelen )override; 
		virtual bool							InitializeModel							();
		virtual bool							CalculateRadiance						( const double** radiance, int* numwavelens, int* numlinesofsight )override;
		virtual bool							CalculateStokesVector					( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight) override;
		virtual bool							GetWeightingFunctions					( const double** wf, int* numwavel, int* numlinesofsight, int* numwf ) override;
//		virtual bool							SetPropertyScalar						( const char* propertyname, double value ) override;							//!< [NOT Supported in Base class]
//		virtual bool							SetPropertyArray						( const char* propertyname, const double* value, int numpoints ) override;	//!< [NOT Supported in Base class]
//		virtual bool							SetPropertyObject						( const char* propertyname, nxUnknown* object ) override;					//!< [NOT Supported in Base class]
//		virtual bool							SetPropertyString						( const char* propertyname, const char* str) override { return false;}
//		virtual bool							GetProperty								( const char* propertyname, const double** value, int* numpoints ) override;
};



/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_MC		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKEngine_Stub_MC : public ISKEngine_Stub
{
	private:
		SKTRAN_Engine_MC_V21                    m_engine;
		SKTRAN_AtmosphericOpticalState_V21		m_opticalstate;
		SKTRAN_Specifications_MC				m_specs;
		SKTRAN_LineOfSightArray_V21				m_linesofsight;
		bool									m_isdirty;
		std::vector<double>						m_wavelen;
		nx2dArray<double>						m_radiance;
		nx2dArray<double>                       m_variance;
		nx2dArray<double>						m_secondary;
		nx2dArray<double>						m_secondaryVariance;
		nx3dArray<double>						m_airMassFactor;
		nx3dArray<double>						m_airMassFactorVariance;
		bool                                    m_storeStokes;
		bool									m_storeAMF;
		bool									m_storeSecondary;
		nx2dArray<skRTStokesVector>             m_radiancePol;
		nx2dArray<ISKStokesVector>				m_radiancepolarized;
		size_t									m_numordersofscatter;
		size_t									m_numthreads;
		bool									m_updateclimatology;
		bool									m_geometryisconfigured;

	
	private:
		bool									GetBasisHelio( GEOGRAPHIC_BASIS* basis, size_t losindex);
		bool									GetBasisGeo( GEOGRAPHIC_BASIS* basis, size_t losindex);
		bool									MakeDefaultOpticalState ( );				// to be removed once climatology interface is done
		bool                                    MakeScalarSetFunctions  ( );
		bool                                    MakeVectorSetFunctions  ( );
		bool                                    MakeVectorGetFunctions  ( );
		bool									MakeStringSetFunctions  ( );

	public:
												ISKEngine_Stub_MC						();
		virtual 							   ~ISKEngine_Stub_MC						() override;
		virtual bool							AddLineOfSight							( double mjd,  const nxVector& observer, const nxVector& lookvector, int* losindex  ) override;
		virtual bool							AddSpecies								( const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty) override;
		virtual bool							AddEmission								( const EMISSION_HANDLE& species,    ISKEmission_Stub* emission ) override;
		virtual bool							SetAlbedo								( double albedo ) override;
		virtual bool							SetBRDF									( ISKBrdf_Stub* brdf ) override;
		virtual bool							SetPolarizationMode						( int polarizationmode) override;
		virtual bool							SetAtmosphericState						( ISKClimatology_Stub* climatology ) override;
		virtual bool							SetWavelengths							( const double* wavelen, int numwavelen ) override; 
		virtual bool							InitializeModel							( ) override;
		virtual bool							CalculateRadiance						( const double** radiance, int* numwavelens, int* numlinesofsight ) override;
		virtual bool							CalculateStokesVector					( const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight) override;
		virtual bool							GetWeightingFunctions					( const double** wf, int* numwavel, int* numlinesofsight, int* numwf ) override;
};


/*-----------------------------------------------------------------------------
*					ISKEngine_Stub_TIR		2018-6-5*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKEngine_Stub_TIR : public ISKEngine_Stub
{
private:
    SKTRAN_TIR_Engine						m_engine;
    SKTRAN_TIR_AtmosphericOpticalState		m_opticalstate;
    SKTRAN_TIR_Specs_User					m_specs;
    SKTRAN_LineOfSightArray_V21				m_linesofsight;
    size_t									m_numthreads;
    bool									m_cacheisvalid;
    std::vector<double>						m_wavelen;
    nx2dArray<double>						m_radiance;
    bool									m_modelisconfigured;
    nx3dArray<double>						m_wfbuffer;

    bool									m_usecache;

    std::map< nxString, std::function<bool(double)> >				m_scalarsetfunctions;
    std::map< nxString, std::function<bool(const double*, int)> >	m_vectorsetfunctions;
    std::map< nxString, std::function<bool(nxUnknown*)> >			m_objectsetfunctions;
    std::map< nxString, std::function<bool(double*)> >				m_scalargetfunctions;
    std::map< nxString, std::function<bool(int)> >					m_vectorgetfunctions;
    std::map< nxString, std::function<bool(const char*)> >			m_stringsetfunctions;
    std::vector<double>												m_getpropertybuffer;

private:
    bool									MakeScalarSetFunctions();
    bool									MakeVectorSetFunctions();
    bool									MakeStringSetFunctions();
    bool									MakeVectorGetFunctions();
    bool									MakeScalarGetFunctions();
    bool									GetPropertyScalar(const char* propertyname, double* value);
    bool									ParseCommandAndIndex(const nxString& input, nxString& command, int& index);
    bool									CheckModelNotInitialized(const char* propertystr) const;

public:
    ISKEngine_Stub_TIR();
    virtual 							   ~ISKEngine_Stub_TIR()override;
    virtual bool							AddLineOfSight(double mjd, const nxVector& observer, const nxVector& lookvector, int* losindex) override;
    virtual bool							AddSpecies(const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty) override;
    virtual bool							AddWeightingFunctionSpecies(CLIMATOLOGY_HANDLE& species);
    virtual bool							AddEmission(const EMISSION_HANDLE& species, ISKEmission_Stub* emissionobject) override;
    virtual bool							SetAlbedo(double albedo)override;
    virtual bool							SetBRDF(ISKBrdf_Stub* brdf) override;
    virtual bool							SetPolarizationMode(int polarizationmode) override;
    virtual bool							SetAtmosphericState(ISKClimatology_Stub* climatology)override;
    virtual bool							SetWavelengths(const double* wavelen, int numwavelen)override;
    virtual bool							InitializeModel();
    virtual bool							CalculateRadiance(const double** radiance, int* numwavelens, int* numlinesofsight)override;
    virtual bool							CalculateStokesVector(const ISKStokesVector** radiancep, int* numwavelens, int* numlinesofsight) override;
    virtual bool							GetWeightingFunctions(const double** wf, int* numwavel, int* numlinesofsight, int* numwf) override;
    virtual bool							SetPropertyScalar(const char* propertyname, double value) override;							//!< [NOT Supported in Base class]
    virtual bool							SetPropertyArray(const char* propertyname, const double* value, int numpoints) override;	//!< [NOT Supported in Base class]
    virtual bool							SetPropertyObject(const char* propertyname, nxUnknown* object) override;					//!< [NOT Supported in Base class]
    virtual bool							SetPropertyString(const char* propertyname, const char* str) override;
    virtual bool							GetProperty(const char* propertyname, const double** value, int* numpoints) override;
};


/*-----------------------------------------------------------------------------
 *					class ISKGeodetic_Stub_std						2014- 5- 7*/
/** A stub for the ISKCgeodetic class
**/
/*---------------------------------------------------------------------------*/

class ISKGeodetic_Stub_std : public ISKGeodetic_Stub
{
	private:
		nxGeodetic		m_geoid;

	private:
		void				MakeStringSetFunctions();

	public:
							ISKGeodetic_Stub_std		( );
		virtual			   ~ISKGeodetic_Stub_std		( ) override;
      	virtual bool		FromGeodetic				( double latitude, double longitude,  double Height ) override;	
      	virtual bool		FromGeocentric				( const nxVector& geocentric  ) override;
		virtual bool		FromTangentPointLocation	( const nxVector& r, const nxVector& lookv ) override;
		virtual bool		FromTangentAltitude			( double required_height, const nxVector& spacecraftlocation, const nxVector& boresightplane, nxVector* requiredlookvector)override;
      	virtual nxVector	GeodeticWest				( ) override;
      	virtual nxVector	GeodeticSouth				( ) override;
      	virtual nxVector	GeodeticUp					( ) override;
		virtual nxVector	Location					( ) override;
		virtual double		GeodeticLongitude			( ) override;
		virtual	double		GeodeticLatitude			( ) override;
		virtual double		Height						( ) override;
		virtual bool		GetShellHeightLocation		( double H, const nxVector& observerposition, const nxVector& look, nxVector* entrypoint, nxVector* exitpoint) override;
		//virtual bool		SetPropertyScalar			( const char* propertyname, double value );	
		virtual nxVector	OsculatingSpheroidCenter    ( ) override;
		virtual double		OsculatingSpheroidRadius    ( ) override;
};

