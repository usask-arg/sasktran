
/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_BaseEngine		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKOpticalProperty_Stub_Base : public ISKOpticalProperty_Stub
{
	private:
		skOpticalProperties*				m_opticalproperty;

	protected:
	public:
											ISKOpticalProperty_Stub_Base	( skOpticalProperties* optprop);
		virtual 						   ~ISKOpticalProperty_Stub_Base	() override;
		virtual nxUnknown*					RawObjectPointer				() override { return m_opticalproperty;}
		virtual bool						SetAtmosphericState				( ISKClimatology_Stub* atmosphere) override;
		virtual bool						SetLocation						( const GEODETIC_INSTANT& pt ) override;
		virtual bool						InternalClimatology_UpdateCache	( const GEODETIC_INSTANT& pt ) override;
		virtual bool						CalculateCrossSections			( double wavenumber,		                    double* absxs,  double* extxs, double* scattxs) override;
		virtual bool						CalculateCrossSectionsArray		( const double * wavenumber, int numwavenumber, double *absxs,  double* extxs, double* scattxs) override;
		virtual bool						CalculatePhaseMatrix			( const double * wavenumber, const double * cosscatterangle, double * phasematrix ) override;
		virtual bool						AddUserDefined					( double temperature,  double* wavelen_nm, int numwave, double* crosssection, int numcross) override;
		virtual bool                        AddUserDefinedPressure(double* pressure, int numpressure, double* temperature, int numtemperature, double* wavelen_nm, int numwavel, double* crosssection, int numcross, double broadnervmr) override;

};


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Hitran		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKOpticalProperty_Stub_Hitran : public ISKOpticalProperty_Stub_Base
{	
	private:
	skOpticalProperties_HitranChemical*		m_hitranoptprop;

	private:
		void 								MakeSetPropertyFunctions();

	public:
											ISKOpticalProperty_Stub_Hitran	( skOpticalProperties_HitranChemical* optprop);
		virtual							   ~ISKOpticalProperty_Stub_Hitran	() override;
};



/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Aerosol		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKOpticalProperty_Stub_Aerosol : public ISKOpticalProperty_Stub_Base
{
	private:
		skOpticalProperties_AerosolProfile*	m_aerosol_optprop;

		std::vector<double>					m_refrac_wavel;
		std::vector<double>					m_refrac_data;

	private:
		void								MakeSetPropertyFunctions();

	public:
											ISKOpticalProperty_Stub_Aerosol	( skOpticalProperties_AerosolProfile* aerosol_optprop);
		virtual							   ~ISKOpticalProperty_Stub_Aerosol	() override;
};


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_Baum		 2016- 9- 23*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKOpticalProperty_Stub_Baum : public ISKOpticalProperty_Stub_Base
{
	private:
		skOpticalProperties_BaumIceCrystals2014*		m_baum_optprop;
		skClimatology_Constant*							m_effectivesize;

	private:
		void								MakeSetPropertyFunctions();

	public:
											ISKOpticalProperty_Stub_Baum( skOpticalProperties_BaumIceCrystals2014* baum );
										   ~ISKOpticalProperty_Stub_Baum () override;
};


/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_UserDefined	2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKOpticalProperty_Stub_UserDefined : public ISKOpticalProperty_Stub_Base
{
	private:
		skOpticalProperties_UserDefinedAbsorption* m_useroptprop;

	private:
		void								MakeSetPropertyFunctions();

	public:
											ISKOpticalProperty_Stub_UserDefined	( skOpticalProperties_UserDefinedAbsorption* useroptprop);
		virtual							   ~ISKOpticalProperty_Stub_UserDefined	() override;
		virtual bool						AddUserDefined					    ( double temperature,  double* wavelen_nm, int numwave, double* crosssection, int numcross) override;

};

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_UserDefinedPressure	2021-08-06*/
 /** **/
 /*---------------------------------------------------------------------------*/

class ISKOpticalProperty_Stub_UserDefinedPressure : public ISKOpticalProperty_Stub_Base
{
	private:
		skOpticalProperties_UserDefinedAbsorptionPressure* m_useroptprop;

private:
	void									MakeSetPropertyFunctions();


	public:
		ISKOpticalProperty_Stub_UserDefinedPressure(skOpticalProperties_UserDefinedAbsorptionPressure* useroptprop);
		virtual							   ~ISKOpticalProperty_Stub_UserDefinedPressure() override;
		virtual bool                        AddUserDefinedPressure(double* pressure, int numpressure, double* temperature, int numtemperature, double* wavelen_nm, int numwavel, double* crosssection, int numcross, double broadnervmr) override;

};

class ISKOpticalProperty_Stub_UserDefinedScatterConstantHeight : public ISKOpticalProperty_Stub_Base
{
private:
	skOpticalProperties_UserDefinedScatterConstantHeight* m_useroptprop;
	size_t m_numwavelengths;
	
	void MakeSetPropertyFunctions();

public:
	ISKOpticalProperty_Stub_UserDefinedScatterConstantHeight(skOpticalProperties_UserDefinedScatterConstantHeight* useroptprop);
	virtual ~ISKOpticalProperty_Stub_UserDefinedScatterConstantHeight() override;
};

/*-----------------------------------------------------------------------------
 *					ISKOpticalProperty_Stub_ConvolvedFixedFWHM		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKOpticalProperty_Stub_ConvolvedFixedFWHM : public ISKOpticalProperty_Stub_Base
{
	private:
		skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM*	m_convolvedoptprop;
		skWavelengthToPSF_TableConstant										m_default_measurementdetails;
		skOpticalProperty_AdditionalStateInfo_NotDependent					m_default_atmosphericstateinfo;
		skClimatology_OneTemperatureAndPressure								m_defaultclimatology;


	private:
		bool								SetHighResProperties						( nxUnknown* userobject);
		void								MakeSetPropertyFunctions					();

	public:
											ISKOpticalProperty_Stub_ConvolvedFixedFWHM	( skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM* useroptprop);
		virtual							   ~ISKOpticalProperty_Stub_ConvolvedFixedFWHM	() override;
};


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKSolarSpectrum_Stub_Base : public ISKSolarSpectrum_Stub
{
	protected:
		skSolarSpectrum*			m_solar;

	public:
								ISKSolarSpectrum_Stub_Base			( skSolarSpectrum* solar);
		virtual 			   ~ISKSolarSpectrum_Stub_Base			();
		virtual nxUnknown*		RawObjectPointer					() override { return m_solar;}
		virtual bool			Irradiance							( double wavelen_nm_vacuum, double* irradiance ) override;
		virtual bool			IrradianceArray						( const double* wavelen_nm_vacuum, double* irradiance, int numpoints) override;
		virtual bool			IrradianceAt1AU						( double        wavelen_nm_vacuum, double* irradiance ) override;
		virtual bool			IrradianceAt1AUArray				( const double* wavelen_nm_vacuum, double* irradiance, int numpoints )override;
		virtual bool			NanometerResolutionFWHM				( double        wavelen_nm_vacuum, double* resolution_nm_fwhm) override;
		virtual bool			NanometerResolutionFWHMArray		( const double* wavelen_nm_vacuum, double* resolution_nm_fwhm, int numpoints) override;
		virtual bool			SampleSpacing						( double wavelen_nm_vacuum, double* sample_spacing) override;
		virtual bool			SampleSpacingArray					( const double* wavelen_nm_vacuum, double* resolution_nm_fwhm, int numpoints)override;
		virtual bool			SetSolarDistanceFromMjd				( double mjd ) override;	
		virtual bool			SetSolarDistanceFromAU				( double au ) override;		
		virtual bool			MinValidWavelength					( double* minwavelength_nm) override;
		virtual bool			MaxValidWavelength					( double* minwavelength_nm) override;


};
