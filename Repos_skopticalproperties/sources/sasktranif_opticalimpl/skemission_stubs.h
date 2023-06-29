
/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_BaseEngine		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKEmission_Stub_Base : public ISKEmission_Stub
{
	private:
		skEmission*							m_emission;

	protected:
	public:
											ISKEmission_Stub_Base			( skEmission* optprop);
		virtual 						   ~ISKEmission_Stub_Base			() override;
		virtual nxUnknown*					RawObjectPointer				() override { return m_emission;}
		virtual bool						UpdateLocation					( const GEODETIC_INSTANT& pt, bool isground ) override;
		virtual bool						UpdateCache						( const GEODETIC_INSTANT& pt ) override;
		virtual bool						IsotropicEmission				( double wavenumber, double* isotropicradiance) override;
		virtual bool						IsotropicEmissionArray			( const double * wavenumber, int numwavenumber, double * isotropicradiance, int numisotrop) override;
//		virtual bool						SetPropertyScalar				( const char* propertyname, double value ) override;							//!< [NOT Supported in Base class]
//		virtual bool						SetPropertyArray				( const char* propertyname, const double* value, int numpoints ) override;	//!< [NOT Supported in Base class]
//		virtual bool						SetPropertyObject				( const char* propertyname, nxUnknown* object ) override;					//!< [NOT Supported in Base class]
//		virtual bool						SetPropertyString			    ( const char* propertyname, const char* str) override;

};




/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Tabulated_HeightWavelength		 2015- 3- 11*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKEmission_Stub_Tabulated_HeightWavelength : public ISKEmission_Stub_Base
{
	private:
		skEmission_Tabulated_HeightWavelength*		m_userdefinedemission;
		nx1dArray<double>							m_currentheightarray;
		nx1dArray<double>							m_currentwavelengtharray;

	private:
		void								MakeSetPropertyFunctions();

	public:
											ISKEmission_Stub_Tabulated_HeightWavelength	( skEmission_Tabulated_HeightWavelength* climate);
		virtual 						   ~ISKEmission_Stub_Tabulated_HeightWavelength	() override;
};


/*-----------------------------------------------------------------------------
 *					ISKEmission_Stub_Thermal                      2017- 8- 30*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKEmission_Stub_Thermal : public ISKEmission_Stub_Base
{
	private:
		skEmission_Thermal*					m_thermalemission;

	private:
		void								MakeSetPropertyFunctions();

	public:
											ISKEmission_Stub_Thermal (skEmission_Thermal* emission);
										   ~ISKEmission_Stub_Thermal () override;
};



/*---------------------------------------------------------------------------
 *             Class ISKEmission_Stub_HitranChemical              2020-08-21 */
/** **/
/*---------------------------------------------------------------------------*/

class ISKEmission_Stub_HitranChemical : public ISKEmission_Stub_Base
{
	private:
		skEmission_HitranChemical*			m_thermalemission;

	private:
		void								MakeSetPropertyFunctions();

	public:
											ISKEmission_Stub_HitranChemical (skEmission_HitranChemical* emission);
										   ~ISKEmission_Stub_HitranChemical () override;
};

