
/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base							2016-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKBrdf_Stub_Base : public ISKBrdf_Stub
{
	private:
		skBRDF*								m_brdf;

	public:
											ISKBrdf_Stub_Base				( skBRDF* brdf);
		virtual 						   ~ISKBrdf_Stub_Base				() override;
		virtual nxUnknown*					RawObjectPointer				() override { return m_brdf;}
		virtual bool						BRDF							( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf);
		virtual bool						SetPropertyScalar				( const char* propertyname, double value ) override;							//!< [NOT Supported in Base class]
		virtual bool						SetPropertyArray				( const char* propertyname, const double* value, int numpoints ) override;	//!< [NOT Supported in Base class]
		virtual bool						SetPropertyObject				( const char* propertyname, nxUnknown* object ) override;					//!< [NOT Supported in Base class]
		virtual bool						SetPropertyString			    ( const char* propertyname, const char* str) override;
};


/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_LambertianAlbedo		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKBrdf_Stub_LambertianAlbedo : public ISKBrdf_Stub_Base
{
	private:
		SKTRAN_BRDF_Lambertian*				m_lambertian;

	public:
											ISKBrdf_Stub_LambertianAlbedo	( SKTRAN_BRDF_Lambertian* lambertian);
		virtual 						   ~ISKBrdf_Stub_LambertianAlbedo	() override;
		virtual bool						SetPropertyScalar				( const char* propertyname, double value) override;	//!< [NOT Supported in Base class]
};

/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_LambertianAlbedo		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKBrdf_Stub_Snow_Kokhanovsky2012 : public ISKBrdf_Stub_Base
{
	private:
		SKTRAN_BRDF_Snow_Kokhanovsky2012*	m_snowbrdf;

	public:
											ISKBrdf_Stub_Snow_Kokhanovsky2012	( SKTRAN_BRDF_Snow_Kokhanovsky2012* snow_brdf);
		virtual 						   ~ISKBrdf_Stub_Snow_Kokhanovsky2012	() override;
		virtual bool						SetPropertyArray					( const char* propertyname, const double* value, int numpoints) override;	//!< [NOT Supported in Base class]
};

/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Roujean		 2016- 12- 12*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKBrdf_Stub_Roujean : public ISKBrdf_Stub_Base
{
	private:
		SKTRAN_BRDF_Roujean*				m_roujeanbrdf;

	public:
											ISKBrdf_Stub_Roujean	( SKTRAN_BRDF_Roujean* roujeanbrdf);
		virtual 						   ~ISKBrdf_Stub_Roujean	() override;
		virtual bool						SetPropertyArray	    ( const char* propertyname, const double* value, int numpoints) override;	//!< [NOT Supported in Base class]
		virtual bool						SetPropertyScalar	    ( const char* propertyname, double value) override;	//!< [NOT Supported in Base class]
};

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Roujean_Kernel		            2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_Roujean_Kernel : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_Roujean_Kernel*			m_roujeankernel;

public:
										ISKBrdf_Stub_Roujean_Kernel(SKTRAN_BRDF_Roujean_Kernel* librdf);
	virtual 						   ~ISKBrdf_Stub_Roujean_Kernel() override;
};

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Li		            2017-03-06*/
/** **/
/*---------------------------------------------------------------------------*/

class ISKBrdf_Stub_Li_Kernel : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_Li_Kernel*	m_librdf;

public:
										ISKBrdf_Stub_Li_Kernel		(SKTRAN_BRDF_Li_Kernel* librdf);
	virtual 						   ~ISKBrdf_Stub_Li_Kernel		() override;
	virtual bool						SetPropertyArray	(const char* propertyname, const double* value, int numpoints) override;	//!< [NOT Supported in Base class]
};

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Ross		            2017-03-08*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_Ross_Kernel : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_Ross_Kernel*	m_rossbrdf;

public:
										ISKBrdf_Stub_Ross_Kernel		(SKTRAN_BRDF_Ross_Kernel* rossbrdf);
	virtual 						   ~ISKBrdf_Stub_Ross_Kernel	() override;
};

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Cox_Munk		            2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_Cox_Munk : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_CoxMunk*	m_waterbrdf;

public:
										ISKBrdf_Stub_Cox_Munk(SKTRAN_BRDF_CoxMunk* waterbrdf);
	virtual 						   ~ISKBrdf_Stub_Cox_Munk() override;
	virtual bool						SetPropertyArray(const char* propertyname, const double* value, int numpoints) override;	//!< [NOT Supported in Base class]
};

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Rahman		            2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_Rahman : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_Rahman*	m_rahmanbrdf;

public:
										ISKBrdf_Stub_Rahman(SKTRAN_BRDF_Rahman* rahmanbrdf);
	virtual 						   ~ISKBrdf_Stub_Rahman() override;
	virtual bool						SetPropertyArray(const char* propertyname, const double* value, int numpoints) override;	//!< [NOT Supported in Base class]
};

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_Hapke		            2017-07-28*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_Hapke : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_Hapke*	m_hapkebrdf;

public:
										ISKBrdf_Stub_Hapke(SKTRAN_BRDF_Hapke* hapkebrdf);
	virtual 						   ~ISKBrdf_Stub_Hapke() override;
	virtual bool						SetPropertyArray(const char* propertyname, const double* value, int numpoints) override;	//!< [NOT Supported in Base class]
};

/*-----------------------------------------------------------------------------
*		ISKBrdf_Stub_LinearCombination		  2017-07-31*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_LinearCombination : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_LinearCombination*		m_brdf;

public:
										ISKBrdf_Stub_LinearCombination(SKTRAN_BRDF_LinearCombination* modisbrdf);
	virtual 						   ~ISKBrdf_Stub_LinearCombination() override;
	virtual bool						SetPropertyArray(const char* propertyname, const double* value, int numpoints) override;	//!< [NOT Supported in Base class]
	virtual bool						SetPropertyObject(const char* propertyname, nxUnknown* object) override;
	virtual bool						SetPropertyScalar(const char* propertyname, double value) override;
};

/*-----------------------------------------------------------------------------
*		ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal		  2017-07-31*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal*	m_modisbrdf;

public:
										ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal(SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal* modisbrdf);
	virtual 						   ~ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal() override;
	virtual bool						SetPropertyArray(const char* propertyname, const double* value, int numpoints) override;	//!< [NOT Supported in Base class]
};

/*-----------------------------------------------------------------------------
*		ISKBrdf_Stub_UserDefinedLatLon		  2020-01-10*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_UserDefinedLatLon : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_UserDefinedLatLon*	m_brdf;

	void								MakeScalarSetFunctions();
	void								MakeVectorSetFunctions();
	void								MakeObjectSetFunctions();

	std::vector<double>					m_latitudes;
	std::vector<double>					m_longitudes;

	nx2dArray<skBRDF*>					m_brdfs;
	size_t								m_brdfindex;

public:
	ISKBrdf_Stub_UserDefinedLatLon(SKTRAN_BRDF_UserDefinedLatLon* brdf);
	virtual 						   ~ISKBrdf_Stub_UserDefinedLatLon() override;
	virtual bool						SetPropertyScalar(const char* propertyname, double value) override;
	virtual bool						SetPropertyArray(const char* propertyname, const double* value, int numpoints) override;
	virtual bool						SetPropertyObject(const char* propertyname, nxUnknown* object) override;

};

/*-----------------------------------------------------------------------------
*		ISKBrdf_Stub_UserDefinedLatLon		  2020-01-10*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_Plane : public ISKBrdf_Stub_Base
{
private:
	skBRDF_AlbedoPlane*	m_brdf;

	void								MakeObjectSetFunctions();

public:
    ISKBrdf_Stub_Plane(skBRDF_AlbedoPlane* brdf);
	virtual 						   ~ISKBrdf_Stub_Plane() override;
    virtual bool						SetPropertyObject(const char* propertyname, nxUnknown* object) override;
};


/*-----------------------------------------------------------------------------
*		ISKBrdf_Stub_UserDefinedLatLon		  2020-01-10*/
/** **/
/*---------------------------------------------------------------------------*/
class ISKBrdf_Stub_SpectralVarying : public ISKBrdf_Stub_Base
{
private:
	SKTRAN_BRDF_SpectralVarying*	    m_brdf;

	void								MakeScalarSetFunctions();
	void								MakeVectorSetFunctions();
	void								MakeObjectSetFunctions();

	std::vector<double>					m_wavelengths;

	nx1dArray<skBRDF*>					m_brdfs;
	size_t								m_brdfindex;

public:
	ISKBrdf_Stub_SpectralVarying(SKTRAN_BRDF_SpectralVarying* brdf);
	virtual 						   ~ISKBrdf_Stub_SpectralVarying() override;
	virtual bool						SetPropertyScalar(const char* propertyname, double value) override;
	virtual bool						SetPropertyArray(const char* propertyname, const double* value, int numpoints) override;
	virtual bool						SetPropertyObject(const char* propertyname, nxUnknown* object) override;

};