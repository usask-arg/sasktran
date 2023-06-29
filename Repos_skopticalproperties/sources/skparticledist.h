

/*---------------------------------------------------------------------------
 *					class skRTParticleDist				2003-10-8			*/
/**	\ingroup skpartdist
 *	A base class for calculating parameterized particle distributions. This
 *	is typically used in aerosol calculations. Most end users do not normally
 *	need access to particle distributions as developers will normally use the
 *	correct distribution in the optical properties of the aerosol they are
 *	developing.*/
/*-------------------------------------------------------------------------*/

class skRTParticleDist : public nxUnknown
{
	public:
		enum ENUM_PDIST		{PDIST_2GAMMA, PDIST_3GAMMA, PDIST_BIMODAL, PDIST_LOGNORMAL, PDIST_POWERLAW, PDIST_MONODISPERSE };

	public:
		bool									IsSameDistributionAs				( const skRTParticleDist* other ) const;

	public:
												skRTParticleDist					() {} // m_minradius = 0.001; m_maxradius = 100.0;}
		virtual								   ~skRTParticleDist					()  {}
		virtual bool							CreateClone							(skRTParticleDist** userclone) const = 0;
		virtual bool							GetQuadratureRadii					( double* minradius, double* maxradius) const = 0;
		virtual double							Distribution						( double R )                    = 0;
		virtual bool							SetDistributionParameters			( double A, double B, double C) = 0;
		virtual void							GetDistributionParameters			( double* A, double* B, double* C) = 0;
		virtual size_t							NumDistributionParameters			()  const = 0;
		virtual bool							GetDistributionParameterArray		( double* parameters, size_t maxparams, size_t* numparams ) const = 0;
		virtual bool							GetDistributionParameterSpeciesID	( CLIMATOLOGY_HANDLE* paramids, size_t maxparams, size_t* numparams) const = 0;
		virtual skRTParticleDist::ENUM_PDIST	Type								() const                           = 0;
		virtual nxString						CachingDescriptor					() const { NXASSERT((false)); return nxString ("Undefined_CachingDescriptor");	}
		virtual double							ASA_To_N							(double asa_um2percm3) const;
		virtual double							N_To_ASA							(double N_percm3) const;
};


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal						2003-10-9*/
/** \ingroup skpartdist
 *	This class implements the log normal particle distribution. See Eqn 2.60 in
 *	 Hansen And Travis. Space Sci. Rev., 16. 527-610, 1974. **/
/*---------------------------------------------------------------------------*/

class skRTParticleDist_LogNormal : public skRTParticleDist
{
	private:
		double								m_QTEMP;
		double								m_A2;
		double								m_A3;
		
	private:
		bool								Copy								( const skRTParticleDist_LogNormal& other);

	public:
		double								ModeRadiusMicrons() const;			//!< Returns the  mode radius in microns
		double								ModeWidth() const;					//!< Returns the mode with (typically around 1.6)
		double								Reff() const;						//!< Returns the Effective Radius in microns 
		double								Veff() const;						//!< Returns the effective variance (dimensionless)

	public:
											skRTParticleDist_LogNormal			(){m_QTEMP  = 0.0; m_A2 = 0.0; m_A3 = 0.0;}				
	virtual								   ~skRTParticleDist_LogNormal			(){};
	virtual bool							CreateClone							(skRTParticleDist** userclone) const override;
	virtual double							Distribution						( double R ) override;
	virtual	bool							SetDistributionParameters			( double  RG,  double  MODEWIDTH,   double  UNUSED) override;
	virtual void							GetDistributionParameters			( double* RG, double*  MODEWIDTH,   double* ALWAYSZERO) override;
	virtual size_t							NumDistributionParameters			()  const override { return 2;}
	virtual bool							GetDistributionParameterArray		( double* parameters, size_t maxparams, size_t* numparams ) const override;
	virtual bool							GetDistributionParameterSpeciesID	( CLIMATOLOGY_HANDLE* paramids, size_t maxparams, size_t* numparams) const override;
	virtual skRTParticleDist::ENUM_PDIST	Type								() const  override { return PDIST_LOGNORMAL;}
	virtual nxString						CachingDescriptor					() const override;							
	virtual bool							GetQuadratureRadii					(double* minradius, double* maxradius) const override;
	virtual double							ASA_To_N							(double asa_um2percm3) const override;
	virtual double							N_To_ASA							(double N_percm3) const override;
};


/*-----------------------------------------------------------------------------
 *					class skRTParticleDist_PowerLaw					2003-10-9*/
/**	\ingroup skpartdist
 *	This class implements the log normal particle distribution. See Eqn 2.61 in
 *	Hansen And Jensen. Space Sci. Rev., 16. 527-610, 1974. 
 *  N(R)=CONST*R**(-A) FROM RMIN=B TO RMAX=C AND ZERO OTHERWISE */
 /*-------------------------------------------------------------------------*/

class skRTParticleDist_PowerLaw : public skRTParticleDist
{
	private:
		double								m_A;
		double								m_B;
		double								m_C;
		double								m_CONST;

	private:
		bool									Copy								( const skRTParticleDist_PowerLaw& other);

	public:
												skRTParticleDist_PowerLaw			(){m_A = 0.0; m_B = 0.0; m_C = 0.0; m_CONST = 0.0;}				
		virtual								   ~skRTParticleDist_PowerLaw			(){};
		virtual bool							CreateClone							(skRTParticleDist** userclone) const;
		virtual double							Distribution						( double R );
		virtual	bool							SetDistributionParameters			( double A, double B, double C);
		virtual void							GetDistributionParameters			( double* A, double* B, double* C);
		virtual size_t							NumDistributionParameters			()  const { return 3;}
		virtual bool							GetDistributionParameterArray		( double* parameters, size_t maxparams, size_t* numparams ) const;
		virtual bool							GetDistributionParameterSpeciesID	( CLIMATOLOGY_HANDLE* paramids, size_t maxparams, size_t* numparams) const;
		virtual skRTParticleDist::ENUM_PDIST	Type								() const { return PDIST_POWERLAW;}
		virtual nxString						CachingDescriptor					() const;
		virtual bool							GetQuadratureRadii					(double* minradius, double* maxradius) const;
};

/*---------------------------------------------------------------------------
 *					class skRTParticleDist_2Gamma					2003-10-8*/
/**	\ingroup skpartdist
 *	This class implements the 2 parameter gamma distribution. See Eqn 2.56 in
 *	Hansen And Jensen. Space Sci. Rev., 16. 527-610, 1974. */
/*-------------------------------------------------------------------------*/

class skRTParticleDist_2Gamma : public skRTParticleDist
{
	private:
		double									m_A;
		double									m_B;
		double									m_AB;
		double									m_RBM3;
		double									m_NORM;

	private:
		bool									Copy								( const skRTParticleDist_2Gamma& other);

	public:
												skRTParticleDist_2Gamma				();
		virtual								   ~skRTParticleDist_2Gamma				(){};
		virtual bool							CreateClone							(skRTParticleDist** userclone) const;
		virtual double							Distribution						( double R );
		virtual	bool							SetDistributionParameters			( double A, double B, double);
		virtual void							GetDistributionParameters			( double* A, double* B, double* C);
		virtual size_t							NumDistributionParameters			()  const { return 2;}
		virtual bool							GetDistributionParameterArray		( double* parameters, size_t maxparams, size_t* numparams ) const;
		virtual bool							GetDistributionParameterSpeciesID	( CLIMATOLOGY_HANDLE* paramids, size_t maxparams, size_t* numparams) const;
		virtual skRTParticleDist::ENUM_PDIST	Type								() const { return PDIST_2GAMMA;}
		virtual nxString						CachingDescriptor					() const;
		virtual bool							GetQuadratureRadii					(double* minradius, double* maxradius) const;
};

/*---------------------------------------------------------------------------
 *					class skRTParticleDist_3Gamma					2003-10-8*/
/*	\ingroup skpartdist
 *	This class implements the 3 parameter gamma distribution.
 *
 *	THIS IS THE THREE PARAMETER GAMMA DISTRIBUTION (MODIFIED) DEIRMENDJIAN.
 *	THIS TYPE OF SIZE DISTRIBUTION IS USED FOR ALL OF THE STANDARD WATER CLOUD
 *	AND HAZE DISTRIBUTIONS IN THE INPUT SUBROUTINE. */
/*-------------------------------------------------------------------------*/

class skRTParticleDist_3Gamma : public skRTParticleDist
{
	private:
		double	m_A;
		double	m_B;
		double  m_C;
		double	m_DLB;

	private:
		bool									Copy								( const skRTParticleDist_3Gamma& other);

	public:
												skRTParticleDist_3Gamma				();
		virtual								   ~skRTParticleDist_3Gamma				(){};
		virtual bool							CreateClone							(skRTParticleDist** userclone) const;
		virtual double							Distribution						( double R );
		virtual bool							SetDistributionParameters			( double A, double B, double C);
		virtual void							GetDistributionParameters			( double* A, double* B, double* C);
		virtual bool							GetDistributionParameterSpeciesID	( CLIMATOLOGY_HANDLE* paramids, size_t maxparams, size_t* numparams) const;
		virtual size_t							NumDistributionParameters			()  const { return 3;}
		virtual bool							GetDistributionParameterArray		( double* parameters, size_t maxparams, size_t* numparams ) const;
		virtual skRTParticleDist::ENUM_PDIST	Type								() const { return PDIST_3GAMMA;}
		virtual nxString						CachingDescriptor					() const;
		virtual bool							GetQuadratureRadii					(double* minradius, double* maxradius) const;
};


/*---------------------------------------------------------------------------
 *					class skRTParticleDist_BimodalGamma				2003-10-9*/
/**	\ingroup skpartdist
 *	This class implements the bimodal gamma distribution. See Eqn 2.59 in
 *	Hansen And Jensen. Space Sci. Rev., 16. 527-610, 1974. */
/*-------------------------------------------------------------------------*/

class skRTParticleDist_BimodalGamma : public skRTParticleDist
{
	private:
		double	m_A1B;
		double	m_A2B;
		double	m_DLA1B;
		double	m_DLA2B;
		double	m_RBM3;

	private:
		bool									Copy								( const skRTParticleDist_BimodalGamma& other);
	public:
												skRTParticleDist_BimodalGamma		(){m_A1B=0.0; m_A2B=0.0; m_DLA1B=0.0; m_RBM3=0.0;}				
		virtual								   ~skRTParticleDist_BimodalGamma		(){};
		virtual bool							CreateClone							(skRTParticleDist** userclone) const;
		virtual double							Distribution						( double R );
		virtual	bool							SetDistributionParameters			( double A, double B, double C);
		virtual void							GetDistributionParameters			( double* A, double* B, double* C);
		virtual size_t							NumDistributionParameters			()  const { return 5;}
		virtual bool							GetDistributionParameterArray		( double* parameters, size_t maxparams, size_t* numparams ) const;
		virtual bool							GetDistributionParameterSpeciesID	( CLIMATOLOGY_HANDLE* paramids, size_t maxparams, size_t* numparams) const;
		virtual skRTParticleDist::ENUM_PDIST	Type								() const { return PDIST_BIMODAL;}
		virtual nxString						CachingDescriptor					() const;
		virtual bool							GetQuadratureRadii					(double* minradius, double* maxradius) const;
};
