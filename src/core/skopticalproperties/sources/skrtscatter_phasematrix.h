

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ParticleBase		2013-1-25*/
/** \ingroup aerosolskopticalprop 
 *	A base class for particles that scatter using MIE or other theories.
 *	Most of these calculatiosn require supply a particle distribution
 *	eg log normal, gamma and a refractive index. This paradigm works for both
 *	spherical MIE and non-spherical TMatrix codes.
 */
/*---------------------------------------------------------------------------*/

class skOpticalProperties_ParticleBase : public skOpticalProperties
{
	public:
		virtual							   ~skOpticalProperties_ParticleBase	(){}
		virtual bool						Set_RefractiveIndex					( skRTRefractiveIndex* newri )           = 0;
		virtual bool						Set_ParticleDistribution			( skRTParticleDist*	   newdistribution ) = 0;
		virtual skRTParticleDist*			Get_Distribution					()										 = 0;
		virtual skRTRefractiveIndex*		Get_RefractiveIndex					()										 = 0;
};

/*---------------------------------------------------------------------------
 *					class skOpticalProperties_MieAerosol					2003-10-9 */
/**	\ingroup aerosolskopticalprop
 *	This is a class that calculates the scattering phase matrix for a
 *	particle distribution of spherical particles obeying mie scattering. Most
 *	users will use class #skOpticalProperties_MieAerosolCached
 *
 *	Internally the class stores a table giving the phase matrix at a set
 *	of interpolation cos( scattering angles).   The table is calculated
 *	when the base class sk_ExtinctionPhaseMatrix calls Interpolate_PhaseMatrixTables
 *  The call is not made until the user wants a Phase Matrix or cross-section
 *	value;
 *
 *	Note that the class does not detect cases where the user changes the
 *	underlying paramaters of the associated refractive index or particle
 *	distribution.  It is not anticpated that the user will "mess" with these
 *	parameters between usage.
 *
 *	The class implicitly assumes that the particle distribution is normalized
 *	to unity.
 *
 *	Calculating the lookup table is quite intensive as for each wavelength
 *	it requires an evaluation of the mie scattering code over all the
 *	possible angles and particle radii.
 */
/*--------------------------------------------------------------------------*/

class skOpticalProperties_MieAerosol : public skOpticalProperties_ParticleBase
{
	private:
		std::mutex							m_mutex_xsectionslock;
		bool									m_isdirty;						//!< Flag to check when the optical property needs to be updated
		double									m_internalwavenumber;			//!< Wavenumber for the current internal calculated optical properties
		sk_MieSphericalWiscombeWrapper			m_mie;							//!< The Mie object for calculating phase matrix of a single sphere.
		skRTParticleDist*						m_distribution;					//!< The object describing the particle distribution
		skRTRefractiveIndex*					m_refractiveindex;				//!< Calculates the refractive index at a given wavenumber
		nxGaussQuadrature<>						m_quadrature;					//!< The quadrature object for integrating over radius distribution
		nx2dArray<double>						m_Pij;							//!< double[ 4, num_Cosangles], The internal Pij phase matrix lookup tables calculated by m_mie [P11, P21, P33, P34]
		size_t									m_numlegendre;
		nx2dArray<double>						m_legendre;

	private:
		void								ReleaseRI						();
		void								init							();
		void								ReleaseDistribution				();
		bool								InterpolatePhaseMatrixTables	( double mu, skRTPhaseMatrix* P);
		bool								CalculateCrossSectionsInternal	( double wavenumber, double* absxs, double* extxs, double* scattxs);
//		bool								DeepCopy						( const skOpticalProperties_MieAerosol& other );
											skOpticalProperties_MieAerosol	( const skOpticalProperties_MieAerosol& other );	// dont allow copy constructor
		skOpticalProperties_MieAerosol&		operator =						( const skOpticalProperties_MieAerosol& other );	// Dont allow assignment operator

	public:
											skOpticalProperties_MieAerosol	();
		virtual							   ~skOpticalProperties_MieAerosol	();
		void								SetDirty						() { m_isdirty = true;}
		void								Set_DistributionQuadratureOrder	( int order )	{ SetDirty(); return m_quadrature.SetOrder(order);}
		void								LoadDefaultStratosphere			();
		void								LoadDefaultTroposphere			();
		bool								Set_SymmetricScatterAngles		( double *startangle,	double *endangle, double *resolution, int numsteps) { SetDirty(); return m_mie.Set_SymmetricScatterAngles( startangle, endangle, resolution, numsteps );}
		bool								Set_NumLegendreMoments(size_t nummoments) { m_numlegendre = nummoments; return m_mie.Set_MaxLegendreMoment((int)m_numlegendre); }
		size_t								Get_MaxLegendreMoments() { return m_numlegendre; }

	public:
		virtual bool						Set_RefractiveIndex				( skRTRefractiveIndex* newri );
		virtual bool						Set_ParticleDistribution		( skRTParticleDist*	   newdistribution );
		virtual skRTParticleDist*			Get_Distribution				()				{ return m_distribution;}
		virtual skRTRefractiveIndex*		Get_RefractiveIndex				()				{ return m_refractiveindex;}
	
	public:
		virtual bool						SetAtmosphericState				( skClimatology* neutralatmosphere)  override;							// {return true;}			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						SetLocation						(const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )  override;	// {return true;}			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						InternalClimatology_UpdateCache		( const GEODETIC_INSTANT& /*pt*/)  override { return true;}
		virtual bool						CalculateCrossSections			( double wavenumber, double* absxs, double* extxs, double* scattxs)  override;																	//!< Calculate cross-sections at the specified wave-number.
		virtual bool						IsScatterer						() const    override { return true;}																//!< Returns true if this particles scatters radiation
		virtual bool						IsAbsorber						() const    override { return true;}                                                      		//!< Returns true if this particles absorbs radiation radiation
		virtual bool						CalculatePhaseMatrix			( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix )  override;			//!< Calculate the phase matrix at the specified wavenumber and scattering angle
		virtual bool						LegendreCoefficientsP11         ( double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff ) override;

};

/*-----------------------------------------------------------------------------
 *					class skOpticalProperties_MieAerosolCached		2008-4-16*/
/** \ingroup aerosolskopticalprop
 *	A class very similar to skOpticalProperties_MieAerosol except that it
 *	caches Mie scattering calculations in files on the users system.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_MieAerosolCached : public skOpticalProperties_MieAerosol
{
	private:
		struct ThreadData
		{
			public:
//				size_t							m_threadindex;
				bool							m_isdirtycached;				//!< have requested cache properties changed?
				double							m_wavenumber;					//!< The wavenumber of this cache
				double							m_absxs;						//!< The absorption cross-section  
				double							m_extxs;						//!< The extinction cross-section
				double							m_scattxs;						//!< The scattering cross-section
				std::vector<skRTPhaseMatrix>	m_phasematrix;					//!< An array skRTPhaseMatrix[numscatteringangles]
				std::vector<double>				m_legendrep11;					//!< An array of legendre moments for P11

			public:
												ThreadData();
				size_t							NumScatteringAngles() const						{ return m_phasematrix.size();}
				bool							SetNumScatteringAngles( size_t numangles)		{ m_phasematrix.resize(numangles); return true;}
				size_t							NumLegendreMoments() const					    { return m_legendrep11.size(); }
				bool							SetNumLegendreMoments(size_t nummoments)		{ m_legendrep11.resize(nummoments); return true; }
		};

	private:
				std::map<size_t, ThreadData>				m_data;							//!< Thread storage. The integer key is the value returned by omp_get_thread_num
		typedef std::map<size_t, ThreadData>::iterator		iterator;
		nxString										m_cachebasedir;					//!< the base directory for the cache files, loaded once in the constructor
		size_t											m_numangles;					//!< The number of angles used in each threads phasematrix;
		size_t											m_nummoments;

	private:
		size_t									NumDataThreads						() const	{ return m_data.size();}
		bool									LookupUpThreadData					(ThreadData** data );
		bool									CheckDirtyAndUpdate					(double wavenumber, ThreadData* data);
		void									SetDirtyCached						();
		bool									CreateTables						(ThreadData* data);
		nxString								FullCacheName						(ThreadData* data);
		bool									WriteCacheFile						( const char* filename, ThreadData* data );
		bool									ReadCacheFile						( const char* filename, ThreadData* data );
		bool									UpdateTables						(ThreadData* data);
		nxString								LoadDirectoryNameFromRegistry		();
		bool									GetPhaseMatrix						( ThreadData* data, double cosscatterangle, skRTPhaseMatrix* phasematrix );
												skOpticalProperties_MieAerosolCached ( const skOpticalProperties_MieAerosolCached& other );	// Dont allow copy constructor
		skOpticalProperties_MieAerosolCached&	operator =					 ( const skOpticalProperties_MieAerosolCached& other );	// Dotn allow assignment operator

	public:
												skOpticalProperties_MieAerosolCached		();
		virtual								   ~skOpticalProperties_MieAerosolCached		() override;
		bool									Set_RefractiveIndex					( skRTRefractiveIndex* ri );
		bool									Set_ParticleDistribution			( skRTParticleDist* distribution );
		bool									SetNumScatteringAngles				( size_t numangles );


		virtual bool							CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs)  override;																	//!< Calculate cross-sections at the specified wave-number.
		virtual bool							CalculatePhaseMatrix				( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix) override;			//!< Calculate the phase matrix at the specified wavenumber and scattering angle
		virtual bool							LegendreCoefficientsP11             ( double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) override;
};

//
/*---------------------------------------------------------------------------
 *					class skOpticalProperties_IceCrystal				2009-6-3 */
/**	\ingroup aerosolskopticalprop
 *	[THREAD-SAFE]This class calculates the scattering phase matrix and cross sections for a
 *	specified distribution of nonspherical particles.
 *	Internally the class stores a table giving the phase matrix at a set
 *	of interpolation cos( scattering angles).   The table is calculated
 *	when the base class sk_ExtinctionPhaseMatrix calls Interpolate_PhaseMatrixTables
 *  The call is not made until the user wants a Phase Matrix or cross-section
 *	value;
 *	Note that the class does not detect cases where the user changes the
 *	underlying paramaters of the associated refractive index or particle
 *	distribution.  It is not anticpated that the user will "mess" with these
 *	parameters between usage. 
 *
 *	This class provides for using a number of numerical algorithms for 
 *	computing the scattering properties: randomly- or oriented rotationally-
 *	symmetric particles with the T-matrix code, or for arbitrary target geometries
 *	with the Discrete Dipole approximation code. 
 *	NOTE: for the randomly-oriented T-matrix algorithm, computation of the particle 
 *	size distribution is performed within the DLL function 'TMATRIXRANDOM', so the 
 *	size distribution must be selected before computing any particle
 *	scattering properties.
 *
 *	Calculating the lookup table is quite intensive as for each wavelength
 *	it requires an evaluation of the T-matrix/DDA scattering code over all the
 *	possible angles and particle radii. Computation time increases greatly with
 *	effective size parameter and particle aspect ratio.
 * 
 *	The scattering algorithm object must be the first attribute added to the 
 *	scattering particle.  For example, if 'optice' is the scattering particle object
 *	with the associated data types:
 *		skOpticalProperties_IceCrystalCached*	optice;
 *		skRTRefractiveIndex*				riIce;
 *		sk_NonsphericalParticle*			scatIce;
 *		skRTParticleDist*					psdIce;
 *	and they are instantiated as:
 *		optice	= new skOpticalProperties_IceCrystalCached;
 *		riIce	= new skRTRefractiveIndex_ICE;
 *		scatIce	= new sk_TMatrixRandomWrapper;
 *		psdIce	= new skRTParticleDist_2Gamma;
 *		psdIce->SetParameters( radius,width,0 );
 *	then the objects should be attached in the following order:
 *		optice->Set_ScatterAlgorithm( scatIce );
 * 		optice->Set_ParticleDistribution( psdIce );
 * 		optice->Set_RefractiveIndex( riIce );
/*--------------------------------------------------------------------------*/

class skOpticalProperties_IceCrystal : public skOpticalProperties_ParticleBase
{
	private:
		sk_NonsphericalParticle*						m_nonsphere;					//!< The object for calculating phase matrix 
		skRTParticleDist*								m_distribution;					//!< The object describing the particle distribution
		skRTRefractiveIndex*							m_refractiveindex;				//!< Calculates the refractive index at a given wavenumber
		nxGaussQuadrature<>								m_quadrature;					//!< The quadrature object for integrating over radius distribution
		        std::map< size_t, std::vector<skRTPhaseMatrix> >		m_PijArray;						//!< [NumThreads] [num_Cosangles], The Pij phase matrix lookup tables calculated by scattering algorithmn [P11,P12,P13,P14,P21,P22,P23,P24,P31,P32,P33,P34,P41,P42,P43,P44]
		typedef std::map< size_t, std::vector<skRTPhaseMatrix> >::iterator		iterator;						//!< [NumThreads] [num_Cosangles], The Pij phase matrix lookup tables calculated by scattering algorithmn [P11,P12,P13,P14,P21,P22,P23,P24,P31,P32,P33,P34,P41,P42,P43,P44]

	private:
		void								init									();
		bool								LookupUpThreadData						(std::vector<skRTPhaseMatrix>** data );
		void								ReleaseRI								();
		void								ReleaseDistribution						();
		void								ReleaseScatterAlgorithm					();
		bool								InterpolatePhaseMatrixTables			( double mu, skRTPhaseMatrix* P, std::vector<skRTPhaseMatrix>* Pij ) const;
		bool								IntegrateOverSingleRadiusCode			( double wavenum, double* absxs, double* extxs, double* scattxs, std::vector<skRTPhaseMatrix>* Pij  );
		bool								IntegrateOverRandomTmatrixDistribution	( double wavenum, double* absxs, double* extxs, double* scattxs, std::vector<skRTPhaseMatrix>* Pij  );

		bool								CalculateCrossSectionsInternal			( double wavenumber, double* absxs, double* extxs, double* scattxs, std::vector<skRTPhaseMatrix>* Pij );
											skOpticalProperties_IceCrystal			( const skOpticalProperties_IceCrystal& other );	// Dont allow copy constructor
		skOpticalProperties_IceCrystal&		operator =								( const skOpticalProperties_IceCrystal& other );	// Dont allow assignment operator

	public:
											skOpticalProperties_IceCrystal	();
		virtual							   ~skOpticalProperties_IceCrystal	();
		bool								Set_AspectRatio					( double aspectRatio );
		bool								Set_ScatterAlgorithm			( sk_NonsphericalParticle* newscat );
		bool								Set_SymmetricScatterAngles		( double *startangle,	double *endangle, double *resolution, int numsteps) { return m_nonsphere->Set_SymmetricScatterAngles( startangle, endangle, resolution, numsteps );}
		bool								Set_AnyScatterAngles			( nx1dArray<double>& degrees)		{ return m_nonsphere->Set_AnyScatterAngles( degrees );}
		void								Set_DistributionQuadratureOrder	( int order )						{ return m_quadrature.SetOrder(order);}
		sk_NonsphericalParticle*			Get_ScatterAlgorithm			()									{ return m_nonsphere; }
		nxString							Get_DescriptiveParticleString	()									{ return m_nonsphere->Get_DescriptiveParticleString(); }

	public:
		virtual bool						Set_RefractiveIndex				( skRTRefractiveIndex* newri );
		virtual bool						Set_ParticleDistribution		( skRTParticleDist*	   newdistribution );
		virtual skRTParticleDist*			Get_Distribution				()									{ return m_distribution; }
		virtual skRTRefractiveIndex*		Get_RefractiveIndex				()									{ return m_refractiveindex; }

	public:
		virtual bool						SetAtmosphericState			( skClimatology* neutralatmosphere ) override; // {return true;}			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						SetLocation					( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  ) override; // {return true;}			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						InternalClimatology_UpdateCache	( const GEODETIC_INSTANT& /*pt*/) override			{ return true;}
		virtual bool						CalculateCrossSections		( double wavenumber, double* absxs, double* extxs, double* scattxs)  override;																	//!< Calculate cross-sections at the specified wave-number.
		virtual bool						IsScatterer					() const  override								{ return true;}																//!< Returns true if this particles scatters radiation
		virtual bool						IsAbsorber					() const  override								{ return true;}                                                      		//!< Returns true if this particles absorbs radiation radiation
		virtual bool						CalculatePhaseMatrix        ( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix) override;			//!< Calculate the phase matrix at the specified wavenumber and scattering angle
};

/*-----------------------------------------------------------------------------
 *					class skOpticalProperties_IceCrystalCached			2009-6-3
/**	\ingroup aerosolskopticalprop
 *	A class very similar to skOpticalProperties_IceCrystal except that it
 *	caches scattering calculations in files on the users system.
*/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_IceCrystalCached : public skOpticalProperties_IceCrystal
{
	private:
		struct ThreadData
		{
			public:
				bool							m_isdirtycached;				//!< have requested cache properties changed?
				double							m_wavenumber;					//!< The wavenumber of this cache
				double							m_absxs;						//!< The absorption cross-section  
				double							m_extxs;						//!< The extinction cross-section
				double							m_scattxs;						//!< The scattering cross-section
				std::vector<skRTPhaseMatrix>	m_phasematrix;					//!< An array skRTPhaseMatrix[numscatteringangles]

			public:
												ThreadData();
				size_t							NumScatteringAngles() const						{ return m_phasematrix.size();}
				bool							SetNumScatteringAngles ( size_t numangles)		{ m_phasematrix.resize(numangles); return true;}
		};

	private:
			    std::map<size_t, ThreadData>				m_data;							// The thread specific data
		typedef std::map<size_t, ThreadData>::iterator		iterator;
		nxString											m_cachebasedir;					//!< the base directory for the cache files, loaded once in the constructor
		size_t												m_numangles;

	private:
		size_t									NumDataThreads							() const	{ return m_data.size();}
		bool									LookupUpThreadData						(ThreadData** data );
		bool									CheckDirtyAndUpdate						(double wavenumber,  ThreadData* data);
		bool									CreateTables							( ThreadData* data);
		nxString								FullCacheName							( ThreadData* data);
		bool									WriteCacheFile							( const char* filename,  ThreadData* data ) const;
		bool									ReadCacheFile							( const char* filename,  ThreadData* data );
		bool									UpdateTables							( ThreadData* data);
		nxString								LoadDirectoryNameFromRegistry			();
		bool									GetPhaseMatrix							( ThreadData* data, double cosscatterangle, skRTPhaseMatrix* phasematrix );

												skOpticalProperties_IceCrystalCached	( const skOpticalProperties_IceCrystalCached& other );	// dont allow copy constructor
		skOpticalProperties_IceCrystalCached&	operator =								( const skOpticalProperties_IceCrystalCached& other );	// Dotn allow assignment operator

	public:
												skOpticalProperties_IceCrystalCached	();
		virtual								   ~skOpticalProperties_IceCrystalCached	() override;
		void									SetDirtyCached							();
		bool									Set_RefractiveIndex						( skRTRefractiveIndex* ri );
		bool									Set_ParticleDistribution				( skRTParticleDist* newdistribution );
		bool									Set_ScatterAlgorithm					( sk_NonsphericalParticle* newscat );
		bool									SetNumScatteringAngles					( size_t numangles );

		virtual bool							CalculateCrossSections					( double wavenumber, double* absxs, double* extxs, double* scattxs)  override;																	//!< Calculate cross-sections at the specified wave-number.
		virtual bool							CalculatePhaseMatrix					( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix) override;			//!< Calculate the phase matrix at the specified wavenumber and scattering angle
};


/*-----------------------------------------------------------------------------
 *					class skOpticalProperties_AerosolProfile				2008-2-28*/
/** \ingroup aerosolskopticalprop
 *	A base class that allows the optical properties of atmospheric particles
 *	to be calculated. The class combines the MIE/TMATRIX codes for individual
 *	particles with the refractive index and suitable particle distributions
 *	such as log normal and gamma. The class also provides the option for
 *	a climatology that provides the spatial variation of distribution
 *	parameters with height, time and location. (eg variations in mode radius
 *	and mode width can be specified via the climatology).
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_AerosolProfile : public skOpticalProperties
{
	private:
		skOpticalProperties_ParticleBase*	m_mieaerosol;					//!< Object calculates scattering properties and stores them online
		skClimatology*						m_moderadiusclimatology;		//!< Object that supports climatology (typically height profile) of mode radius, mode width and otehr parameters.
		skRTParticleDist*					m_distribution;					//!< Height profile of log normal distributions
		skRTRefractiveIndex*				m_refractiveindex;				//!< The chemical species of the aerosol material
		GEODETIC_INSTANT					m_lastpt;
		bool								m_isdirty;
		bool								m_heightdependent;

	private:
		void								init							();
		bool								CheckDirtyAndUpdate				();
		void								ReleaseResources				();
		void								SetDirty								() { m_isdirty = true;}
//		bool								DeepCopy								( const skOpticalProperties_AerosolProfile& other);
											skOpticalProperties_AerosolProfile		( const skOpticalProperties_AerosolProfile& other );		// dont allow copy constructor
		skOpticalProperties_AerosolProfile& operator = 								( const skOpticalProperties_AerosolProfile& other );		// Dont allow assignment operator

	public:
											skOpticalProperties_AerosolProfile		();
		virtual							   ~skOpticalProperties_AerosolProfile		();
		bool								SetOpticalProperties					( skOpticalProperties_ParticleBase* mieaerosol);
		bool								SetRefractiveIndex						( skRTRefractiveIndex*	ri			 );
		bool								SetParticleDistribution					( skRTParticleDist*		distribution );
		bool								SetParticleSizeClimatology				( skClimatology*		lognormalclimatology );

		bool								SetLogNormalProfileClimatology			( const double* altmeters, const double* moderadius_microns, const double* modewidth, size_t numalt);
		bool								SetLogNormalProfileClimatologyFromFile	( const char* filename  );
		bool								SetGammaProfileClimatology				( const double* altmeters, const double* effectiveradius_microns, const double* rate, size_t numalt);
		bool								GetDistributionParameter				( const CLIMATOLOGY_HANDLE& species, const GEODETIC_INSTANT& pt, double* value);
		bool								ASA_To_ExtinctionPerCm					( double wavenumber, double asa_um2percm3,    double* extinction_percm);
		bool								ExtinctionPerCm_to_ASA					( double wavenumber, double extinctionpercm,  double* asa_um2percm3);

	public:
		virtual bool						SetAtmosphericState					( skClimatology* neutralatmosphere ) override;
		virtual bool						SetLocation							( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  ) override;
		virtual bool						InternalClimatology_UpdateCache		( const GEODETIC_INSTANT& pt) override;
		virtual bool						CalculateCrossSections			( double wavenumber, double* absxs, double* extxs, double* scattxs) override;																	//!< Calculate cross-sections at the specified wave-number.
		virtual bool						CalculatePhaseMatrix			( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix) override;			//!< Calculate the phase matrix at the specified wavenumber and scattering angle
		virtual bool						IsScatterer						() const  override { return true;}																//!< Returns true if this particles scatters radiation
		virtual bool						IsAbsorber						() const  override { return true;}                                                      		//!< Returns true if this particles absorbs radiation radiation
		virtual bool						LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) override;
		virtual bool						SetHeightDependent(bool dependent) { m_heightdependent = dependent; return true; }
		virtual bool						IsHeightDependent() const override { return m_heightdependent; }
};


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileH2SO4		2008-2-29*/
/** \ingroup aerosolskopticalprop
 *	A class used to represent the optical properties of sulphate particles
 *	in the atmosphere. The sulphate particles are represented as log-normal
 *	spheres and use MIE scattering code. Climatologies can be given that
 *	provide mode radius and mode width as a function of position of time. It
 *	is possible to change all of this if needed.
 *
 *	The sulphate refractive index table is only valid from 300 nm to 15 microns.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_AerosolProfileH2SO4 : public skOpticalProperties_AerosolProfile
{
	private:
													skOpticalProperties_AerosolProfileH2SO4 ( const skOpticalProperties_AerosolProfileH2SO4& other );	// Dont allow copy constructor
		skOpticalProperties_AerosolProfileH2SO4&	operator =								( const skOpticalProperties_AerosolProfileH2SO4& other );	// Dont allow assignment operator

	public:
													skOpticalProperties_AerosolProfileH2SO4();
		virtual									   ~skOpticalProperties_AerosolProfileH2SO4();
};


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileDust		2009-5-19*/
/** \ingroup aerosolskopticalprop
 *	A class used to represent the optical properties of dust particles
 *	in the atmosphere. The dust particles are currently represented as
 *	log-normal spheres and uses MIE scattering code. Climatologies can be given that
 *	provide mode radius and mode width as a function of position of time. It
 *	is possible to change all of this if needed.
 *
 *	The dust refractive index table is valid from 200 nm to 15 microns.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_AerosolProfileDust: public skOpticalProperties_AerosolProfileH2SO4
{
	public:
								skOpticalProperties_AerosolProfileDust();
	virtual					   ~skOpticalProperties_AerosolProfileDust() {}
};

/*-----------------------------------------------------------------------------
*					skOpticalProperties_AerosolProfileWater		2017-3-01*/
/** \ingroup aerosolskopticalprop
*	A class used to represent the optical properties of water particles
*	in the atmosphere. The dust particles are currently represented as
*	log-normal spheres and uses MIE scattering code. Climatologies can be given that
*	provide mode radius and mode width as a function of position of time. It
*	is possible to change all of this if needed.
*
*	The water refractive index table is valid from 210 nm to 15.15 microns.
**/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_AerosolProfileWater : public skOpticalProperties_AerosolProfileH2SO4
{
	public:
								skOpticalProperties_AerosolProfileWater();
	virtual					   ~skOpticalProperties_AerosolProfileWater() {}
};


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_AerosolProfileIce		2009-5-19*/
/** \ingroup aerosolskopticalprop
 *	A class used to represent the optical properties of ice particles
 *	in the atmosphere. The dust particles are currently represented as
 *	TMatrix Random particles with a Gamma distribution. Climatologies can be given that
 *	provide gamma effective radius and gamma rate as a function of position of time. It
 *	is possible to change all of this if needed.
 *
 *	The ice refractive index table is valid from 210 nm to 15 microns.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_AerosolProfileIce: public skOpticalProperties_AerosolProfile
{
	private:
												skOpticalProperties_AerosolProfileIce	( const skOpticalProperties_AerosolProfileIce& other );	// Dont allow copy constructor
		skOpticalProperties_AerosolProfileIce&	operator =								( const skOpticalProperties_AerosolProfileIce& other );	// Dont allow assignment operator
	public:
												skOpticalProperties_AerosolProfileIce();
	virtual									   ~skOpticalProperties_AerosolProfileIce() {}
};

/*-----------------------------------------------------------------------------
*					skOpticalProperties_AerosolProfileIce_Mie		2017-03-01*/
/** \ingroup aerosolskopticalprop
*	A class used to represent the optical properties of ice particles
*	in the atmosphere. The dust particles are currently represented as
*	Mie scatterers with a log normal particle size distribution. Climatologies can be given that
*	provide mode radius and mode width as a function of position of time. It
*	is possible to change all of this if needed.
*
*	The ice refractive index table is valid from 210 nm to 15 microns.
**/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_AerosolProfileIce_Mie : public skOpticalProperties_AerosolProfileH2SO4
{
	private:
		skOpticalProperties_AerosolProfileIce_Mie(const skOpticalProperties_AerosolProfileIce& other);	// Dont allow copy constructor
	public:
												skOpticalProperties_AerosolProfileIce_Mie();
	virtual									   ~skOpticalProperties_AerosolProfileIce_Mie() {}
};


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_SimpleRayleigh					2008-4-4*/
/** \ingroup airskopticalprop
 *	\deprecated
 *	The simple Rayleigh scattering cross-section. This is legacy code that
 *	provides backward compatibility with older versions of Sasktarn
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_SimpleRayleigh : public skOpticalProperties
{
	private:
											skOpticalProperties_SimpleRayleigh	( const skOpticalProperties_SimpleRayleigh& other);		// Dont allow Copy constructor
		skOpticalProperties_SimpleRayleigh&	operator =							( const skOpticalProperties_SimpleRayleigh& other);		// Dont allow assignment operator
	public:
											skOpticalProperties_SimpleRayleigh	();
		virtual							   ~skOpticalProperties_SimpleRayleigh	(){};
	public:
		virtual bool						SetAtmosphericState				( skClimatology* neutralatmosphere)								override; // { return true;}			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						SetLocation						( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )		override; // { return true;}			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						CalculateCrossSections			( double wavenumber, double* absxs, double* extxs, double* scattxs)			override;																		//!< Calculate cross-sections at the specified wave-number.
		virtual bool						CalculatePhaseMatrix			( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)	override;				//!< Calculate the phase matrix at the specified wavenumber and scattering angle
		virtual bool						InternalClimatology_UpdateCache		( const GEODETIC_INSTANT& /*pt*/)	override	{ return true;}
		virtual bool						IsScatterer						() const						override    { return true;}
		virtual bool						IsAbsorber						() const						override    { return false;}
        virtual bool                        LegendreCoefficientsPolarized(double wavenumber, double *a1, double *a2, double *a3,
                                                                          double *a4, double *b1, double *b2, int usermaxcoeff,
                                                                           int &opticalmaxcoeff) override;
        virtual bool                        LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) override;
};


/*---------------------------------------------------------------------------
 *					class skOpticalProperties_RayleighDryAir				2003-10-23 */
/**	\ingroup airskopticalprop
 *	This class calculates the corss-sections phase matrix elements for Rayleigh scattering due
 *	to dry air.  It follows the algorithm presented by Bates 1984 and
 *	exactly replicates his cross-section calculations to the 4 significant digits in his Table 1.
 *  Note that this code calculates the cross-section per "air" molecule.  The number
 *	of air molecules per cm3 is simply the total number of all molecules per cm3. The
 *	cross-section is pre-weighted to account for the different gas ratios in the mix.
 *	This class does not cache any phase matrices as its relatively efficient to calculate them
 *	on the fly.
 *
 *	\par
 *	The only difference between my calculation and Bates is that I make an account
 *  for the tiny fraction of gas that is not N2, O2, Ar or CO2.  I assume the small
 * 	residual gas is similar in properties to Argon.  Bates does not describe what
 *	he does with the small residual (it is implied that it is ignored in his paper).
 *
 *	\par Thread Safety
 *	The code was modified in February 2014 to be thread safe for CalculateCrossSections
 *	and CalculatePhaseMatrix.
 *
 *	\par References:
 *	-# Rayleigh Scattering by Air: D.R. Bates. Planet. Space. Sci., 32, 6, 785-790, 1984
 *  -# The Refractive indices and Verdet constants of the inert gases: A. Dalgarno and
 *	   A.E. Kingston., Proc.  Roy. Soc. A., 259, 424-429, 1960.
 *  -# Interferometric determination of the refractive index of carbon dioxide in the
 *     ultra violet region. A. Bideua-Mehu, Y. Guern, R. Adjean and A. Johannin-Gilles.
 *	   Optics Communications, 9, 432-434, 1973.
 */
/*-------------------------------------------------------------------------*/

class skOpticalProperties_RayleighDryAir : public skOpticalProperties
{
	private:
		struct RayleighWavelength_TLS								// Thread safe storage
		{
			double								m_wavenumber;		//!< Thread safe storage for the last wavelength processed by this wavelength 
			double								m_xsection;			//!< Thread safe storage for storing the cross-section from the last calculation
			double								m_delta;			//!< Thread safe storage for the Delta variable
			double								m_deltaprime;		//!< Thread safe storage for the delta prime value
		};

	private:
		double										m_O2mix;			//!< Fraction of O2 in dry air mixture
		double										m_N2mix;			//!< Fraction of N2 in dry air mixture
		double										m_CO2mix;			//!< Fraction of CO2 in dry air mixture
		double										m_Armix;			//!< Fraction of Argon in dry air mixture
		double										m_XXmix;			//!< Fraction of miscellaneous items in the dry air
		        std::map< size_t, RayleighWavelength_TLS>			m_threadstate;
		typedef std::map< size_t, RayleighWavelength_TLS>::iterator	iterator;

	private:
		skOpticalProperties_RayleighDryAir&	operator=							( const skOpticalProperties_RayleighDryAir& other ); // = delete;	// Dont allow assignment operator
											skOpticalProperties_RayleighDryAir	( const skOpticalProperties_RayleighDryAir& other ); // = delete;	// dont allow copy constructor
		bool								LookupUpThreadData					(RayleighWavelength_TLS** data );
		virtual void						AddDepolarization					( double Fk, double fraction, skOpticalProperties_RayleighDryAir::RayleighWavelength_TLS* state) const;
		bool								Interpolate_PhaseMatrixTables		( double mu, skRTPhaseMatrix* P, skOpticalProperties_RayleighDryAir::RayleighWavelength_TLS* state) const;
		virtual bool						CalculateCrossSections				( double wavenum, double* absxs, double* extxs, double* scattxs, skOpticalProperties_RayleighDryAir::RayleighWavelength_TLS* state) const; 

	public:
											skOpticalProperties_RayleighDryAir	();
		virtual							   ~skOpticalProperties_RayleighDryAir	() override;
//		double								Delta								(size_t threadindex) const		{ return m_threadstate[threadindex].m_delta;}
//		double								DeltaPrime							(size_t threadindex) const		{ return m_threadstate[threadindex].m_deltaprime;}
//		double								Depol								(size_t threadindex) const;

	public:
		virtual bool						SetAtmosphericState					( skClimatology* neutralatmosphere )  override;
		virtual bool						SetLocation							( const GEODETIC_INSTANT& pt, bool* crosssectionschanged     )  override;
		virtual bool						InternalClimatology_UpdateCache		( const GEODETIC_INSTANT& pt)                                                                     override;
		virtual bool						CalculateCrossSections				( double wavenumber,  double* absxs, double* extxs, double* scattxs)  override;
        virtual bool						CalculateCrossSectionsArray			( const double * wavenumber, int numwavenumber, double *absxs,  double* extxs, double* scattxs) override;
		virtual bool						CalculatePhaseMatrix				( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)  override;
		virtual bool						IsScatterer                         () const override { return true;}
		virtual bool						CalculateP11                        ( double wavenumber, std::pair<double, size_t> cosscatterandindex, double& p11 ) override;
		virtual bool						IsAbsorber							() const  override    { return false;}
		virtual bool						LegendreCoefficientsP11             (double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) override;
        virtual bool                        LegendreCoefficientsPolarized       (double wavenumber, double* a1, double* a2, double* a3, double* a4, double* b1, double* b2, int usermaxcoeff, int& opticalmaxcoeff) override;

};

/*---------------------------------------------------------------------------
 *					class skOpticalProperties_RayleighDryAir_Inelastic				2003-10-23 */
 /**	\ingroup airskopticalprop
  *	This class calculates separates the elastic (Cabannes) and inelastic (Raman)
  * parts of Rayleigh scattering.
  *
  * CalculateCrossSections and CalculatePhaseMatrix are almost identical to those
  * in skOpticalProperties_RayleighDryAir, except that the King's correction
  * factor has been modified to correspond to only the Cabannes scattering.
  * See the depolarization factors defined in Table 1 of Chance and Spurr. This reduces
  * the scattering cross sections by a few percent, and makes the phase function 
  * slightly less isotropic. The extinction cross section is identical to that returned
  * by skOpticalProperties_RayleighDryAir.
  *
  * CalculateInelasticCrossSections returns the difference between the full 
  * Rayleigh scattering cross section and the Cabannes scattering cross section described 
  * above. The weights of individual lines are computed with Equation (10) from Chance and
  * Spurr for a constant excitation wavelength. Equation (10) is also used for 
  * CalculateInelasticCrossSections_Reverse, but with a constant outgoing wavelength.
  * State and transition data is taken from the supplamentary material to "Spectroscopy 
  * and Radiative Transfer of Planetary Atmospheres" by Chance & Martin 
  * (https://global.oup.com/booksites/content/9780199662104/data/).
  *
  * CalculateInelasticPhaseMatrix calculates the nearly isotropic phase matrix for
  * Raman scattering. According to Equation (2) in Chance and Spurr the phase function
  * is a function of the depolarization ratio, which is 6/7 for Ramann scattering. It 
  * is assumed that this applies to the full phase matrix as well, so that Equations
  * 2.15 and 2.16 from Hansen and Travis can be used with 6/7 as the depolarization ratio.
  * 
  *	\par References:
  *	-# Rayleigh Scattering by Air: D.R. Bates. Planet. Space. Sci., 32, 6, 785-790, 1984
  * -# The Refractive indices and Verdet constants of the inert gases: A. Dalgarno and
  *	   A.E. Kingston., Proc.  Roy. Soc. A., 259, 424-429, 1960.
  * -# Interferometric determination of the refractive index of carbon dioxide in the
  *    ultra violet region. A. Bideua-Mehu, Y. Guern, R. Adjean and A. Johannin-Gilles.
  *	   Optics Communications, 9, 432-434, 1973.
  * -# Light Scattering in Planetary Atmospheres. J.E. Hansen and L.D. Travis. Space
  *    Science Reviews, 16, 527-610, 1974.
  * -# Ring effect studies: Rayleigh scattering, including molecular parameters for
  *    rotational Raman scattering, and the Fraunhofer spectrum. K.V. Chance and R.J.D.
  *    Spurr. Applied Optics, 36(21), 5224-5230, 1997.
  */
  /*-------------------------------------------------------------------------*/

class skOpticalProperties_RayleighDryAir_Inelastic : public skOpticalProperties_Inelastic
{
private:
	struct RayleighWavelength_TLS								// Thread safe storage
	{
		double								m_wavenumber;		//!< Thread safe storage for the last wavelength processed by this wavelength 
		double								m_xsection;			//!< Thread safe storage for storing the cross-section from the last calculation
		double								m_delta;			//!< Thread safe storage for the Delta variable
		double								m_deltaprime;		//!< Thread safe storage for the delta prime value
		bool								m_cabannes;
	};

private:
	double										m_O2mix;			//!< Fraction of O2 in dry air mixture
	double										m_N2mix;			//!< Fraction of N2 in dry air mixture
	double										m_CO2mix;			//!< Fraction of CO2 in dry air mixture
	double										m_Armix;			//!< Fraction of Argon in dry air mixture
	double										m_XXmix;			//!< Fraction of miscellaneous items in the dry air
	std::map< size_t, RayleighWavelength_TLS>			m_threadstate;
	typedef std::map< size_t, RayleighWavelength_TLS>::iterator	iterator;

private:
	double										m_O2partition;		// partition function reciprocal for O2
	double										m_N2partition;		// partition function reciprocal for N2

private:
	skClimatology*								m_backgroundatmosphere;
	double										m_temperaturekelvin;

private:
	bool										LookupUpThreadData(RayleighWavelength_TLS** data);
	virtual void								AddDepolarization(double Fk, double fraction, skOpticalProperties_RayleighDryAir_Inelastic::RayleighWavelength_TLS* state) const;
	bool										Interpolate_PhaseMatrixTables(double mu, skRTPhaseMatrix* P, skOpticalProperties_RayleighDryAir_Inelastic::RayleighWavelength_TLS* state) const;
	virtual bool								CalculateCrossSections(double wavenum, double* absxs, double* extxs, double* scattxs, skOpticalProperties_RayleighDryAir_Inelastic::RayleighWavelength_TLS* state, bool cabannes) const;

	bool										CalculatePartitionFunctions( );
	double										CalculatePolarizabilityAnisotropyN2(double excitationwavenumber);
	double										CalculatePolarizabilityAnisotropyO2(double excitationwavenumber);

	bool										CalculateInelasticCrossSections_Classical(double wavenumin, double* inelxs);
	bool										CalculateInelasticCrossSections_Classical(double wavenumin, size_t lineidx, double* wavenumout, double* inelxs);
	bool										CalculateInelasticCrossSections_Quantum(double wavenumin, double* inelxs);
	bool										CalculateInelasticCrossSections_Quantum(double wavenumin, size_t lineidx, double* wavenumout, double* inelxs);

	bool										CalculateInelasticCrossSections_ClassicalReverse(double wavenumout, double* inelxs);
	bool										CalculateInelasticCrossSections_ClassicalReverse(double wavenumout, size_t lineidx, double* wavenumin, double* inelxs);
	bool										CalculateInelasticCrossSections_QuantumReverse(double wavenumout, double* inelxs);
	bool										CalculateInelasticCrossSections_QuantumReverse(double wavenumout, size_t lineidx, double* wavenumin, double* inelxs);


public:
	struct sk_raman_partition_data
	{
		int N;
		int J;
		int gN;
		double Eterm;
	};

	struct sk_raman_xs_data
	{
		double deltaE;
		int J;
		int gN;
		double Eterm;
		double cPT;
	};

public:
										skOpticalProperties_RayleighDryAir_Inelastic();
	virtual							   ~skOpticalProperties_RayleighDryAir_Inelastic() override;
	//		double								Delta								(size_t threadindex) const		{ return m_threadstate[threadindex].m_delta;}
	//		double								DeltaPrime							(size_t threadindex) const		{ return m_threadstate[threadindex].m_deltaprime;}
	//		double								Depol								(size_t threadindex) const;

public:
	virtual bool						SetAtmosphericState(skClimatology* neutralatmosphere)  override;
	virtual bool						SetLocation(const GEODETIC_INSTANT& pt, bool* crosssectionschanged)  override;
	virtual bool						InternalClimatology_UpdateCache(const GEODETIC_INSTANT& pt)  override;
	virtual bool						CalculateCrossSections(double wavenumber, double* absxs, double* extxs, double* scattxs) override;																//!< Calculate cross-sections at the specified wave-number.
	virtual bool						CalculatePhaseMatrix(double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix)  override;
	virtual bool						IsScatterer() const  override { return false; }
	virtual bool						IsAbsorber() const  override { return false; }

	virtual size_t						NumInelasticLines() const override;
	virtual bool						CalculateInelasticCrossSections(double wavenumin, double* inelxs) override;
	virtual bool						CalculateInelasticCrossSections(double wavenumin, size_t lineidx, double* wavenumout, double* inelxs) override;
	virtual bool						CalculateInelasticCrossSections_Reverse(double wavenumout, double* inelxs) override;
	virtual bool						CalculateInelasticCrossSections_Reverse(double wavenumout, size_t lineidx, double* wavenumin, double* inelxs) override;
	virtual bool						CalculateInelasticPhaseMatrix(double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix) override;
};

/*-----------------------------------------------------------------------------
 *					class skOpticalProperties_HenyeyGreenstein		2008-11-17*/
/** \ingroup aerosolskopticalprop
 *	A class that calculates the Henyey-Greenstein phase matrices.  The calculation
 *	is only scalar and does not provide any polarization: only P(1,1) of the phase
 *	matrix is non-zero.
 * 
 *	\par Support for Scattering Cross Section
 *	The Henyey Greenstein algorithm does not specify the scattering cross-section. This
 *	implentation provides this information by extracting extinction, scattering and
 *	cross-section information from another instance of class skOpticalProperties. The
 *	user is free to choose the instance and may, for example, choose a user defined table
 *	or the Mie scattering cross-sections which when coupled to this object provide HG phase matrices.
 *	It is important that the the supplied instance returns true for the "IsScatterer" method otherwise
 *	if this calculation is to have any effect on the radiative transfer calculation.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_HenyeyGreenstein : public skOpticalProperties
{
	private:
		skOpticalProperties*					m_extinction;		//!< Object that calculate sthe extinction.
		double									m_g;				//!< The Henyey-Greenstein aymmetry factor.

	private:
//		bool									DeepCopy							( const skOpticalProperties_HenyeyGreenstein& other );
												skOpticalProperties_HenyeyGreenstein( const skOpticalProperties_HenyeyGreenstein& other );	// Dont allow copy constructor
		skOpticalProperties_HenyeyGreenstein&	operator=							( const skOpticalProperties_HenyeyGreenstein& other );	// Dotn allow assignment operator

	public:
												skOpticalProperties_HenyeyGreenstein	();
		virtual								   ~skOpticalProperties_HenyeyGreenstein	();
		bool									SetCrossSectionObject			( skOpticalProperties* object );
		bool									Set_HG_AsymmetryFactor			( double g );

	public:
		virtual bool							SetAtmosphericState				( skClimatology* neutralatmosphere) override; 					//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool							SetLocation						( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  ) override; 					//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool							InternalClimatology_UpdateCache		( const GEODETIC_INSTANT& pt)  override;
		virtual bool							IsScatterer						() const  override;
		virtual bool							IsAbsorber						() const  override;
		virtual bool							CalculateCrossSections			( double wavenumber, double* absxs, double* extxs, double* scattxs)  override;																//!< Calculate cross-sections at the specified wave-number.
		virtual bool							CalculatePhaseMatrix			( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix) override;
};




/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_GaussQuadrature		2009-11-10*/
/** \internal
 *	A helper class used to integrate the absorption, extinction and scattering 
 *	cross-sections across a guassian point spread function using a 50 point
 *	gaussian quadrature.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_Convolved_GaussQuadrature 
{
	private:
		nxGaussQuadratureBase			m_gauss;

		double							m_lambda0;
		double							m_sigma;
		double							m_normalizer;
		size_t							m_numsteps;

		double							m_absorptioncm2;
		double							m_extinctioncm2;
		double							m_scattercm2;

	private:
		bool							Initialize					( double lambda0, double instpsf_fwhm, double quadrature_resolution_nm  );
		bool							GetIntegrationRange			( double* minwavelen_nm, double* maxwavelen_nm );
		double							Distribution				( double lambda );
		bool							UpdateCrossSectionSummation	( double lambda, double weight, skOpticalProperties* highresextinction, double*sumcheck);
		bool							ConvolveTrapezoid			( skOpticalProperties* highresextinction, double* sumcheck);
		bool							ConvolveGaussQuad			( skOpticalProperties* highresextinction, double* sumcheck);

	public:
										skOpticalProperties_Convolved_GaussQuadrature();
		double							Extinctioncm2	() const { return m_extinctioncm2;}
		double							Absorptioncm2	() const { return m_absorptioncm2;}
		double							Scattercm2		() const { return m_scattercm2;}
		bool							Convolve		( double lambda0, double instpsf_fwhm, double quadratureresolution_nm, bool usegaussianquadrature, skOpticalProperties* highresextinction);					// Integrate when a function can evaluate Y at location X
};

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved		2009-11-6*/
/** \ingroup skopticalpropmisc
 *	Note this class was one of the first convolving classes developed and has been
 *	superceeded by other classes. Take a look at class #skOpticalProperties_ConvolvedDiscreteWavelenCachedState
 *	if you are planning to use measured cross-sections as that class has much tighter integration
 *	that is not available in this class (e.g automatically adjusting PSF with wavelength and temperature)
 *		
 *	This class is still useful if you just want to tinker around. It will take a
 *	high resolution cross-section object and convolve it with a guassian point
 *	spread function on the the fly.  The code convolves the extinction, absorption
 *	and scattering cross-sections using a point spread function width defined in a
 *	derived class through virtual method Get_FWHM.  This class is nice if you want to
 *	experiment with the quadrature.
 *
 *	The class does not attempt to generate convolutions of phase matrices. The
 *	phase matrix calculation is performed using the "high resolution" cross-section 
 *	object only at the central wavelength of interest (i.e. a monochromatic calculation).
 *
 *	\par See Also
 *	 #skOpticalProperties_ConvolvedDiscreteWavelenCachedState
 *	 #skOpticalProperties_Convolved_FixedFWHM
 *	 #skOpticalProperties_Convolved_GaussQuadrature
 **/
/*---------------------------------------------------------------------------*/

class  skOpticalProperties_Convolved : public skOpticalProperties
{
	private:
		std::map< size_t, skOpticalProperties_Convolved_GaussQuadrature>			m_quadrature;			//!< We need one quadrature per thread 
		typedef std::map< size_t, skOpticalProperties_Convolved_GaussQuadrature>::iterator	iterator;		//!< We need one quadrature per thread 
		skOpticalProperties*							m_highresextinction;			//!< The object that supplies the high resolution cross-section spectra.
		double											m_quadratureresolution_nm;		//!< The resolution (in nanometers) to be used for the quiadrarture integration (set by user)
		bool											m_usegaussianquad;				//!< True if we want to use guassian quadrarture rather than trapezoidal (to integrate over the gaussian point spread function)
		skClimatology*									m_backgroundatmosphere;

	private:
		bool											LookupUpThreadData					(skOpticalProperties_Convolved_GaussQuadrature** data );
		void											CommonInit();
		bool											CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs, skOpticalProperties_Convolved_GaussQuadrature* quadrature ) const;
														skOpticalProperties_Convolved		( const skOpticalProperties_Convolved& other ); // Dont allow copy constructor
		skOpticalProperties_Convolved&					operator =							( const skOpticalProperties_Convolved& other );	// Dont allow assignment operator

	public:
														skOpticalProperties_Convolved		();
		virtual										   ~skOpticalProperties_Convolved		();
		bool											SetHighResolutionCrossSectionObject	( skOpticalProperties* extinctionobject );
		bool											SetQuadratureResolution			    ( double quadresolution_nm )	{ m_quadratureresolution_nm = quadresolution_nm; return (m_quadratureresolution_nm > 0);}
		bool											UseGaussianQuadrature				( bool usegauss )				{ m_usegaussianquad = usegauss; return true;}
	
	public:
		virtual	bool									Get_FWHM							( double lambda_nm, double* fwhm_nm ) const = 0;

	public:
		virtual bool									SetAtmosphericState					( skClimatology* neutralatmosphere)  override;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool									SetLocation							( const GEODETIC_INSTANT& pt,  bool* crosssections_changed )  override;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool									InternalClimatology_UpdateCache			( const GEODETIC_INSTANT& pt)  override;
		virtual bool									CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs)  override;			//!< Calculate cross-sections at the specified wave-number.
		virtual bool									IsScatterer							() const  override;						//!< Returns true if this particles scatters radiation
		virtual bool									IsAbsorber							() const  override;						//!< Returns true if this particles absorbs radiation radiation
		virtual bool									CalculatePhaseMatrix				( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)  override;
};

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_FixedFWHM				2009-11-10*/
/** \ingroup skopticalpropmisc
 *	A derived class that allows users to generate convolved cross-sections
 *	from high resolution cross-sections using a Guassian PSF with a constant
 *	(with wavelength) FWHM. The full width half max is specified in nanometers.
 **/
/*---------------------------------------------------------------------------*/

class  skOpticalProperties_Convolved_FixedFWHM : public skOpticalProperties_Convolved
{
	private:
		double										m_fwhm_nm;								//!< The full width half max expressed in nm.

	private:
		bool										DeepCopy								( const skOpticalProperties_Convolved_FixedFWHM& other);
													skOpticalProperties_Convolved_FixedFWHM	( const skOpticalProperties_Convolved_FixedFWHM& other);		// Dont allow copy constructor
		skOpticalProperties_Convolved_FixedFWHM&	operator =								( const skOpticalProperties_Convolved_FixedFWHM& other);		// Dont allow assignment
	public:
													skOpticalProperties_Convolved_FixedFWHM	( );
													skOpticalProperties_Convolved_FixedFWHM	( skOpticalProperties* highresextinction, double fwhm_nm, double quad_resol_nm );
		virtual									   ~skOpticalProperties_Convolved_FixedFWHM	( ) override;
		bool										SetFWHM									( double fwhm ) { m_fwhm_nm = fwhm; return (fwhm > 0.0);}
		virtual	bool								Get_FWHM								( double lambda_nm, double* fwhm_nm ) const override;
};



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra		2009-11-10*/
/** \ingroup skopticalpropmisc
 *	A class that creates absoprption spectra for one molecule by combining
 *	several cross-section objects togther. This was built so we could use
 *	one set of O3 cross-sections which only cover one part of the spectrum
 *	and supplement other regions with another cross-section object.
 *	It exists and that is about all that can be said as it is not used very
 *	frequently.
 **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_MultipleOverlappingSpectra : public skOpticalProperties
{
	private:
		class entry_struct
		{
			public:
				double						m_lowerwavelength;
				double						m_upperwavelength;
				skOpticalProperties*		m_crosssection;
			
			private:
				void						ReleaseResources();

			public:
											entry_struct();
										   ~entry_struct();
			   bool							SetValues	 ( double lowerwavelength_nm, double upperwavelength_nm, skOpticalProperties*	crosssection);
		};

		class entry_locator
		{
			private:
				double	m_wavelen;

			public:
											entry_locator( double wavelen )																{ m_wavelen = wavelen;}
				bool						operator()	 ( const skOpticalProperties_MultipleOverlappingSpectra::entry_struct& entry ) const	{ return (m_wavelen >= entry.m_lowerwavelength) && (m_wavelen < entry.m_upperwavelength);}
		};




	private:
		        std::list< entry_struct >					m_entries;		
		typedef std::list< entry_struct >::iterator			iterator;
		typedef std::list< entry_struct >::const_iterator	const_iterator;

	private:
		bool												FindBoundingEntry								( double wavenumber, entry_struct** entry);
															skOpticalProperties_MultipleOverlappingSpectra	( const skOpticalProperties_MultipleOverlappingSpectra& other );	// dont allow copy constructor
		skOpticalProperties_MultipleOverlappingSpectra&		operator =										( const skOpticalProperties_MultipleOverlappingSpectra& other );	// dont allow assignment operator

	public:
															skOpticalProperties_MultipleOverlappingSpectra	();
		virtual											   ~skOpticalProperties_MultipleOverlappingSpectra	();
		bool												AddEntry										( double startlambda_nm, double endlambda_nm, skOpticalProperties* crosssection );

	public:
		virtual bool										SetAtmosphericState					( skClimatology* neutralatmosphere)  override;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool										SetLocation							( const GEODETIC_INSTANT& pt,  bool* crosssections_changed )  override;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool										InternalClimatology_UpdateCache			( const GEODETIC_INSTANT& pt)  override;
		virtual bool										CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs)  override;			//!< Calculate cross-sections at the specified wave-number.
		virtual bool										CalculatePhaseMatrix				( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix )  override;
		virtual bool										IsScatterer							() const  override;						//!< Returns true if this particles scatters radiation
		virtual bool										IsAbsorber							() const  override;						//!< Returns true if this particles absorbs radiation radiation
};
