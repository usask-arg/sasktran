


/*-----------------------------------------------------------------------------
 *					class SKTRAN_TableSolarTransmissionProfile_V21		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableSolarTransmissionProfile_V21
{
	private:
		std::shared_ptr< const SKTRAN_RayFactory_Base>				m_rayfactory;
		double														m_cossza;
		double														m_slon;
		std::shared_ptr<const SKTRAN_GridDefRayTracingShells_V21>	m_heights;
		std::vector< SKTRANSO_RayInternalGeometry >					m_rays;	

	private:
		void										ReleaseResources	();
		bool										AllocateRays		( size_t numrays );

	public:
													SKTRAN_TableSolarTransmissionProfile_V21	();
												   ~SKTRAN_TableSolarTransmissionProfile_V21	();
													SKTRAN_TableSolarTransmissionProfile_V21	( const SKTRAN_TableSolarTransmissionProfile_V21&  moveother);
		void										operator =									( const SKTRAN_TableSolarTransmissionProfile_V21&  moveother);
		
		double										CosSZA										() const				{ return m_cossza;}
		double										SLON										() const				{ return m_slon;}
		double										Height										(size_t hindex) const	{ return m_heights->At( (SKTRAN_GridIndex)hindex);}
		size_t										NumRays										() const				{ return m_rays.size();}
		SKTRANSO_RayInternalGeometry*				RayAtVar									( size_t hidx );
		const SKTRANSO_RayInternalGeometry*			RayAt										( size_t hidx ) const;
		const SKTRAN_GridDefRayTracingShells_V21*	HeightsPtr									() const { return m_heights.get();}
		bool										CreateProfile								( size_t profileidx, const SKTRAN_SpecsInternal_V21* specs );

};

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableSolarTransmission		2010-2-19*/
/** @ingroup LegacySO
 *	The parent class for all Solar Transmission tables implemented in
 *	Sasktran. It is typically used by the internal diffuse rays and
 *	ground points to calculate the optical depth of the sun
 *	at any location in the atmopshere. It is used calculate single scatter terms.
 *
 *	This is a purely virtual class and is really an interface definition.
 *	The Solar transmission table is used calculate the optical depth of the sun
 *	at any location in the atmopshere. The geometry aspect sets up all of the
 *	various wavelength independent aspects, (e.g. ray tracing and interpolation).
 *
 *	This class is constant when wavelength independent ray geometry is used and
 *	one instance of this class can be shared between many simultaneous threads.
 *
 *	\par See Also
 *	#SKTRAN_TableSolarTransmissionOptical_V2
**/
/*---------------------------------------------------------------------------*/

class SKTRANSO_TableSolarTransmission : public SKTRANSO_JindexTableBase
{
	private:
		SKTRAN_TableSolarTransmissionOptical_V21*				m_opticaltable;
		std::vector<SKTRAN_TableSolarTransmissionProfile_V21>	m_profiles;						//!< The array of solar transmisison profiles
		std::vector<size_t>										m_pointsindex;					//!< The starting point index for each profile
		size_t													m_numpoints;
		const SKTRAN_GridDefCosSZA_V21*							m_cosszagrid;					//!< The cos(SZA) location of each profile
		const SKTRAN_GridDefSLON_V21*							m_slongrid;						//!< The SLON location of each profile
		bool													m_notrequiredforsinglescatter;

	private:
		void									ReleaseResources				();
		bool									FetchSLONandSZAGrids			(const SKTRAN_SpecsInternal_Diffuse_V21*    diffusespecs );
		bool									PointIndexToProfileAndHeight	( size_t pointindex, size_t* profileindex, size_t* heightindex) const;
		bool									AllocateProfiles				( size_t numprofiles, const SKTRAN_SpecsInternal_V21* specs );
		bool									ConfigureGeometry_Stage1		( const SKTRAN_SpecsInternal_V21* specs );
		bool									ConfigureProfileIndexing		( );


	protected:
		bool										ProfileSZAandSLON			( size_t profileidx, double* sza, double* slon) const;
		const SKTRAN_GridDefRayTracingShells_V21*	ProfileHeightsPtr			( size_t profileidx ) const;

	public:
												SKTRANSO_TableSolarTransmission();
		virtual								   ~SKTRANSO_TableSolarTransmission();


		size_t									NumPoints						() const { return m_numpoints;}

		bool									RayAtPoint						( size_t pointindex, const SKTRANSO_RayInternalGeometry** ray ) const;

		bool									NotRequiredForSingleScatter		() const				{ return m_notrequiredforsinglescatter;}
		const SKTRAN_GridDefCosSZA_V21*			SZAGrid							() const				{ return m_cosszagrid;}
		const SKTRAN_GridDefSLON_V21*			SLONGrid						() const				{ return m_slongrid;}

		bool									IsDefined						() const				{ return (m_numpoints > 0);}


		bool									ConfigureGeometry_Stage2MT		( size_t pointindex,
																				  SKTRANSO_Quadrature_TLS_V21*				quadrature,
																				  const SKTRAN_SpecsInternal_V21*		specs
																				);

		bool									AttachOpticalToGeometry			(  );

		bool									ConfigureOptical				( bool										singlescatter,
																				  const SKTRAN_TableOpticalProperties_V21*	optprop,
																				  SKTRAN_ThreadManager *					threadmanager
																				);

		bool									ConvertJindexToPosition			( const SKTRANSO_JIndex* entry,
																						  size_t* positionindex
																				) const;

		bool									ConfigureOptical_Stage2MT		( size_t pointindex,
																				  SKTRANSO_Quadrature_TLS_V21* threadquadrature
																				);

		virtual const SKTRAN_StokesScalar*		ConvertJIndexToRadiancePtr		( const SKTRANSO_JIndex*		entry,
																				  ENUM_SKTRAN_JSOURCE			jsource
																				) const override;
	public:


		virtual bool							ConfigureGeometry				( const SKTRAN_SpecsInternal_V21*		specs,
																				  SKTRAN_ThreadManager*					threadmanager
																				 );


};



/*-----------------------------------------------------------------------------
 *					class SKTRANSO_TableGroundPointDiffuse		2010-2-19*/
/** @ingroup LegacySO
 *	The diffuse ground point table calculates the diffuse radiance emanating
 *	from the ground/lowest layer. The basic premise of the class is to distribute
 *	a set of points across the ray tracing region and to calculate
 *	the downward and upward flux at each point. Rays traced in the model can then determine
 *	the ground scatter signal by interpolating the upward flux at each
 *	ground point (as well as applying a simplified albedo BRDF).  Note that
 *	single scatter ground calculations are not handled by this class.
 *
 *	The biggest complication for this class is that we must match each ground point to
 *	a diffuse point on a one-to-one basis. This is done so we can easily calculate the
 *	incoming, downward, radiance over the hemisphere at a point. The original
 *	Sasktran versions explicitly tied the diffuse ground points to the bottom point of 
 *	each diffuse profiles in teh diffuse points table.  While this is okay and is still done
 *	it might make more sense in some scenarios to scatter a few extra 
 *	diffuse ground points around so we can better handle diffuse rays that hit the ground
 *	a long way from any of the diffuse profile.
 *
 *	I have tried to move all of the code that explicitly ties the ground points to
 *	diffuse profiles in the table of diffuse points: that connection is now made 
 *	in user specific implementations outside of the main engine project. In the future
 *	we can develop ground point tables that provide their own internal distribution of
 *	diffuse points.
 *
 *	\par Internal Jindex Tables to Calculate Upward Flux
 *	The table of diffuse ground points needs to calculate the upward flux at all
 *	of the points in the table. This requires an integral of the downward component
 *	of radiance over the \f$2\pi \f$ unit hemisphere at each point in the table. This integral
 *	is performed by the SKTRAN_Quadrature_TLS object in a multi-threaded environment. The
 *	integral is performed in the following stages:
 *
 *		- ConfigureGeometry
 *		- CreateJIndexTables
 *		- ScatterIncomingRadiance
 *
 *	\par Interpolation of Upward Flux in tables
 *	The tables primary function is to provide the estimate of diffuse upward flux at any point
 *	(on the ground) in the ray tracing region. This functionality is provided via methods
 *	
**/
/*---------------------------------------------------------------------------*/

class SKTRANSO_TableGroundPointDiffuse : public SKTRANSO_JindexTableBase
{
	private:
		SKTRAN_GroundPointDiffuseGeometry_V21*			m_groundpoints;			
		size_t											m_numgroundpoints;
		SKTRAN_TableGroundPointDiffuseOptical_V21*		m_opticaltable;

	private:
		void											ReleaseResources					();
	
	public:
														SKTRANSO_TableGroundPointDiffuse	();
		virtual										   ~SKTRANSO_TableGroundPointDiffuse	();

		bool											ConfigureGeometry					( const SKTRAN_SpecsInternal_V21*			specs,
																							  SKTRAN_ThreadManager*						threadmanager
																							);

		bool											AttachOpticalToGeometry				();

		bool											ConfigureOptical					( bool										singlescatter,
																							  const SKTRAN_TableOpticalProperties_V21*	opticalprops,
																							  SKTRAN_ThreadManager*						threadmanager );

		SKTRAN_GroundPointDiffuseGeometry_V21*			PointAtVar							(size_t idx);
		const SKTRAN_GroundPointDiffuseGeometry_V21*		PointAt								(size_t idx) const;
		size_t											NumPoints							()	const			{ return m_numgroundpoints;}


		bool											CreateJIndexTables_HemisphereIntegral			( SKTRAN_ThreadManager* threadmanager);
		bool											ScatterIncomingRadiance							( SKTRAN_ThreadManager* threadmanager );
		bool											ScatterIncomingRadiance_MT						( size_t pointindex );
		bool											ConfigureGeometry_Stage1						( const SKTRAN_SpecsInternal_V21*	specs );
		bool											ConfigureOptical_Stage2MT						( size_t pointindex, const SKTRAN_TableOpticalProperties_V21* optprop, SKTRANSO_Quadrature_TLS_V21* threadquadrature );
		virtual const SKTRAN_StokesScalar*				ConvertJIndexToRadiancePtr						( const SKTRANSO_JIndex* entry, ENUM_SKTRAN_JSOURCE jsource ) const;

	public:

		virtual size_t		NumberOfGroundPointsRequired		( const SKTRAN_SpecsInternal_V21*	specs
																) const = 0;


		virtual bool		ConfigureGeometry_Stage2MT			( size_t									pointindex,
																  SKTRANSO_Quadrature_TLS_V21*					threadquadrature
																) = 0;
};


/*-----------------------------------------------------------------------------
 *					class SKTRANSO_TableDiffusePoints		2010-2-19*/
/** @ingroup LegacySO
**/
/*---------------------------------------------------------------------------*/

class SKTRANSO_TableDiffusePoints: public SKTRANSO_JindexTableBase
{
	private:
		std::vector<SKTRAN_DiffuseProfileGeometry_V21*>			m_profiles;						//!< The array [m_numprofiles] of diffuse profiles in this table
		std::vector<size_t>										m_pointsindex;					//!< The starting point index for each profile
		const SKTRAN_GridDefSLON_V21*							m_slongrid;						//!< The solar longitude (in radians) of each diffuse profile (often all 0)
		const SKTRAN_GridDefCosSZA_V21*							m_cosszagrid;					//!< The cosine of solar zenith angle of each diffuse profile
		SKTRAN_TableDiffusePointsOptical_V21*					m_opticaltable;


	private:
		void											ReleaseResources				();
		bool											AllocateProfiles				( size_t numprofiles );
		bool											DeleteProfiles					();

		
		bool											FetchSLONandSZAGrids			( const SKTRAN_SpecsInternal_Diffuse_V21*    diffusespecs );
	
		bool											CreateGrid_Stage1				( std::shared_ptr<const SKTRAN_CoordinateTransform_V2> 		  coords, 
																						  const SKTRAN_SpecsInternal_RayTracing_V21* raytracingspecs,
																						  const SKTRAN_SpecsInternal_Quadrature_V21* quadraturespecs,
																						  const SKTRAN_SpecsInternal_Diffuse_V21*    diffusespecs
																						);



	public:
														SKTRANSO_TableDiffusePoints	();
		virtual										   ~SKTRANSO_TableDiffusePoints	();

		bool											IsDefined								() const { return !m_profiles.empty();}
		size_t											NumProfiles								() const { return m_profiles.size();}
		const SKTRAN_GridDefSLON_V21*					SLonGrid								() const { return m_slongrid;}
		const SKTRAN_GridDefCosSZA_V21*					CosSZAGrid								() const { return m_cosszagrid;}
		size_t											NumDiffusePoints						() const { return m_pointsindex[m_profiles.size()];}
		bool											PointIndexToProfileAndHeight			( size_t pointindex, size_t* profileindex, size_t* heightindex ) const;
		const SKTRAN_TableDiffusePointsOptical_V21*		OpticalTable							() const { return m_opticaltable;}
		SKTRAN_TableDiffusePointsOptical_V21*			OpticalTableVar							()       { return m_opticaltable;}
		const SKTRAN_DiffusePointGeometry_V21*			PointAt									( size_t profileidx, size_t radiusindex )	const;
		const SKTRAN_DiffuseProfileGeometry_V21*			ProfileAt								( size_t profileidx ) const { return m_profiles[profileidx];}
		SKTRAN_DiffuseProfileGeometry_V21*				ProfileAtVar							( size_t profileidx )       { return m_profiles[profileidx];}
		bool											LookupDiffusePoint						( size_t pointindex, SKTRAN_DiffusePointGeometry_V21** point );
		bool											CreateJIndexTables_RayIntegral			( SKTRAN_ThreadManager* threadmanager  );
		bool											AttachOpticalToGeometry					( );




	// --- virtual that this class creates.

	public:
		virtual size_t									MaxInterpolatedPointsPerJ				() const = 0;

		virtual bool									ConfigureGeometry						( const SKTRAN_SpecsInternal_V21*				specs,
																								  SKTRAN_ThreadManager*							threadmanager );

		virtual bool									ConfigureOptical						( bool										singlescatter,
																								  const SKTRAN_TableOpticalProperties_V21*	optprop,
																								  SKTRAN_ThreadManager *					threadmanager
																								);

		virtual const SKTRAN_StokesScalar*				ConvertJIndexToRadiancePtr				( const SKTRANSO_JIndex*	entry,
																								  ENUM_SKTRAN_JSOURCE		jsource
																								) const;

		virtual bool									ConfigureGeometry_Stage2MT				( size_t						pointindex,
																								  SKTRANSO_Quadrature_TLS_V21*		threadquadrature
																								);
};


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS		2010-6-18*/
/** A class that creates a single scatter table for a line of sight ray.
 *	The purpose is to provide a quick and efficient single scatter method
 *	for lines of sight rays. This allows us to skip initialization of the
 *	solar transmission table for pure single scatter calculations.
 *
 *	This strategy works well when there are only a few lines of sight (< 100?).
 *	Each line of sight typically generates between 50 and 100 rays to the sun while
 *	the solar transmission table creates around 5000 - 10000 lines of sight to the sun.
 *	The single scatter speed advantage disappears if the array of lines of sight
 *	(from the observer) generates more rays to the sun than the solar transmission table.
 *	
 **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_TableRayLOS : public SKTRANSO_JindexTableBase
{
	private:
		std::weak_ptr< const SKTRAN_RayFactory_Base>						m_losrayfactory;
		SKTRANSO_RayLOSGeometry_V21*										m_ray;									// The ray that we are tabulating
		SKTRAN_GridDefBase_V2												m_distancefromobserver;					// The grid specifying distancefromobserver	[m_numraystosun] of distance of rays to sun from observer									
		std::vector< SKTRANSO_RayInternalGeometry> 							m_raystosun;							// The array m_raystosun					[m_numraystosun    ] of rays to the sun
		std::vector< SKTRAN_StokesScalar >									m_incomingsolarradiance;				// The array m_incomingsolarradiance		[m_numraystosun    ] of incoming solar radiance (one for each element of m_raystosun

	private:
		bool											ConfigureEntry				( size_t						index,
																					  double						s,
																					  SKTRANSO_Quadrature_TLS_V21*		quadrature
																					);

		bool												AllocateStorage					( size_t numpoints);

		void												ReleaseStorage					();

		size_t												NumRaysToSun					() const { return m_raystosun.size();}
		const std::weak_ptr<const SKTRAN_RayFactory_Base>&	LOSRayFactory					() const { return m_losrayfactory;}



	protected:
		bool										ConfigureRayToSunLocations		( std::vector< SKTRAN_Distance>&	locations,
																					  SKTRANSO_Quadrature_TLS_V21*					quadrature,
																					  SKTRANSO_RayLOSGeometry_V21*				losray
																					);

	public:	
													SKTRANSO_TableRayLOS( std::weak_ptr< const SKTRAN_RayFactory_Base> rayfactory );
		virtual 								   ~SKTRANSO_TableRayLOS();

		virtual	 bool								ConfigureGeometry				( SKTRANSO_Quadrature_TLS_V21* quadrature,	 SKTRANSO_RayLOSGeometry_V21* losray ) = 0;

		bool										ConfigureOptical				( SKTRANSO_Quadrature_TLS_V21 *	   quadrature );

		virtual bool								InterpolateTable				( const HELIODETIC_POINT&			location,
																					  const HELIODETIC_UNITVECTOR&		look,
																						    SKTRANSO_JIndex*			vertexdescriptortable, 
																						    size_t						maxpts, 
																						    size_t*						npts,
																						    double						weight ) const override;

		virtual const SKTRAN_StokesScalar*			ConvertJIndexToRadiancePtr		( const SKTRANSO_JIndex*		entry,
																					  ENUM_SKTRAN_JSOURCE			jsource
																					) const override;




};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_TableRayLOSFactory					2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableRayLOSFactory : public nxUnknown
{
	public:
								SKTRAN_TableRayLOSFactory() {};
		virtual				   ~SKTRAN_TableRayLOSFactory() {}
		virtual bool			CreateInternalSingleScatterTable	(  SKTRANSO_TableRayLOS** internalsinglescattertable, std::weak_ptr< const SKTRAN_RayFactory_Base> rayfactory ) const = 0;
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableScatteringMatrixOptical_V2		2010-3-2*/
/** **/
/*---------------------------------------------------------------------------*/

/*
class SKTRAN_TableScatteringMatrixOptical_V2
{
	public:
		virtual bool										AttachToGeometry()							  = 0;
		virtual SKTRAN_ScatteringMatrixProfileOptical_V21*	ScatteringProfileAtDiffuseIndex( size_t	idx ) = 0;
		
		virtual bool										LookupScatterPoint	( size_t pointindex, SKTRAN_ScatteringMatrixPointOptical_V21** point ) = 0;
		virtual bool										NumPoints			( ) const = 0;

};
*/

/*-----------------------------------------------------------------------------
 *					class SKTRAN_TableScatteringMatrixGeometry_V2		2010-3-2*/
/** The scattering matrix table stores all of the information used to scatter
 *  incoming rays to outbound source functions for each point in the diffuse
 *	points table. The scattering matrix table is so closely coupled to the
 *	diffuse points table that it is internally created by a derived instance
 *	of the diffuse points table.
 *
 *	The purpose of the class is to provide the incoming ray drirections and
 *	Derived instances of the class are expecinitialize
 **/
/*---------------------------------------------------------------------------*/

/*
class SKTRAN_TableScatteringMatrixGeometry_V2: public nxUnknown
{
	public:
		enum OPTIONSENUM { OPTION_SAME_OUTBOUND_SPHERE_EVERYWHERE = 1,
						   OPTION_SAME_INCOMING_RAYS_EVERYWHERE   = 2 };

		virtual bool		CreateOpticalTable	 ( SKTRAN_TableScatteringMatrixOptical_V2**  opticaltable
												 ) const = 0;


		virtual bool		IsOptionTrue		 ( SKTRAN_TableScatteringMatrixGeometry_V2::OPTIONSENUM options
												 ) const = 0;


		virtual bool		LookupScatteringGrids( size_t                                          profileindex,		//!< Profile index of a point in the diffuse points table
												   size_t                                          heightindex,			//!< Height  index of a point in the diffuse points table
												   const SKTRAN_GridDefDiffuseIncomingZenith_V21**  incomingzenith,		//!< return the incoming zenithangles
												   const SKTRAN_GridDefDiffuseIncomingAzimuth_V21** incomingazimuth,		//!< return the incoming azimuth angles
												   const SKTRAN_UnitSphere_V2**	   unitsphere							//!< return the outbound unit sphere directions
												 ) const = 0;

		virtual bool		LookupScatterPoint	( size_t pointindex, SKTRAN_ScatteringMatrixPointGeometry_V21** point ) = 0;
		virtual bool		NumPoints			( ) const = 0;
};
*/
