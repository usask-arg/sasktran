

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePointsOptical_V21		2007-12-18*/
/** The table of diffuse profiles for a specific wavelength **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableDiffusePointsOptical_V21
{
	private:
		SKTRANSO_TableDiffusePoints*					m_geometry;						//!< The geometry table that created this instance (do not modify after constructor)
		std::vector<SKTRAN_DiffuseProfileOptical_V21* >	m_profiles;						//!< The array of diffuse profiles

	private:
		void											ReleaseDiffuseArray();
		bool											AllocateDiffuseArrays ();

	public:
														SKTRAN_TableDiffusePointsOptical_V21	( SKTRANSO_TableDiffusePoints*	geometry );
		virtual										   ~SKTRAN_TableDiffusePointsOptical_V21	();
		bool											AttachToGeometry						( );
		const SKTRAN_DiffuseProfileOptical_V21*			ProfileAt								( size_t idx )	const { return m_profiles[idx];}
		const SKTRAN_DiffusePointOptical_V21*			PointAt									( size_t profileidx, size_t radiusindex )	const {return ProfileAt(profileidx)->At(radiusindex);}
		size_t											NumProfiles								() const { return m_profiles.size();}		//!< The number of diffuse profiles;
//		const SKTRAN_StokesScalar*						LookUpJindex							( const SKTRANSO_JIndex* index, ENUM_SKTRAN_JSOURCE tablesource ) const;

	public:
		bool											ConfigureOptical						(  bool singlescatter, const SKTRAN_TableOpticalProperties_V21*	optProp, SKTRAN_ThreadManager* threadmanager  );
//		bool											InitializeFirstOrderIncomingRadiances	();
//		bool											ScatterIncomingRadiance					( bool ignorehighaltdiffuse );
//		bool											CalculateIncomingRadiance				( bool ignorehighaltdiffuse );
		size_t											NumDiffusePoints						() const;
		bool											LookupDiffusePoint						( size_t pointindex, SKTRAN_DiffusePointOptical_V21** point );
};



/*-----------------------------------------------------------------------------
 *					SKTRAN_GroundPointDiffuseGeometry_V21		2008-2-6*/
/** A class that represents a ground point. 
 **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_GroundPointDiffuseOptical : public SKTRANSO_JindexTableBase
{
	private:
		SKTRAN_JValueTable_V21							m_diffuseupwardflux;		//!< Table that evaluates the upward flux from the incoming radiances
		SKTRAN_StokesScalar								m_FupMS;					//!< The upward flux at this point 
		const SKTRAN_GroundPointDiffuseGeometry_V21*	m_geometrypoint;			//!< The geometry describing this ground point and the incoming radys.

	private:
		void											ReleaseResources				();


	public:
														SKTRANSO_GroundPointDiffuseOptical	(const SKTRAN_GroundPointDiffuseGeometry_V21* geometry);
		virtual										   ~SKTRANSO_GroundPointDiffuseOptical	();
		bool											AttachToGeometry					() ;
		bool											ConfigureOptical					( const SKTRAN_TableOpticalProperties_V21* optProp, SKTRANSO_Quadrature_TLS_V21* quadrature );
		bool											UpdateDiffuseUpwardFlux				();
		const SKTRAN_StokesScalar*						DiffuseUpwardFlux					() const { return &m_FupMS;}
		const SKTRAN_GroundPointDiffuseGeometry_V21*	Geometry							() const { return m_geometrypoint;}


		virtual bool									InterpolateTable					( const HELIODETIC_POINT&       location,
																							  const HELIODETIC_UNITVECTOR&	look,
																							  SKTRANSO_JIndex*				vertexdescriptortable, 
																							  size_t						maxpts, 
																							  size_t*						npts,
																							  double						weight ) const override { NXASSERT((false)); return false;}

		virtual const SKTRAN_StokesScalar*				ConvertJIndexToRadiancePtr			( const SKTRANSO_JIndex*		entry,
																							  ENUM_SKTRAN_JSOURCE			jsource
																							) const override;

};	


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21				2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableGroundPointDiffuseOptical_V21
{
	private:
		SKTRANSO_TableGroundPointDiffuse*		m_geometrytable;
		SKTRANSO_GroundPointDiffuseOptical**	m_groundpoints;
		size_t									m_numgroundpointsresvd;

	private:
		bool									AllocateStorage								( size_t numpoints );
		void									ReleaseResources							( );


	public:
												SKTRAN_TableGroundPointDiffuseOptical_V21	( SKTRANSO_TableGroundPointDiffuse* groundpointtable);
											   ~SKTRAN_TableGroundPointDiffuseOptical_V21	( );
		bool									AttachToGeometry							( );
		bool									ConfigureOptical_Stage2MT					( size_t pointindex, const SKTRAN_TableOpticalProperties_V21* optprop, SKTRANSO_Quadrature_TLS_V21* threadquadrature );
		const SKTRAN_StokesScalar*				LookupJValue_DiffuseUpwardFlux				( const SKTRANSO_JIndex* index ) const;
		const SKTRAN_StokesScalar*				LookupJValue_DiffuseIncomingRays			( const SKTRANSO_JIndex* index ) const;
		bool									ScatterIncomingRadiance						( );
		size_t									NumPoints									( ) const;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmissionOptical_V2		2007-11-20*/
/**	A class used to calculate the solar transmission over the range of
 *	altitudes and solar zenith angles encountered in the radiative transfer
 *	problem. This class manages the wavelength dependent aspects of the calculation
 *	and is resident inside one of the algorithms processing threads. It complements the 
 *	corresponding geometry class. 
 *
 *	\par See Also
 *	#SKTRANSO_TableSolarTransmission
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableSolarTransmissionOptical_V21
{
	private:
		SKTRANSO_TableSolarTransmission*			m_geometrytable;
		SKTRAN_StokesScalar* 						m_transmission;			// An array [m_numangles, m_numradii] of transmission from the Sun to the location of the shells [m_numangles, m_numradii];
		size_t										m_numpoints;

	private:
		void										ReleaseResources							();
		bool										AllocateMemory								( size_t numpoints );

	public:
													SKTRAN_TableSolarTransmissionOptical_V21	(SKTRANSO_TableSolarTransmission* geometrytable);												//!< Default constructor
		virtual									   ~SKTRAN_TableSolarTransmissionOptical_V21	();

	public:
		bool										AttachToGeometry							();
		const SKTRAN_StokesScalar*					RayAt										( size_t position ) const;
		bool										ConfigureOptical_Stage2						( size_t pointindex, SKTRANSO_Quadrature_TLS_V21* threadquadrature );

};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_TableLinesOfSightOptical_V21		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableLinesOfSightOptical_V21
{
	private:
		size_t											m_numrays;
		size_t											m_numraysreserved;
		SKTRAN_RayLOSOptical_V21*						m_rayoptical;
		SKTRAN_TableLinesOfSight_V21*					m_geometry; 


	private:
		bool											AllocateRays							( size_t n );

	public:
														SKTRAN_TableLinesOfSightOptical_V21		(SKTRAN_TableLinesOfSight_V21* geometry);
		virtual										   ~SKTRAN_TableLinesOfSightOptical_V21		();
		void											ReleaseResources						();
		size_t											NumRays									() const { return m_numrays;}
		SKTRAN_RayLOSOptical_V21*						RayAtVar								( size_t index) { NXASSERT((index < m_numrays)); return m_rayoptical+index;}
		bool											AttachToGeometry						( );
		bool											ConfigureOptical						( bool singlescatter, const SKTRAN_TableOpticalProperties_V21*  optprop, SKTRAN_ThreadManager* threadmanager );
		bool											CalculateObserverIncomingRadiance		( SKTRAN_StokesScalar* buffer, size_t maxrays);

};
