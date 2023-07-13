/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints		2007-12-14*/
/** A class used to hold the diffuse profile table. This is essentially a 2-D
 *	array of points distributed in altitude and solar zenith angle. The distribution
 *	of points is defined by the SKTRAN_GridDefDiffuseHeights_V21 and the SZARAnge
 *	and SZADelta parameters in method Configure.
 *
 *	Each point in the table defines a unit sphere and manages the incoming and
 *	outgoing distribution of diffuse radiation.  The distribution of incoming
 *	radiances are defined on the unit sphere by the SKTRAN_GridDefDiffuseIncomingAzimuth_V21 and
 *	SKTRAN_SpecificationsMultiScatter_V2 grids.  The outgoing radiances are defined
 *	by SKTRAN_UnitSphere_WithLookupTable_V2.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableDiffusePoints_2D_Height_SZA : public SKTRANSO_TableDiffusePoints
{
	private:
		const SKTRAN_GridDefDiffuseHeights_V21*				m_radii;


	private:
		size_t												FindingBoundingAltitudeIndices	( const HELIODETIC_POINT&	location,
																							  SKTRAN_GridIndex*			index,
																							  double*					weight,
																							  size_t					maxvertices
																							) const;


		bool												InterpolateTable				( const HELIODETIC_POINT&			location,
																							  const HELIODETIC_UNITVECTOR&		look,
																							  SKTRANSO_JIndex*					vertexdescriptortable,
																							  size_t							maxpts,
																							  size_t*							npts,
																							  double							weight
																							) const;


	public:
															SKTRAN_TableDiffusePoints_2D_Height_SZA	();
		virtual											   ~SKTRAN_TableDiffusePoints_2D_Height_SZA	();

		size_t												FindingBoundingLocationIndices	( const HELIODETIC_POINT&	location,
																						      SKTRAN_GridIndex*			index,
																						      double*					weight,
																						      size_t					maxvertices
																							) const;

		virtual size_t										MaxInterpolatedPointsPerJ		() const { return 12;}		// 4 spatial points by 3 vertices per point

		virtual bool										ConfigureGeometry				( const SKTRAN_SpecsInternal_V21*		specs,
																							  SKTRAN_ThreadManager*					threadmanager ) override;


};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_TableGroundPointDiffuse_Colocated		2008-2-6*/
/** The table of ground points are used to evalaute the upwelling radiation
 *	coming from the ground. They are used for the radiance calculation at the
 *	end point of a ray that hits the ground.
 *
 *	The current version of SASKTRAN demands that the ground point table have
 *	the same angular division as the diffuse profile table and that the diffuse
 *	profile table include points at the ground.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableGroundPointDiffuse_Colocated : public SKTRANSO_TableGroundPointDiffuse
{
	private:
		const SKTRAN_TableDiffusePoints_2D_Height_SZA*				m_diffusetable;

	private:
		void														ReleaseResources				();
		size_t														FindingBoundingLocationIndices	( const HELIODETIC_POINT& location, SKTRAN_GridIndex* index, double* weight, size_t maxvertices) const;
		const SKTRAN_TableDiffusePoints_2D_Height_SZA*				DiffusePointsTable				() const { return m_diffusetable;}

	public:
																	SKTRAN_TableGroundPointDiffuse_Colocated( const SKTRAN_TableDiffusePoints_2D_Height_SZA* diffusetable);
		virtual													   ~SKTRAN_TableGroundPointDiffuse_Colocated();

	public:


		virtual size_t												NumberOfGroundPointsRequired				( const SKTRAN_SpecsInternal_V21*	specs  ) const;

		virtual bool												ConfigureGeometry_Stage2MT					( size_t					pointindex,
																												  SKTRANSO_Quadrature_TLS_V21* threadquadrature
																												);

		virtual bool												InterpolateTable							( const HELIODETIC_POINT&		location,
																												  const HELIODETIC_UNITVECTOR&	look,
																												  SKTRANSO_JIndex*				vertexdescriptortable,
																												  size_t						maxpts,
																												  size_t*						npts,
																												  double						weight ) const override;

};



/*-----------------------------------------------------------------------------
 *					class SKTRAN_TableSolarTransmission_2D_Height_SZA		2007-11-16*/
/** @ingroup common
 *	An implementation of the Solar Transmission Table class,
 *	#SKTRANSO_TableSolarTransmission. This class is used
 *	when it is sufficient to only consider the variation of solar optical depth with
 *	solar zenith angle. It should not be used when atmospheric conditions vary with
 *	solar longituide.
 *
 *	The table has a fairly high startup cost as it has a fine spacing of vertical
 *	profiles (default is 0.5 degrees in SZA) and each profile typically has 100
 *	vertical points, each of which sends a ray to the sun.  The table is indexed
 *	by solar zenith angle and height.
 *
 *	\par See Also
 *	#SKTRAN_TableSolarTransmissionOptical_V2
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableSolarTransmission_2D_Height_SZA : public SKTRANSO_TableSolarTransmission
{
	private:
		const SKTRAN_GridDefCosSZA_V21*					m_szagrid;					//!< The cos(SZA) location of each profile
		const SKTRAN_GridDefRayTracingShells_V21*		m_altitudegrid;				//!< The altitude grid, which is assumed to be used by all heights

	private:
		void											ReleaseResources	();

	public:
														SKTRAN_TableSolarTransmission_2D_Height_SZA();
		virtual 									   ~SKTRAN_TableSolarTransmission_2D_Height_SZA();
		size_t											FindingBoundingLocationIndices	(  const HELIODETIC_POINT& location, SKTRAN_GridIndex* index, double* weight, size_t maxvertices) const;
		size_t											FindingBoundingAltitudeIndices	(  const HELIODETIC_POINT& location, SKTRAN_GridIndex* index, double* weight, size_t maxvertices) const;


	public:

		virtual bool											ConfigureGeometry			( const SKTRAN_SpecsInternal_V21* specs, SKTRAN_ThreadManager* threadmanager ) override;

		virtual bool											InterpolateTable			( const HELIODETIC_POINT&		location,
																						      const HELIODETIC_UNITVECTOR&	look,
																						      SKTRANSO_JIndex*				vertexdescriptortable,
																						      size_t						maxpts,
																						      size_t*						npts,
																						      double						weight ) const override;

//		virtual bool									ConvertJindexToPosition			( const SKTRANSO_JIndex* entry, size_t* positionindex ) const;
};


/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height		2007-11-20*/
/** A class used to cache the optical properties of the atmosphere. This class
 *	implements optical properties which are only a function of altitude.
 *	NOte the classes SKTRAN_QuadratureScatteringMatrixxxxxx are also affected by
 *	by this.
 **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_InternalEmissionPropertiesTable_1D_Height : public SKTRANSO_InternalEmissionPropertiesTable
{
	private:
		double													m_wavelen;
		double													m_mjd;
		bool													m_isempty;					//!< Flag used to indicate that there are no atmospheric emissions
		HELIODETIC_POINT										m_location;
		size_t													m_numheighttoindex;			//!< Number of elements in the height to index table (typically 10,000)
		size_t*													m_heighttoindextable;		//!< A table [m_numheighttoindex] that quickly converts a radius to an index with specified resolution
		double													m_heightindexresolution;	//!< The resolution o fthe radius to index rsolution, typically 10 meters
		double													m_minheight;
		const SKTRAN_GridDefOpticalPropertiesRadii_V21*			m_altitudegrid;				//!< The grid that defines the spatial grid for specifying optical properties
		std::vector<SKTRAN_StokesScalar>						m_radiance;					//!< Array [numshells] of isotropic radiance in atmosphere.
		SKTRAN_StokesScalar										m_groundradiance;			//!< The radiance at the ground surface.


	private:
		void										ReleaseResources			 ();
		void										ReleaseObjects				 ();
		bool										Allocate					 ( size_t numcells );
		bool										ConfigureAltitudeToIndexTable();
		bool										IndexOfPointBelowOrEqual	 ( double h0, size_t* lowindex ) const;
		bool										IndexOfPointEqualOrAbove	 ( double h0, size_t* hihindex ) const;
		size_t										NumShells					 () const { return (m_altitudegrid     == NULL)? 0 : m_altitudegrid->NumAltitudes();}


	public:
													SKTRANSO_InternalEmissionPropertiesTable_1D_Height		();
		virtual 								   ~SKTRANSO_InternalEmissionPropertiesTable_1D_Height		();
		virtual	bool								IsEmpty									() const override			{ return m_isempty;}
		virtual void								SetEmpty								(bool value) override       {m_isempty = value;}
		virtual bool								ConfigureGeometry						( const SKTRAN_SpecsInternal_Base* specs ) override;
		virtual bool								ConfigureOptical						( double wavelen, const SKTRAN_CoordinateTransform_V2* coords, SKTRAN_AtmosphericEmission* emissionstate ) override;
		virtual double								GetIsotropicRadianceInAtmosphere		( const HELIODETIC_POINT& point ) const override;
		virtual double								GetIsotropicGroundRadiance				( const HELIODETIC_POINT& point ) const override;

		virtual bool								InterpolateTable						( const HELIODETIC_POINT&			location,				// Location where interpolation required
																						      const HELIODETIC_UNITVECTOR&		look,					// Look direction where interpolation required
																								    SKTRANSO_JIndex*			vertexdescriptortable, 
																								    size_t						maxpts, 
																								    size_t*						npts,
																								    double						weight ) const override;

		virtual const SKTRAN_StokesScalar*			ConvertJIndexToRadiancePtr		( const SKTRANSO_JIndex*		entry, 
																					  ENUM_SKTRAN_JSOURCE			jsource ) const override;

};




/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableRayLOS		2010-6-21*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableRayLOS_Legacy : public SKTRANSO_TableRayLOS
{

	public:
																SKTRAN_TableRayLOS_Legacy   ( std::weak_ptr< const SKTRAN_RayFactory_Base> rayfactory);
		bool													ConfigureGeometry			( SKTRANSO_Quadrature_TLS_V21* quadrature,	 SKTRANSO_RayLOSGeometry_V21* losray );
};



class SKTRAN_TableRayLOSFactory_Legacy : public  SKTRAN_TableRayLOSFactory
{
	public:
										SKTRAN_TableRayLOSFactory_Legacy();
		virtual						   ~SKTRAN_TableRayLOSFactory_Legacy();
		virtual	bool 					CreateInternalSingleScatterTable( SKTRANSO_TableRayLOS** internalsinglescattertable,std::weak_ptr< const SKTRAN_RayFactory_Base> rayfactory  ) const override;
};

