
/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineInterface_V2		2008-3-6*/
/** **/
/*---------------------------------------------------------------------------*/




/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffusePointOptical_V21		2007-12-17*/
/** This class represents a single point on a single profile used in the table
 *	of diffuse profiles.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_DiffusePointOptical_V21 : public nxUnknown
{
	private:
		const SKTRAN_DiffusePointGeometry_V21*				m_geometry;					//!< The difuuse geometry table, this has the grid for the outgoing unit sphere as well as the incoming azimuth and zenith grid definitions
		size_t												m_numlosIn;					//!< The number of lines of sight coming into this point
		std::vector<SKTRAN_RayInternalDiffuseOptical_V2>	m_losIn;					//!< The array [m_numlosIn] of lines of sight coming into this point
		std::vector<SKTRAN_StokesScalar>					m_incomingradiance;			//!< The array [m_numlosIn] of incoming radiance, generated from integrating along 
//		SKTRAN_StokesScalar*								m_outdiffuseforwardJ;		//!< The array [m_numlosIn] of outgoing source functions due to forward delta function scatter at this point, may be NULL
		std::vector<SKTRAN_StokesScalar>					m_outgoingJ;				//!< The array [m_numoutgoingrad] of outgoing source functions at this point
		size_t												m_numoutgoingJ;				//!< The number of outgoing source functions at this point. This is set by the unit sphere
		SKTRAN_ScatteringMatrixPointOptical_V21				m_scatterobject;			//!< Object to scatter the incoming rays to the outbound directions

	private:
		void										ReleaseLos	();
		bool										AllocateLos	();
		void										ReleaseOutgoingJ				();
		bool										AllocateOutgoingJ				( const SKTRAN_UnitSphere_V2* unitsphere);

	public:
													SKTRAN_DiffusePointOptical_V21	( const SKTRAN_DiffusePointGeometry_V21* geometry );
		virtual									   ~SKTRAN_DiffusePointOptical_V21	();


		bool										ConfigureOptical				( const SKTRAN_TableOpticalProperties_V21*	optprop,
																				      SKTRANSO_Quadrature_TLS_V21*                 quadrature );

		bool										AttachToGeometry				( );
//		bool										SetScatteringMatrix				( SKTRAN_ScatteringMatrixPointOptical_V21* scatterobject );

		SKTRAN_RayInternalDiffuseOptical_V2*		IncomingLineOfSight				( size_t idx )						{ return (&m_losIn.at(idx));}
		const SKTRAN_RayInternalDiffuseOptical_V2*	IncomingLineOfSight				( size_t idx )				const	{ return (&m_losIn.at(idx));}
		const SKTRAN_StokesScalar*					IncomingRadianceArray			()							const	{ return (&m_incomingradiance.front());}
		const SKTRAN_StokesScalar*					OutgoingJPtr					( size_t unitsphereindex )	const	{ return (&m_outgoingJ.at(unitsphereindex));}
		SKTRAN_StokesScalar*						OutgoingJPtr					( size_t unitsphereindex )			{ return (&m_outgoingJ.at(unitsphereindex));}
		SKTRAN_ScatteringMatrixPointOptical_V21*		ScatteringMatrixObject			()								{ return &m_scatterobject;}	
		size_t										NumIncomingLos					() const							{ return m_numlosIn; }
		size_t										NumOutgoingRad					() const							{ return m_numoutgoingJ; }
		const SKTRAN_DiffusePointGeometry_V21*		Geometry						() const							{ return m_geometry;}
		void										ClearOutgoingJ					();
		bool										CalculateIncomingRadiances		( bool ignorehighaltdiffuse );
		bool										ScatterIncomingRadiance			( bool ignorehighaltdiffuse );
		bool										InitializeFirstOrderIncomingRadiances();
};

/*-----------------------------------------------------------------------------
 *					class DiffusePointGeometry						2007-11-7*/
/** A single point in a diffuse profile. The point is at a given heliodetic
 *	location and stores all of the geometry information. The point stores 
 *	an array of incoming rays and an array of outbound source functions.
 *	It also has a copy of the scattering matrix details.
 *
 *	\par See Also
 *	#SKTRAN_DiffuseProfileGeometry_V21
 *	#SKTRANSO_TableDiffusePoints
 *	
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_DiffusePointGeometry_V21
{

	private:
		HELIODETIC_POINT														m_location;					//!< The location of this point in heliodetic coordinates
		SKTRAN_GridIndex														m_diffuseheightindex;		//!< The height index of this point in the parent diffuse points table
		SKTRAN_GridIndex														m_diffuseprofileindex;		//!< The profile index of this point in the parent diffuse points table.
		std::shared_ptr< const SKTRAN_RayFactory_Base>							m_diffuserayfactory;		//!< The ray factory used to create the diffuse rays.
		std::vector< SKTRANSO_RayInternalGeometry >								m_losGeometryIn;			//!< The array of incoming line of sights (numazimuths-1 * numzeniths-1)
		const SKTRAN_UnitSphereLatLonGrid*										m_incomingunitvectors;		//!< The grid defining the incoming unit vectors (away from this point) 

		SKTRAN_ScatteringMatrixPointGeometry_V21			m_scatteringmatrix;
		SKTRAN_DiffusePointOptical_V21*						m_opticalpoint;					//!< The optical point associated with this point.
		bool												m_ishighaltitudediffuse;		//!< An option used to turn off high altitude diffuse processing

	private:
		void												ReleaseResources					();
		bool												AllocateLOSArray					();
		bool												ConfigureLOSArray					( SKTRANSO_Quadrature_TLS_V21* quadraturespecs );

	public:
															SKTRAN_DiffusePointGeometry_V21		();
		virtual											   ~SKTRAN_DiffusePointGeometry_V21		();
//															SKTRAN_DiffusePointGeometry_V21		( SKTRAN_DiffusePointGeometry_V21& moveother);
//		void												operator =							( SKTRAN_DiffusePointGeometry_V21&& moveother);

		// ----- Multi-threaded functions
		bool												ConfigureGeometry_Stage2MT			( SKTRANSO_Quadrature_TLS_V21* quadrature );
		bool												CreateJIndexTables_MT				( SKTRANSO_Quadrature_TLS_V21* quadrature );

		// ----- Single threaded functions


		bool												ConfigureGeometry_Stage1			( std::shared_ptr<const SKTRAN_RayFactory_Base>&  diffuserayfactory,
																							      //const SKTRAN_SpecsInternal_RayTracing_V21*	  raytracingspecs,
																								  const SKTRAN_SpecsInternal_Diffuse_V21*		  diffusespecs,
																								  const HELIODETIC_POINT&						  point,
																								  SKTRAN_GridIndex								  diffuseprofileindex,
																								  SKTRAN_GridIndex								  diffuseradiusindex
																								);

		bool												SetScatteringMatrix					( SKTRAN_ScatteringMatrixPointGeometry_V21* scattermatrix );

//		bool												ZenithAndAzimuthOfIncomingLos		( size_t idxoflosin, double* zenith, double*azimuth )	const; 
//		bool												ZenithAndAzimuthOfOutboundSource	( size_t idxofsrc,   double* zenith, double*azimuth )	const;
		bool												CreateOpticalPoint					();
		size_t												NumLOSIn							()				const { return m_losGeometryIn.size();}
		const SKTRANSO_RayInternalGeometry&					LOSAt								(size_t idx)	const { return m_losGeometryIn.at(idx);} 
		SKTRANSO_RayInternalGeometry*						LOSAtVar							(size_t idx)		  { return &m_losGeometryIn.at(idx);} 
		const HELIODETIC_POINT&								Location							()				const { return m_location;}
		HELIODETIC_UNITVECTOR								LocalZenith							()				const;

//		size_t												NumIncomingZenith					()				const { return m_zenithsIn->NumAngles();}
//		size_t												NumIncomingAzimuth					()				const { return m_azimuthsIn->NumAngles();}
		size_t												NumIncomingRays						()				const;
		const SKTRAN_UnitSphereLatLonGrid*					IncomingUnitVectors					()				const { return m_incomingunitvectors;}

//		bool												IncomingSolidAngleArray				( const double** domega, size_t* numangles )const;
		const SKTRAN_ScatteringMatrixPointGeometry_V21*		ScatteringMatrix					() const { return &m_scatteringmatrix;}
		SKTRAN_DiffusePointOptical_V21*						OpticalPoint						()       { return m_opticalpoint;}
		const SKTRAN_DiffusePointOptical_V21*				OpticalPoint						() const { return m_opticalpoint;}
		size_t												ProfileIndex						() const { return m_diffuseprofileindex;}
		bool												IsHighAltitude						() const { return m_ishighaltitudediffuse;}
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileGeometry_V21		2007-12-12*/
/** An array of diffuse points that line at different altitudes at a fixed
 *	solar zenith angle.  Each profile is an array of points at the same solar
 *	zenith angle.  The lowest point is always the ground and must equal the
 *	
 **/
/*---------------------------------------------------------------------------*/


class SKTRAN_DiffuseProfileGeometry_V21
{
	private:
		const SKTRANSO_TableDiffusePoints*			m_parent;						//!< The parent Table that owns this diffuse profile.
		std::vector<SKTRAN_DiffusePointGeometry_V21>	m_diffusepoints;				//!< The array [m_numpoints] of diffuse points 
//		size_t											m_numpoints;					//!< The number of points excluding the ground point
		HELIODETIC_POINT								m_location;						//!< The location of this profile
		SKTRAN_GridIndex								m_profileidx;					//!< The index of the profile in the array of profiles, provided by the Table Of governing table.

	private:
		void											ReleaseResources				();
		bool											AllocatePoints					(size_t numradii);


	public:
														SKTRAN_DiffuseProfileGeometry_V21();
		virtual										   ~SKTRAN_DiffuseProfileGeometry_V21();
		SKTRAN_GridIndex								ProfileIndex					()				const { return m_profileidx;}	
		const HELIODETIC_POINT&							Location						()				const { return m_location;}
		size_t											NumPoints						()				const { return m_diffusepoints.size();}
		const SKTRAN_DiffusePointGeometry_V21*			At								(size_t idx)	const { return &(m_diffusepoints.at(idx));}
		SKTRAN_DiffusePointGeometry_V21*				AtVar							(size_t idx)		  { return &(m_diffusepoints.at(idx));}

		const SKTRANSO_TableDiffusePoints*			Parent							()				const { return m_parent;}

		bool											ConfigureGeometry_Stage1		( const  HELIODETIC_POINT&							 location,
																						  SKTRAN_GridIndex							         profileidx,
																						  std::shared_ptr<const SKTRAN_CoordinateTransform_V2> 	coords, 
																						  std::shared_ptr<const SKTRAN_RayFactory_Base>			diffuserayfactory,
																						  const  SKTRAN_SpecsInternal_Diffuse_V21*		     diffusespecs,
																						  const  SKTRANSO_TableDiffusePoints*		 parent
																						);

		bool											CreateOpticalProfile			( SKTRAN_DiffuseProfileOptical_V21** optprofile );
};



/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileOptical_V21		2007-12-14*/
/** This is a profile of diffuse points at a given solar zenith angle.
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_DiffuseProfileOptical_V21
{
	private:
		SKTRAN_DiffuseProfileGeometry_V21*				m_geometry;			//!< The geometry (wavelength independent aspects) of this vertical profile
		std::vector<SKTRAN_DiffusePointOptical_V21*>		m_points;			//!< The points storing the wavelength dependent parts of this profile

	private:
		void											ReleaseResources();

	public:
														SKTRAN_DiffuseProfileOptical_V21	(SKTRAN_DiffuseProfileGeometry_V21* geometry);
													   ~SKTRAN_DiffuseProfileOptical_V21	();
		bool											AttachToGeometry				();
		size_t											NumPoints						() const { return m_points.size(); }
		const SKTRAN_DiffusePointOptical_V21*			At								(size_t idx)	const { return m_points[idx];}
		SKTRAN_DiffusePointOptical_V21*					AtVar							(size_t idx)		  { return m_points[idx];}
//		bool											ScatterIncomingRadiance			(bool ignorehighaltdiffuse);
//		bool											CalculateIncomingRadiance		(bool ignorehighaltdiffuse);
//		bool											InitializeFirstOrderIncomingRadiances();
};




/*-----------------------------------------------------------------------------
 *					SKTRAN_GroundPointDiffuseGeometry_V21		2008-2-6*/
/** A class that represents the geometry of a ground point. The ground
 *	points require special treatment as they typically scatter incoming
 *	radiation using a lambertian or BRDF.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_GroundPointDiffuseGeometry_V21
{
	private:
		SKTRANSO_JIndexArray							m_diffusedownwardfluxtable;		//!< Table the evaluates the downward flux from the incoming radiances
		SKTRAN_GridIndex								m_pointindex;					//!< The index of this point in the Ground Points Table
		const SKTRAN_DiffusePointGeometry_V21*			m_diffusepoint;					//!< The diffuse point that is coincident with this point. They are both at the lowest shell
		SKTRANSO_GroundPointDiffuseOptical*			m_opticalpoint;					//!< The (one and only) optical point that corresponds to this geometry point

	private:
		void											ReleaseResources						();

	public:
														SKTRAN_GroundPointDiffuseGeometry_V21	();
		virtual										   ~SKTRAN_GroundPointDiffuseGeometry_V21	();
		HELIODETIC_UNITVECTOR							LocalZenith								() const;
		bool											ConfigureGeometry						(  SKTRAN_GridIndex pointindex, const SKTRAN_DiffusePointGeometry_V21* diffusepoint);
		const SKTRANSO_JIndexArray*					DownwardFluxJIndexTable					() const { return &m_diffusedownwardfluxtable;}
		const SKTRAN_DiffusePointGeometry_V21*			DiffuseGeometryPoint					() const { return m_diffusepoint;}
		      SKTRANSO_GroundPointDiffuseOptical*		OpticalPoint							()       { return m_opticalpoint;}
		const SKTRANSO_GroundPointDiffuseOptical*		OpticalPoint							() const { return m_opticalpoint;}					//!< The (one and only) optical point that corresponds to this geometry point
		bool											CreateIncomingJIndexTableForDownwardFlux();       
};



class SKTRAN_TableLinesOfSight_V21;
class SKTRAN_TableLinesOfSightOptical_V21;



/*-----------------------------------------------------------------------------
 *					class SKTRAN_TableLinesOfSight_V21		2008-4-21*/
/** \internal
 *	Stores the lines of sight used in the model. The class was modified in
 *	June 2013 to adjust the rays so they reported the same tangent altitude
 *	in the osculating sphere system as the they did in the real earth geoid.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableLinesOfSight_V21
{
	private:
		std::shared_ptr< const SKTRAN_RayFactory_Base>					m_rayfactory;
		std::shared_ptr< const SKTRAN_RayFactory_Base>					m_solarrayfactory;
		std::vector< SKTRANSO_RayLOSGeometry_V21 >						 m_raygeometry;
		SKTRAN_TableLinesOfSightOptical_V21*							m_opticaltable;
		SKTRAN_LineOfSightArray_V21										m_observerlinesofsight;

	private:
		bool										AllocateRays									( size_t n, const SKTRAN_TableRayLOSFactory* singlescatterfactory );
		bool										TranslateLimbViewingGeometryToOsculatingSphere	( SKTRAN_LineOfSightEntry_V2*	entry, const SKTRAN_CoordinateTransform_V2* coordinates);
	public:
													SKTRAN_TableLinesOfSight_V21		();
													SKTRAN_TableLinesOfSight_V21		(SKTRAN_TableLinesOfSight_V21&& moveother);
		virtual									   ~SKTRAN_TableLinesOfSight_V21		();
		void										ReleaseResources					();

		bool										SetLinesOfSight						( const SKTRAN_LineOfSightArray_V21&	 observerlinesofsight,
																						  const SKTRAN_CoordinateTransform_V2*	coordinates,
																						  SKTRAN_RayTracingRegionManager*	    raytracingregionmanager);

		bool										ConfigureGeometry					( bool								 singlescatter, 
																						  const SKTRAN_SpecsInternal_V21*	 modelspecs,
																						  SKTRAN_ThreadManager*				 threadmanager
																						);

		bool										ConfigureOptical					( bool										singlescatter, 
																						  const SKTRAN_TableOpticalProperties_V21*   optprop,
																						  SKTRAN_ThreadManager*						threadmanager );

		size_t												NumRays								()				const		{ return m_raygeometry.size();}
		const SKTRANSO_RayLOSGeometry_V21*					RayAt								(size_t idx)	const		{ return &m_raygeometry.at(idx);}
		SKTRANSO_RayLOSGeometry_V21*						RayAtVar							(size_t idx)				{ return &m_raygeometry.at(idx);}
//		const SKTRAN_TableLinesOfSightOptical_V21*			OpticalTable						() const					{ return m_opticaltable;}
		SKTRAN_TableLinesOfSightOptical_V21*				OpticalTableVar						()							{ return m_opticaltable;}
		SKTRAN_LineOfSightArray_V21*						LinesOfSight						() 							{ return &m_observerlinesofsight;}
		const SKTRAN_LineOfSightArray_V21*					LinesOfSightConst					() const 					{ return &m_observerlinesofsight;}
		const SKTRAN_RayFactory_Base*						RayFactory					()									{ return m_rayfactory.get();}

};




/*-----------------------------------------------------------------------------
 *					class SKTRAN_ThreadManager		2010-3-22*/
/** This a relatively small class that manages the muti-threading activity
 *	within Sasktran. This class is closely coupled with class SKTRAN_ThreadInstanceBoost.
 *	This class is the "manager" which is called from the main processing thread
 *	running in Sasktran. The class notifies all of its "child " threads that there
 *	is some processing to be done and fires them off. The manager class waits for
 *	all of the child thread classes to finish their work before exiting.
 *
 *	This class represents a significant increase in complexity over the
 *	original Sasktran which only multithreaded in wavelength. This class now
 *	multi-threads through several areas (12 at the last count) of the code. I 
 *	did not want to waste CPU cycles creating and destroying threads so I
 *	create all of the childs threads just once and then put them to sleep
 *	waiting for a notification from the parent. 
 *
 *	The multi-threading is based upon each child thread processing one point or
 *	one ray depending upon context. This keeps the syncronization of objects fairly straight
 *	forward. (Multiple threads working on one ray would require a
 *	lot of syncronization objects and associated complexity). The child thread
 *	asks the the manage for the next point or ray to process. The point or ray
 *	is retrieved from the appropriate Sasktran table and processed. Each threads loops around
 *	until the manager decrees that all the poiints are processed.
 *
 *	The manager class has a specific function for each of the (12 at the
 *	last count) areas of Sasktran that multi-threaded.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_ThreadManager : public nxUnknown
{
	private:

		// These variables are just temporary variables that are accessible from the worker threads.
		// I prefer to use this method as it lets me avoid do a whole bunch of typecasting thus increasing
		// overall reliability.

		SKTRANSO_TableDiffusePoints*					m_diffusetable_geometry;				
		SKTRAN_TableDiffusePointsOptical_V21*			m_diffusetable_optical;
		SKTRANSO_TableSolarTransmission*				m_solartable_geometry;
		SKTRANSO_TableGroundPointDiffuse*				m_groundpointtable_geometry;
		SKTRAN_TableLinesOfSight_V21*					m_lostable_geometry;
		SKTRAN_TableLinesOfSightOptical_V21*			m_lostable_optical;
		const SKTRAN_TableOpticalProperties_V21*		m_optprop;
		const SKTRAN_LineOfSightArray_V21*				m_observerspecs;
		const SKTRAN_SpecsInternal_V21*					m_specifications;
		bool											m_singlescatter;
		bool											m_ignorehighalt;

	private:
		size_t										m_numpoints;				// The number of points or rays to process by the thread array

	private:
		bool										(SKTRAN_ThreadManager::*	m_threadfunction)  (SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);

	private:
		bool										TLS_DiffusePointsTable_ConfigureGeometryStage2		 ( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);
		bool										TLS_DiffusePointsTable_CreateJIndices				 ( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);
		bool										TLS_DiffusePointsTable_ConfigureOpticalStage2		 ( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);
		bool										TLS_DiffusePointsTable_EvaluateIncomingRays		     ( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex );
		bool										TLS_DiffusePointsTable_ScatterIncomingRadiance		 ( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex );
		bool										TLS_DiffusePointsTable_FirstOrderInitialize          ( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex );

		bool										TLS_SolarTransmissionTable_ConfigureGeometryStage2	( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);
		bool										TLS_SolarTransmissionTable_ConfigureOpticalStage2	( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);

		bool										TLS_GroundPointsTable_ConfigureGeometryStage2		( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);
		bool										TLS_GroundPointsTable_CreateJIndexTables			( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);
		bool										TLS_GroundPointsTable_ConfigureOpticalStage2		( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);
		bool										TLS_GroundPointsTable_ScatterIncomingRadiance		( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);

		bool										TLS_LinesOfSightTable_ConfigureGeometryStage2		( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);
		bool										TLS_LinesOfSightTable_ConfigureOpticalStage2		( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex);

	protected:
		size_t										NumPoints							() const				{ return m_numpoints;} 
		void										ResetNumPoints						( size_t numpoints)		{ m_numpoints = numpoints;}



	public:
		/** Called by each child worker thread to execute the action requested. The instance of SKTRAN_ThreadManager
		 *	will save the requested action before starting up the child threads with a call to #NotifyWorkerThreadsAndWaitForEnd.
		 *	NotifyWorkerThreads results in each worker thread calling this function which executes in the context of the
		 *	worker thread, provided via #threadquadrature. The derived instance of SKTRAN_ThreadManager will use this
		 *	function to "switch" to the desired thread safe task at hand. The function will typically loop 
		 *	through each point in a table (in tandem with other worker threads) and process that point. 
	     *
		 *	\param threadquadrature
		 *	The quadrature object designated to this worker thread. The quadrature object usually holds all the
		 *	necessary thread specific buffers and context necessary to process the action.
		 *
		 *	\return
		 *	true if success otherwise false.
		**/
		
		bool								ThreadExecuteAction								( SKTRANSO_Quadrature_TLS_V21*	threadquadrature, size_t pointindex );

		/** Called by the main processing to complete stage 2 processing of the geometry diffuse points table.
		 *	The stage 2 processing uses several worker threads to plow through the processing. Returns true if the 
		 *  processing is successful, error otherwise.
		 *
		 *	\param table
		 *	The table being processed. This table will be modified upon exit.
		**/

		bool								DiffusePointsTable_ConfigureGeometryStage2		( SKTRANSO_TableDiffusePoints*	table  );

		/** Called by the main processing to create the JindexTables for each ray of each point in the geometry diffuse points table.
		 *	The stage 2 processing uses several worker threads to plow through the processing. Returns true if the 
		 *  processing is successful, error otherwise.
		 *
		 *	\param table
		 *	The table being processed. This table will be modified upon exit.
		**/

		bool								DiffusePointsTable_CreateJIndices				( SKTRANSO_TableDiffusePoints*	table  );

		/** Called by the main processing to complete stage 2 processing of the optical diffuse points table.
		 *	The stage 2 processing uses several worker threads to plow through the processing. Returns true if the 
		 *  processing is successful, error otherwise.
		 *
		 *	\param table
		 *	The table being processed. This table will be modified upon exit.
		**/

		bool								DiffusePointsTable_ConfigureOpticalStage2		( SKTRAN_TableDiffusePointsOptical_V21*			table,
																							  const SKTRAN_TableOpticalProperties_V21*		optprop  );

		bool								DiffusePointsTable_EvaluateIncomingRays          ( SKTRAN_TableDiffusePointsOptical_V21*			table,
																							   bool ignorehighalt
																							  );

		bool								DiffusePointsTable_ScatterIncomingRadiance       ( SKTRAN_TableDiffusePointsOptical_V21*			table,
																							   bool ignorehighalt
																							  );

		bool								DiffusePointsTable_FirstOrderInitialize          ( SKTRAN_TableDiffusePointsOptical_V21*			table);


		/** Called by the main processing to complete stage 2 processing of the Solar Transmission table. This
		 *	usually calculates the extinction from a point in the atmosphere to the sun. Each point in the
		 *	table is handled by several worker threads. Returns true if the processing is successful, error otherwise.
		 *
		**/

		bool								SolarTransmissionTable_ConfigureGeometryStage2	( SKTRANSO_TableSolarTransmission*			solartable,
																							  const SKTRAN_SpecsInternal_V21*				specs );

		bool								SolarTransmissionTable_ConfigureOpticalStage2	( SKTRANSO_TableSolarTransmission*			solartable,
																							  const SKTRAN_TableOpticalProperties_V21*		optprop );


		bool								GroundPointsTable_ConfigureGeometryStage2		( SKTRANSO_TableGroundPointDiffuse*			groundpointstable);

		bool								GroundPointsTable_CreateJIndexTables			( SKTRANSO_TableGroundPointDiffuse*			groundpointstable);

		bool								GroundPointsTable_ConfigureOpticalStage2		( SKTRANSO_TableGroundPointDiffuse*			groundpointstable,
																							  const SKTRAN_TableOpticalProperties_V21*		optprop);

		bool								GroundPointsTable_ScatterIncomingRadiance		( SKTRANSO_TableGroundPointDiffuse*			groundpointstable );



		bool								LinesOfSightTable_ConfigureGeometryStage2		( SKTRAN_TableLinesOfSight_V21*					lostable,
																							  const SKTRAN_LineOfSightArray_V21*				observerspecs,
																							  const SKTRAN_SpecsInternal_V21*				specs,
																							  bool											singlescatter
																							);

		bool								LinesOfSightTable_ConfigureOpticalStage2		( SKTRAN_TableLinesOfSightOptical_V21*			lostable,
																							  bool											singlescatter
																						    );


	public:
													SKTRAN_ThreadManager	();
		virtual									   ~SKTRAN_ThreadManager	();

	public:
		virtual bool								SetOpticalProps						( const SKTRAN_TableOpticalProperties_V21* optprop) = 0;
		virtual bool								CreateThreads						( size_t numthreads, const SKTRAN_SpecsInternal_V21* modelspecifications, const SKTRAN_EngineDiffuseTables* modeltables ) = 0;
		virtual bool								NotifyWorkerThreadsAndWaitForEnd	( size_t numpoints )  = 0;
		virtual bool								CloseThreadsAndRelease				( ) = 0;
};

extern SKTRAN_ThreadManager* SKTRAN_ThreadManagerBoost_CreateNewInstance();
extern SKTRAN_ThreadManager* SKTRAN_ThreadManagerOpenMP_CreateNewInstance();





/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineDiffuseTables		2010-4-5*/
/** This is a class that stores all of the tables used for diffuse
 *	calculations. These tables are never referenced if we are doing single
 *	scatter calculations.
 **/
/*---------------------------------------------------------------------------*/

 class SKTRAN_EngineDiffuseTables
{
	private:
//		bool										m_geometryisdirty;
		SKTRANSO_TableDiffusePoints*				m_diffusetable;				
		SKTRANSO_TableSolarTransmission*			m_solartransmissiontable;	
		SKTRANSO_TableGroundPointDiffuse*			m_diffusegroundpointtable;
		SKTRANSO_InternalEmissionPropertiesTable*	m_emissiontable;


	private:
//		bool										GeometryIsUnchanged				() const { return !m_geometryisdirty;}

	public:
		bool										DiffuseTablesAlreadyExist			();
//		void										SetDirtyGeometry					() { m_geometryisdirty = true;}
		void										ReleaseResources					();
		
		bool										CreateEmptyDiffuseTables			(  const SKTRAN_SpecsInternal_V21*		modelspecifications );

		
		bool										ConfigureDiffuseGeometryTables	    ( bool									singlescatter,
																						  SKTRAN_SpecsInternal_V21*				modelspecifications,
																						  SKTRAN_ThreadManager*					threadmanager );

		bool										ConfigureOpticalTables				( bool										singlescatter,
																						  const SKTRAN_TableOpticalProperties_V21*	opticalproptable,
																						  SKTRAN_ThreadManager*						threadmanager );

		bool										ConfigureOpticalEmissionTables		(  double									wavelen,
																						   const SKTRAN_CoordinateTransform_V2*		coords,
																						   SKTRAN_AtmosphericEmission*				atmosphericemissions);


	public:
															SKTRAN_EngineDiffuseTables	();
														   ~SKTRAN_EngineDiffuseTables	();
		const SKTRANSO_TableDiffusePoints*					DiffusePointsTable			() const	{ return m_diffusetable;}
		      SKTRANSO_TableDiffusePoints*					DiffusePointsTableVar		()			{ return m_diffusetable;}
		const SKTRANSO_TableSolarTransmission*				SolarTransmissionTable		() const 	{ return m_solartransmissiontable;}	
		      SKTRANSO_TableSolarTransmission*				SolarTransmissionTableVar	()			{ return m_solartransmissiontable;}	
		const SKTRANSO_TableGroundPointDiffuse*				DiffuseGroundPointTable		() const	{ return m_diffusegroundpointtable;}
		      SKTRANSO_TableGroundPointDiffuse*				DiffuseGroundPointTableVar	()			{ return m_diffusegroundpointtable;}
		const SKTRANSO_InternalEmissionPropertiesTable*		EmissionsTable				() const	{ return m_emissiontable;}
		      SKTRANSO_InternalEmissionPropertiesTable*		EmissionsTableVar			()			{ return m_emissiontable;}
 };


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine		2007-12-14*/
/** A top level class used to hold all of the various, geometry based tables
 *	used in the SASKTRAN model. This is a class seen by the user and acts
 *	as a interface between the user and the radiative transfer calculations.
 *	The class has one function that does work, #ConfigureGeometry. Its sole
 *	purpose is to calculate the required model geometry using the input
 *	specifications. The model geometry is used as a template by all of the
 *	optical threads which are actually calculating radiances.
 *
 *	The engine is passed to most (if not all) all of the geometry based objects created
 *	by this object.  This allows one object in the model to refer to
 *	another object in the model. Thus the engine is similiar to a container
 *	for global objects.
 *
 **/
/*---------------------------------------------------------------------------*/
/**
class SKTRANSO_Engine
{
	private:
		SKTRAN_ThreadManager*						m_threadmanager;
		SKTRAN_EngineDiffuseTables					m_tables;
		SKTRAN_TableLinesOfSight_V21				m_linesofsight;
		SKTRAN_SpecsInternal_V21					m_modelspecifications;
		SKTRAN_TableOpticalProperties_V21*			m_opticalpropertiestable;
		bool										m_updateclimatologycache;

	private:
		void										ReleaseResources								();
		bool										CalculateDiffuseRadiance						();
		bool										CreateOpticalPropertyTables						();
		bool										ConfigureOpticalTables							( bool singlescatter, const SKTRAN_TableOpticalProperties_V21*  opticalproperties);
		bool										CheckGeometryTables								( bool singlescatter );
		bool										SetLinesOfSight									( const SKTRAN_LineOfSightArray_V21*		  los );
		bool										InitializeFirstOrderIncomingRadiancesAndScatter	();	
		bool										CalculateIncomingRadianceAndScatter				( bool ignorehighaltdiffuse);
		bool										ScatterIncomingRadianceAtGround					( bool ignorehighaltdiffuse);

		bool										CalculateOpticalPropertiesTable_SingleThreaded	( double							  wavelen,
																									  SKTRAN_AtmosphericOpticalState_V21*  opticalstate,
																									  bool								  userupdateclimatology );
	public:
		bool										GetRayTracingGrid						( std::vector<double>* shells ) const;
		GEODETIC_INSTANT							ReferencePoint							() const;

	public:
													SKTRANSO_Engine						();
		virtual									   ~SKTRANSO_Engine						();

		bool										ConfigureModel							(       SKTRANSO_SpecificationsUser&	      modelspecifications,
																							  const SKTRAN_LineOfSightArray_V21&	  linesofsight,
																							  size_t                              numthreads );


		bool										CalculateRadiance						( std::vector<SKTRAN_StokesScalar>*		losradiance,
																							  double								wavelen,
																							  size_t								numordersofscatter,
																							  SKTRAN_AtmosphericOpticalState_V21*	opticalstate,
																							  bool									updateclimatology = false,
																							  SKTRAN_DiagnosticInterface*			diag = NULL);
};
*/



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine		2012-8-20*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_Engine : public SKTRAN_Engine_Base
{
	private:
		SKTRAN_AtmosphericOpticalState_V21*			m_lastopticalstateobject;
		SKTRAN_AtmosphericEmission*					m_lastatmosphericemissionobject;
		bool										m_isfirsttimeafterinit;

		SKTRAN_EngineDiffuseTables					m_tables;
		SKTRAN_TableLinesOfSight_V21				m_linesofsight;
		SKTRAN_SpecsInternal_V21					m_modelspecifications;
		SKTRAN_TableOpticalProperties_V21*			m_opticalpropertiestable;
		SKTRAN_ThreadManager*						m_threadmanager;
		bool										m_updateclimatologycache;
		
	private:
		bool										GetRayTracingGrid								( std::vector<double>* shells ) const;
		bool										CalculateOpticalPropertiesTable_SingleThreaded	( double							  wavelen,
																									  SKTRAN_AtmosphericOpticalState_V21* opticalstate,
																									  bool								  userupdateclimatology );
		bool										CalculateDiffuseRadiance						();
		bool										SetLinesOfSight									( const SKTRAN_LineOfSightArray_V21*  los );
		void										ReleaseResources								();
		bool										CreateOpticalPropertyTables						();
		bool										CalculateIncomingRadianceAndScatter				( bool ignorehighaltdiffuse);
		bool										ScatterIncomingRadianceAtGround					( bool ignorehighaltdiffuse);
		bool										InitializeFirstOrderIncomingRadiancesAndScatter	();	
		bool										ConfigureOpticalTables							( bool singlescatter, const SKTRAN_TableOpticalProperties_V21*  opticalproperties );
		bool										CheckGeometryTables								( bool singlescatter );
		bool										ConfigureOpticalEmissionTables					( double wavelen, SKTRAN_AtmosphericEmission*	atmosphericemissions);
		bool										CalculateOpticalPropertiesTable					( double										wavelen,
																									  SKTRAN_AtmosphericOpticalState_V21*			opticalstate,
																									  bool											userupdateclimatology
																									);

	public:
													SKTRANSO_Engine						();
		virtual									   ~SKTRANSO_Engine						();
		const SKTRAN_EngineDiffuseTables*			DiffuseTables						() const	{ return &m_tables;}
		const SKTRAN_TableLinesOfSight_V21*			LinesOfSight						() const	{ return &m_linesofsight;}
		const SKTRAN_SpecsInternal_V21*				InternalSpecs						() const	{ return &m_modelspecifications;}
		const SKTRAN_TableOpticalProperties_V21*	OpticalProperties					() const	{ return m_opticalpropertiestable;}
		GEODETIC_INSTANT							ReferencePoint						() const;

	public:

		virtual bool								ConfigureModel						(       SKTRAN_SpecsUser_Base&	      modelspecifications,
																							  const SKTRAN_LineOfSightArray_V21&  linesofsight,
																							  size_t                              numthreads ) override;

		virtual bool								CalculateRadiance						( std::vector<SKTRAN_StokesScalar>*		losradiance,
																							  double								wavelen,
																							  size_t								numordersofscatter,
																							  SKTRAN_AtmosphericOpticalState_V21*	opticalstate,
																							  std::vector<skRTStokesVector>*        losvector = nullptr,
																							  bool									updateclimatology = false,
																							  SKTRAN_DiagnosticInterface*			diag = NULL) override;
};



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine		2012-8-20*/
/** **/
/*---------------------------------------------------------------------------*/

// class SKTRANSO_Engine : public SKTRANSO_Engine
//{
//
//
//	public:
//
//	public:
//													SKTRANSO_Engine							();
//		virtual									   ~SKTRANSO_Engine							();
//
//
//};
//