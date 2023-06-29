
class SKTRANSO_RayInternalGeometry;

/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayStorage_InternalJindex		2008-2-7*/
/** A class used to store the line of sight geometry.  This class uses the
 *	methods from SKTRANSO_RayStorageBaseGeometry to perform the raytracing code 
 *	but provides the tables for calculating the source functions due to
 *	solar (single) scatter, diffuse (multiple) scatter and ground albedo
 *	at quadrature points along the ray.
 *
 *	\par The Diffuse J
 *	Class member m_diffuseJ stores indexes and weights into the outbound 
 *	rays of the table of diffuse profile. Each entry indexes one outbound ray
 *	at one point in the diffuse point table.
 *
 *	\par Solar Transmission J
 *	Class member m_singlescatterJ stores the indexes need to calculate
 *	the transmission of the atmosphere from the sun to a quadrature point
 *	along the ray.
 *
 *	\par Ground Point Diffuse Upwelling
 *	Class member m_groundpointJ stores the indexes required to calculate the
 *	albedo scattering terms off the ground.  It interpolates the
 *	table of ground points in the table of diffuse profiles to obtain higher order
 *	scattering off the ground.
 *
 *

 **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_RayStorage_InternalJindex
{
	private:
		SKTRANSO_JIndexArray						m_diffuseJ;						//!< The table indexing the diffuse source functions for each quadrature point along the ray
		SKTRANSO_JIndexArray						m_singlescatterJ;				//!< The table indexing the solar transmission to each quadrature point along the ray.
		SKTRANSO_JIndexArray						m_groundpointMSJ;				//!< The table indexing the multiple scatter term from the ground point at the end of the ray (if any).
		SKTRANSO_JIndexArray						m_groundpointSSJ;				//!< The table indexing the single scatter term from the ground point at the end of the ray (if any).
		SKTRANSO_JIndexArray						m_emissionJ;					//!< The table indexing the photochemical and thermal emission along the ray (if any).
		SKTRANSO_JIndexArray						m_groundppointemissionJ;		//!< The table indexing the photochemical and thermal emission from the ground point if any

	protected:
		bool										CreateSingleScatterJIndexTables		( SKTRANSO_Quadrature_TLS_V21*	quadrature, const SKTRANSO_RayInternalGeometry* ray);
		bool										CreateMultipleScatterJIndexTables	( SKTRANSO_Quadrature_TLS_V21*	quadrature, const SKTRANSO_RayInternalGeometry* ray, bool singlescatter);

	public:
													SKTRANSO_RayStorage_InternalJindex		();
		virtual									   ~SKTRANSO_RayStorage_InternalJindex		();

		bool										CreateJIndexTables						( SKTRANSO_Quadrature_TLS_V21*	quadrature, const SKTRANSO_RayInternalGeometry* ray, bool singlescatter);
		const SKTRANSO_JIndexArray&					DiffuseJIndex							() const	{ return m_diffuseJ;} 
		const SKTRANSO_JIndexArray&					GroundPointMSJIndex						() const	{ return m_groundpointMSJ;}
		const SKTRANSO_JIndexArray&					SingleScatterJIndex						() const	{ return m_singlescatterJ;}
		const SKTRANSO_JIndexArray&					GroundPointSSJIndex						() const	{ return m_groundpointSSJ;}
		const SKTRANSO_JIndexArray&					AtmosphericEmissionJIndex				() const	{ return m_emissionJ;}
		const SKTRANSO_JIndexArray&					GroundPointEmissionJIndex				() const	{ return m_groundppointemissionJ;}
		SKTRANSO_JIndexArray&						DiffuseJIndexVar						()			{ return m_diffuseJ;} 
		SKTRANSO_JIndexArray&						GroundPointMSJIndexVar					()			{ return m_groundpointMSJ;}
		SKTRANSO_JIndexArray&						SingleScatterJIndexVar					()			{ return m_singlescatterJ;}
		SKTRANSO_JIndexArray&						GroundPointSSJIndexVar					()			{ return m_groundpointSSJ;}
		SKTRANSO_JIndexArray&						AtmosphericEmissionJIndexVar			()			{ return m_emissionJ;}
		SKTRANSO_JIndexArray&						GroundPointEmissionJIndexVar			()			{ return m_groundppointemissionJ;}
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Base		 2014- 11- 25*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_RayInternalGeometry
{
	private:
		SKTRANSO_RayStorage_InternalJindex					m_indices;				// The indices associated with this ray
		std::unique_ptr< SKTRAN_RayOptical_Base>			m_ray;					// The ray trajectory, includes the geometry and optical/wavelength component of the ray

	public:
															SKTRANSO_RayInternalGeometry(){}
		virtual											   ~SKTRANSO_RayInternalGeometry(){};
															SKTRANSO_RayInternalGeometry(const SKTRANSO_RayInternalGeometry&  moveother)  { NXASSERT(( !moveother.m_ray )); }
		void												operator=					(const SKTRANSO_RayInternalGeometry&  moveother ) { NXASSERT(( !moveother.m_ray )); };

		const SKTRANSO_RayStorage_InternalJindex*			Indices			() const											{ return &m_indices;}
		SKTRANSO_RayStorage_InternalJindex*					IndicesVar		()													{ return &m_indices;}
		const SKTRAN_RayStorage_Base*						Storage			() const											{ return m_ray->Storage();}
		SKTRAN_RayStorage_Base*								StorageVar		()													{ return m_ray->StorageVar();}
		const SKTRAN_RayOptical_Base*						Ray				() const											{ return m_ray.get();}
		SKTRAN_RayOptical_Base*								RayVar			()													{ return m_ray.get();}
		bool												AssignRay		( std::unique_ptr< SKTRAN_RayOptical_Base> ray )	{ m_ray = std::move(ray); return true;}
		virtual const SKTRANSO_JindexTableBase*					InternalSolarTransmissionTable					() const			{ return nullptr;}

};



/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayLOSGeometry_V21		2010-6-16*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_RayLOSGeometry_V21 : public SKTRANSO_RayInternalGeometry
{
	private:
		SKTRANSO_TableRayLOS*							m_singlescattertable;


	public:
													SKTRANSO_RayLOSGeometry_V21						();
		virtual									   ~SKTRANSO_RayLOSGeometry_V21						();
		bool										SetInternalSolarTransmissionTable				( SKTRANSO_TableRayLOS* table );
		bool										ConfigureInternalSolarTransmissionTableGeometry	( SKTRANSO_Quadrature_TLS_V21* quadrature );	
		bool										ConfigureInternalSolarTransmissionTableOptical	( SKTRANSO_Quadrature_TLS_V21* quadrature  );
		virtual const SKTRANSO_JindexTableBase*			InternalSolarTransmissionTable					() const override; // { return m_singlescattertable;}
};


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayInternalOptical		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_RayInternalOptical
{
	private:
		const SKTRANSO_RayInternalGeometry*		m_geometryray;
		double									m_totaltransmission;

	public:
												SKTRANSO_RayInternalOptical();
		virtual								   ~SKTRANSO_RayInternalOptical();
		bool									AttachToGeometry					( const SKTRANSO_RayInternalGeometry* geometry );
		bool									ConfigureOptical					( SKTRANSO_Quadrature_TLS_V21* quadrature, bool totaltransmissiononly, bool usecachedtransmission  );
		const SKTRANSO_RayInternalGeometry*		GeometryRay							() const { return m_geometryray;}
		SKTRAN_StokesScalar						TotalTransmissionAlongRay			() const { return  SKTRAN_DBL_TO_STOKES_SCALAR(m_totaltransmission);}

};


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayInternalDiffuseOptical_V2		2007-11-23*/
/** A class used to represent a single line of sight with given optical
 *	properties.  This class uses a line of sight geometry as the basi
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayInternalDiffuseOptical_V2 : public SKTRANSO_RayInternalOptical
{
	private:
		const SKTRANSO_RayInternalGeometry*			m_internaldiffusegeometry;
		SKTRAN_JValueTable_V21						m_diffuseJ;							//!< The table indexing the 2nd order + higher diffuse source functions for each quadrature point along the ray
		SKTRAN_JValueTable_V21						m_groundpointMSJ;					//!< The table indexing the linear ground point	at the end of the ray (if any).
		SKTRAN_StokesScalar							m_singlescatterJ;					//!< The incoming signal from single scattered sunlight, calculated during ConfigureOptical
	public:
													SKTRAN_RayInternalDiffuseOptical_V2		();
		virtual									   ~SKTRAN_RayInternalDiffuseOptical_V2		();

		const SKTRANSO_RayInternalGeometry*			InternalDiffuseGeometry				() const { return m_internaldiffusegeometry;}
		bool										AttachToGeometry					( const SKTRANSO_RayInternalGeometry* los );

		virtual bool								ConfigureOptical					( SKTRANSO_Quadrature_TLS_V21*					quadrature,
																						  bool										singlescatter,
																						  bool										usecachedtransmission,
																						  bool										usecachedcellfactors );

		SKTRAN_StokesScalar							GetFirstOrderIncomingRadiance		() const	{ return m_singlescatterJ;}

		bool										CalculateTotalRadianceAtOrigin		( SKTRAN_StokesScalar* radiance );

};



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayInternalDiffuseOptical_V2		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayLOSOptical_V21 : public SKTRAN_RayInternalDiffuseOptical_V2
{
	private:
		SKTRANSO_RayLOSGeometry_V21*					m_losgeometryray;

	public:
													SKTRAN_RayLOSOptical_V21					();
		virtual									   ~SKTRAN_RayLOSOptical_V21					();
		bool										AttachToGeometry							( SKTRANSO_RayLOSGeometry_V21* los);
		virtual bool								ConfigureOptical							( SKTRANSO_Quadrature_TLS_V21*					quadrature,
																								  bool										singlescatter,
																								  bool										usecachedtransmission,
																								  bool										usecachedcellfactors );

};


