/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_RayTracing_V21_General	2010-4-28*/
/**   @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_RayTracing_V21 : public nxUnknown
{
	private:
		std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21>	m_raytracingshells;

	public:
															SKTRAN_SpecsInternal_RayTracing_V21();
		virtual											   ~SKTRAN_SpecsInternal_RayTracing_V21();
		bool												ConfigureRayTracingShellAlts	( const double*  alts_meters, size_t npts );

	public:
		virtual std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21>	RayTracingShells () const { return m_raytracingshells;}
		virtual       size_t										MaxShellsAlongRay() const;

};





/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Diffuse_V21		2010-3-3*/
/**  @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_Diffuse_V21: public nxUnknown
{
	private:

	public:
															    SKTRAN_SpecsInternal_Diffuse_V21() {}
		virtual												   ~SKTRAN_SpecsInternal_Diffuse_V21() {}

		virtual const SKTRAN_GridDefScatterAngle_V21*			ScattterAngleGrid				() const = 0;

		virtual const SKTRAN_GridDefCosSZA_V21*					DiffuseProfileCosSZA			( ) const = 0;									//!< Return the cosine of Solar zenith angle of each diffuse profile
		virtual const SKTRAN_GridDefSLON_V21*					DiffuseProfileSLON				( ) const = 0;								//!< Return the solar longitude of each diffuse profile
		virtual const SKTRAN_GridDefDiffuseHeights_V21*			DiffuseHeights					( size_t profileidx) const = 0;					//!< Return the height profile of points for each diffuse profile

		virtual const SKTRAN_GridDefCosSZA_V21*					SolarTransmissionCosSZA			( ) const = 0;									//!< Return the cosine of Solar zenith angle of each diffuse profile
		virtual const SKTRAN_GridDefSLON_V21*					SolarTransmissionSLON			( ) const = 0;									//!< Return the solar longitude of each diffuse profile
		virtual const SKTRAN_TableRayLOSFactory*				LOSSingleScatterTableFactory	( ) const = 0;
	    
//		virtual const SKTRAN_GridDefDiffuseIncomingZenith_V2*	IncomingZenith					( size_t profileindex, size_t heightindex ) const = 0;
//	    virtual const SKTRAN_GridDefDiffuseIncomingAzimuth_V2*	IncomingAzimuth					( size_t profileindex, size_t heightindex ) const = 0;
	    virtual const SKTRAN_UnitSphereLatLonGrid*				IncomingUnitSphere				( size_t profileindex, size_t heightindex ) const = 0;
		virtual const SKTRAN_UnitSphere_V2*						OutboundUnitSphere				( size_t profileindex, size_t heightindex ) const = 0;
		virtual bool											FirstTimeInitializeDiffuseObjects() = 0;

		virtual double 											MaxDiffuseAltitude				( size_t profileindex) const = 0;

		virtual bool											CreateEmptyDiffuseTables		( SKTRANSO_TableDiffusePoints**	   diffusepointstable,
																								  SKTRANSO_TableSolarTransmission**  solartranstable,
																								  SKTRANSO_TableGroundPointDiffuse** groundpointstable,
																								  SKTRANSO_InternalEmissionPropertiesTable**		 emissionstable
																								) const = 0;

};



/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsInternal_GroundPoints_V21		2010-3-30*/
/**  @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_GroundPoints_V21: public nxUnknown
{
	private:
		SKTRAN_AlbedoBRDF_V2*							m_albedobrdf;

	public:
														SKTRAN_SpecsInternal_GroundPoints_V21   ( SKTRAN_AlbedoBRDF_V2* albedobrdf );
		virtual										   ~SKTRAN_SpecsInternal_GroundPoints_V21	();

		virtual const SKTRAN_AlbedoBRDF_V2*				AlbedoBRDF								() const;
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_GroundPoints_V21_Lambertian		2010-4-29*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_GroundPoints_V21_Lambertian : public SKTRAN_SpecsUser_GroundPoints_V21
{
	public:
																SKTRAN_SpecsUser_GroundPoints_V21_Lambertian	();
		virtual												   ~SKTRAN_SpecsUser_GroundPoints_V21_Lambertian	();
		virtual bool											CreateInternalSpecs( const SKTRAN_SpecsInternal_GroundPoints_V21** internalspecs ) const ; 

};



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_Quadrature_V21		2010-3-3*/
/**  @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_Quadrature_V21: public nxUnknown
{

	public:
																SKTRAN_SpecsInternal_Quadrature_V21() {}
		virtual												   ~SKTRAN_SpecsInternal_Quadrature_V21() {}

		virtual const SKTRAN_Quadrature_Factory_V21*				QuadratureFactory					() const = 0;
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21		2010-4-14*/
/**   @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21: public nxUnknown
{
	public:
																	SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21() {}
		virtual													   ~SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21() {}

		virtual const SKTRAN_GridDefOpticalPropertiesRadii_V21*		OpticalPropertiesGrid						  ()  const  = 0;
		virtual bool										CreateEmptyOpticalPropertiesTable			  ( SKTRAN_TableOpticalProperties_V21**	  opticalpropertiestable ) const = 0; 
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_Specifications_Tables_V2		2010-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_V21 : public SKTRAN_SpecsInternal_Base
{
		const SKTRAN_SpecsInternal_RayTracing_V21*				m_raytracingspecs;
		      SKTRAN_SpecsInternal_Diffuse_V21*					m_diffusespecs;
		const SKTRAN_SpecsInternal_GroundPoints_V21*			m_groundpointspecs;
		const SKTRAN_SpecsInternal_Quadrature_V21*				m_quadraturespecs;
		const SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21*	m_opticaltablespecs;
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>	m_coordinatesystem;
		std::shared_ptr< SKTRAN_RayFactory_Base>				m_rayfactory_transmissiononly;
		std::shared_ptr< SKTRAN_RayFactory_Base>				m_rayfactory_scattered;
		std::shared_ptr< SKTRAN_RayFactory_Base>				m_rayfactory_lineofsights;

	private:
		void													ReleaseResources			 ();

	public:
																SKTRAN_SpecsInternal_V21	 ();
		virtual												   ~SKTRAN_SpecsInternal_V21	 ();
		bool													Initialize					 (  const SKTRANSO_SpecificationsUser* userspecifications );
		const SKTRAN_SpecsInternal_RayTracing_V21*				RayTracingSpecs				 () const { return m_raytracingspecs;}
		const SKTRAN_SpecsInternal_Diffuse_V21*					DiffuseSpecs				 () const { return m_diffusespecs;}
			  SKTRAN_SpecsInternal_Diffuse_V21*					DiffuseSpecsVar				 ()		  { return m_diffusespecs;}
		const SKTRAN_SpecsInternal_GroundPoints_V21*			GroundPointSpecs			 () const { return m_groundpointspecs;}
		const SKTRAN_SpecsInternal_Quadrature_V21*				QuadratureSpecs				 () const { return m_quadraturespecs;}
		const SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21*	OpticalTableSpecs			 () const { return m_opticaltablespecs;}
		const SKTRAN_CoordinateTransform_V2*					CoordinateSystemPtr			 () const { return m_coordinatesystem.get();}
		std::shared_ptr< const SKTRAN_RayFactory_Base>			RayFactoryTransmissionOnly	 () const { return m_rayfactory_transmissiononly;}
		std::shared_ptr< const SKTRAN_RayFactory_Base>			RayFactoryScattered			 () const { return 	m_rayfactory_scattered;}
		std::shared_ptr< const SKTRAN_RayFactory_Base>			RayFactoryLOS				 () const { return m_rayfactory_lineofsights;}
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>	CoordinateSystemObject		 () const { return m_coordinatesystem;}
};


