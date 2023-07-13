
class SKTRAN_SpecsInternal_RayTracing_V21;
class SKTRAN_SpecsInternal_Diffuse_V21;
class SKTRANSO_TableDiffusePoints;
class SKTRANSO_TableSolarTransmission;
class SKTRANSO_TableGroundPointDiffuse;
class SKTRAN_CoordinateTransform_V2;
class SKTRAN_SpecsInternal_GroundPoints_V21;
class SKTRAN_SpecsInternal_Quadrature_V21;
class SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21;


/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsUser_RayTracing_V21				2010-5-19*/
/** @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_RayTracing_V21
{
	private:
		std::vector<double>										m_raytracingshells;

	public:
		bool													ConfigureRayTracingShellAlts	( const double*  alts_meters, size_t npts );
		double													GroundAltitude					() const { return m_raytracingshells.front();}				
		double													TOAAltitude						() const { return m_raytracingshells.back();}				

	public:
																SKTRAN_SpecsUser_RayTracing_V21	();
		virtual												   ~SKTRAN_SpecsUser_RayTracing_V21	(){};
		virtual bool 											CreateInternalSpecs				( const SKTRAN_SpecsInternal_RayTracing_V21**   internalspecs ) const;
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsUser_Diffuse_V21				2010-5-19*/
/**  @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_Diffuse_V21
{
	public:
		virtual												   ~SKTRAN_SpecsUser_Diffuse_V21(){}
		virtual bool 											CreateInternalSpecs				(       SKTRAN_SpecsInternal_Diffuse_V21**   internalspecs,
																								  std::shared_ptr< const SKTRAN_CoordinateTransform_V2>&       coordinatesystem
																								) const = 0;

};



/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsUser_GroundPoints_V21_Lambertian		2010-5-20*/
/**  @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_GroundPoints_V21
{
	public:
		virtual			   ~SKTRAN_SpecsUser_GroundPoints_V21(){}
		virtual bool		CreateInternalSpecs							( const SKTRAN_SpecsInternal_GroundPoints_V21** internalspecs ) const  = 0; 

};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsUser_Quadrature_V21		2010-5-20*/
/**  @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_Quadrature_V21
{
	public:
		virtual												   ~SKTRAN_SpecsUser_Quadrature_V21() {}
		virtual bool											CreateInternalSpecs( const SKTRAN_SpecsInternal_Quadrature_V21** internalspecs ) const  = 0; 
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_OpticalPropertiesGrid_V21		2010-5-20*/
/**  @ingroup specs
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_OpticalPropertiesGrid_V21
{
	public:
		virtual			   ~SKTRAN_SpecsUser_OpticalPropertiesGrid_V21() {}
		virtual bool		CreateInternalSpecs( const SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21** internalspecs ) const  = 0; 
};




/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser		2010-3-3*/
/** @ingroup LegacySO
**/
/*---------------------------------------------------------------------------*/

class SKTRANSO_SpecificationsUser : public SKTRAN_SpecsUser_Base
{
	public:

		virtual const SKTRAN_SpecsUser_RayTracing_V21*						RayTracingSpecs			() const = 0;
		virtual const SKTRAN_SpecsUser_Diffuse_V21*							DiffuseSpecs			() const = 0;
		virtual const SKTRAN_SpecsUser_GroundPoints_V21*					GroundPointSpecs		() const = 0;
		virtual const SKTRAN_SpecsUser_Quadrature_V21*						QuadratureSpecs			() const = 0;
		virtual const SKTRAN_SpecsUser_OpticalPropertiesGrid_V21*			OpticalGridSpecs		() const = 0;
		virtual bool														CreateCoordinateSystem	                    (std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* usercoords, double groundht, double toaht ) const = 0;
		virtual bool														UpdateUndefinedParametersFromLinesOfSight	( const SKTRAN_LineOfSightArray_V21& linesofsight ) = 0;
		virtual bool														AtmosphericEmissionsAreEnabled () const = 0;
		virtual SKTRAN_RayTracingRegionManager*								RayTracingRegionManagerVar	   () = 0;



};






