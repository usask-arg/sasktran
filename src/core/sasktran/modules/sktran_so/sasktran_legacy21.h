#if !defined(SASKTRANV21_LEGACY_H_INCLUDED)
#define SASKTRANV21_LEGACY_H_INCLUDED


//#if defined(_MSC_VER)
//#if !defined(SKTRAN_TABLESUNIFORM_LIB)								// Must #define NXBASE_CORE_LIB in compiler switches when building nxbase, so nxbase compiles and links properly
//	#pragma message("Automatically linking the sasktran legacy configuration classes v2.1 libraries into this project")
//	#if defined (_MANAGED)									// Is Visual C++ generating MSIL with /CLR
//		#if defined (NXDEBUG)
//			#pragma comment( lib, "sktran_TablesUniform21dclr")
//		#else
//			#pragma comment( lib, "sktran_TablesUniform21rclr")
//		#endif
//	#else
//		#if defined (NXDEBUG)
////			#pragma comment( lib, "sktran_TablesLegacy21_debug")
//		#else
////			#pragma comment( lib, "sktran_TablesLegacy21_release")
//		#endif
//	#endif
//#else
//	#pragma message("NOT automatically linking the sasktran legacy configuration v2.1 classes libraries into this project")
//#endif
//#endif



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height		2010-4-29*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height : public SKTRAN_SpecsUser_OpticalPropertiesGrid_V21
{
	private:
		std::vector<double>								m_heights;

	public:
														SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height(){};
		virtual										   ~SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height(){};
		bool											ConfigureOpticalPropertyShells	( const double* alts, size_t npts );

	public:
		virtual bool									CreateInternalSpecs( const SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21** internalspecs ) const;
};

class SKTRAN_TableRayLOSFactory_Legacy;

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_Diffuse_Legacy		2010-4-27*/
/** These are the diffuse specifciations used by the user. Configure the diffuse specications
*	as part of the larger engine specification and then call "engine".ConfigureModel
 *
 *	\par Usage
 *	This class was created to make a friendlier user interface. The class is
 *	configured and then placed as a parameter in the call to engine.ConfigureModel. After that
 *	the class instance can be destroyed as the specifications have been copied to an internal
 *	diffuse specifications class. This avoids the user from having to worry about object
 *	lifetime issues, nxUnknown and AddRef and Release which created all sorts of memory
 *	leaks and access violations in earlier code.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_Diffuse_Legacy : public SKTRAN_SpecsUser_Diffuse_V21
{
	private:
		bool													m_useUserDefinedSZATransmissiongrid;	//!< True if the user want to overide the default SZA grid used for the Solar transmission table.
		std::vector<double>										m_szatransmission;						//!< The grid of cos(sza) used to determine the placing of the solar transmission table grid.
		double													m_degreesPerSolarTransSza;				//!< The number of SZA degrees for each solar transmission entry
	

		//double												m_degreesPerDiffuseSZAProfile;			//!< The angular separation (degrees) of diffuse profiles in the Solar Zenith Angle  grid;= 0.015;
		size_t													m_numDiffuseProfiles;
		double													m_maxdiffusealtitudealongincomingrays;	//!< Ignore diffuse scattering contribution to "diffuse rays" from points above this altitude
		double													m_scatteringresolution_degrees;			//!< The degrees between the entires in the internal scattering matrix tables.
		std::vector<double>										m_incomingazimuth;						//!< The azimuths that will be used to create the diffuse grid
		std::vector<double>										m_diffusealtitudes;						//!< The altitudes used to generate the diffuse profiles. (this might be tweaked up when create the internal diffuse specifications)
		bool													m_useUserDefinedincomingzenith;			//!< Use the user defined incoming zeniths rather groundRes, horizonRes and atmosRes
		std::vector<double>										m_userdefinedincomingzenith;			//!< The user defined incoming zenith angles (in degrees).
		size_t													m_groundRes;							//!< The number of incoming zenith angles used to cover the ground ((usually coarse to medium)
		size_t													m_horizonRes;							//!< The number of incoming zenith angles used to cover the horizon region (usually fine scale)
		size_t													m_atmosRes;								//!< The number of incoming zenith angles used to cover the zenith directions looking upward.
		const SKTRAN_RayTracingRegionManager*					m_rayregionmanager;
		SKTRAN_UnitSphere_V2*									m_unitsphere;

		bool													m_use_losinternalsinglescattertables;	//!< If true then use internal single scatter tables for each line of sight. Helps provide speed up for single scatter.
		bool													m_configuredForMC;
		bool													m_shiftedHorizon;						//!< True if using a variable horizon of 20 degree width.



	private:
		bool													IsGroundPoint						( double h ) const;

		bool													ConfigureDefaults					();
		bool													ConfigureEvenSpacedShells			( double minshellheight_meters, double shellwidth, double maxshell);

		bool													MakeDiffuseRadiiGrid				( SKTRAN_GridDefDiffuseHeights_V21**         userdiffuseheights ) const;
		bool													MakeDiffuseSzaGrid					( SKTRAN_GridDefCosSZA_V21**                 userdiffuseSzagrid, SKTRAN_GridDefSLON_V21** userdiffuseslongrid ) const;

		bool													MakeIncomingUnitSpheres				(   SKTRAN_GridDefDiffuseHeights_V21*				internaldiffuseheights, 
																										std::vector< SKTRAN_UnitSphereLatLonGrid*>*		incomingzenithangle,
																									    std::shared_ptr< const SKTRAN_CoordinateTransform_V2>&			coordinatesystem
																									 ) const ;

//		bool													MakeIncomingZenithAngles			( std::vector< SKTRAN_GridDefDiffuseIncomingZenith_V2*>* incomingzenithangle, const SKTRAN_CoordinateTransform_V2* coordinatesystem) const;
		bool													MakeIncomingAzimuthAnglesV2			( SKTRAN_GridDefDiffuseIncomingAzimuth_V21** userincomingazimuthangle ) const;
		bool													MakeScatterAngleGrid				( SKTRAN_GridDefScatterAngle_V21**           userscatteranglegrid ) const;
		bool													MakeSolarTransmissionSZAGrid		( SKTRAN_GridDefCosSZA_V21**                 userSZAgrid, SKTRAN_GridDefSLON_V21** userSLONgrid ) const;
		bool													MakeOutboundUnitSphere				( SKTRAN_UnitSphere_V2**					userunitsphere ) const;
		bool													MakeMaxDiffuseAltitude				( double*									maxJAltitude ) const;
		bool													MakLOSInternalScatterFactory		( SKTRAN_TableRayLOSFactory_Legacy**		lossinglescatterfactory ) const;

		bool													SetUserDefinedIncomingZenithGrid	( SKTRAN_GridDefDiffuseIncomingZenith_V21*   zenithgrid, bool isgroundpoint) const;

		bool													SetIncomingZenithGrid				( bool										isgroundpoint,
																									  double									altitude,
																									  const SKTRAN_GridDefDiffuseHeights_V21*	internaldiffuseheights,
																									  SKTRAN_GridDefDiffuseIncomingZenith_V21*	zenithgrid,
																									  const SKTRAN_RayTracingRegionManager*		rayregionmanager,
																									  const SKTRAN_CoordinateTransform_V2*		coordinatesystem ) const;
		bool													SetIncomingZenithGridShiftedHorizon	( bool										isgroundpoint,
																									  double									altitude,
																									  const SKTRAN_GridDefDiffuseHeights_V21*	internaldiffuseheights,
																									  SKTRAN_GridDefDiffuseIncomingZenith_V21*	zenithgrid,
																									  const SKTRAN_RayTracingRegionManager*		rayregionmanager,
																									  const SKTRAN_CoordinateTransform_V2*		coordinatesystem ) const;

	public:
																SKTRAN_SpecsUser_Diffuse_Legacy	( const SKTRAN_RayTracingRegionManager* rayregionmanager  );
															   ~SKTRAN_SpecsUser_Diffuse_Legacy	();
		bool													ConfigureDiffuseAltitudeResolution			( const double* altitudes, size_t npts );
		bool													ConfigureCosScatteringAngleGrid				( double degrees_resolution );
		bool													ConfigureIncomingZenithResolutions			( size_t groundResolution, size_t horizonResolution, size_t atmosResolution, bool shiftedHorizon );
		bool													ConfigureUserDefinedIncomingZenithAngle		( double* zenithradians, size_t numzen);
		bool													ConfigureIncomingAzimuthResolution			( const double* azimuths, size_t npts );
		//bool													ConfigureDegreesSZAPerDiffuseProfile		( double diffusedeltasza_degrees );
		bool													ConfigureNumberDiffuseProfiles				( size_t numProfiles);
		bool													ConfigureDegreesSZAPerSolarTransmission		( double deltasza_degrees );
		bool													ConfigureMaxJsourceAltitudeAlongDiffuseRays	( double heightm);
		bool													ConfigureUserDefinedTransmissionSZAGrid		( double* cossza, size_t numsza);
		bool													ConfigureOutboundUnitSphere					(  SKTRAN_UnitSphere_V2* sphere );
		bool													UseGlobalSolarTransmissionTableForSingleScatter	( bool useglobal);
		void													ConfigureForMonteCarlo						( bool isUsedForMC ) { m_configuredForMC = isUsedForMC;}
		bool													IsConfiguredForMonteCarlo					() const { return m_configuredForMC; }
		bool													IsProperlyDefined							() const;




		virtual bool 											CreateInternalSpecs				(       SKTRAN_SpecsInternal_Diffuse_V21**   internalspecs,
																								  std::shared_ptr< const SKTRAN_CoordinateTransform_V2> & coordinatesystem) const override;

};



/*-----------------------------------------------------------------------------
 *					class SKTRAN_SpecsUser_QuadratureLegacy_V21		2010-5-20*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsUser_Quadrature_V21_Legacy : public SKTRAN_SpecsUser_Quadrature_V21
{
	public:
																SKTRAN_SpecsUser_Quadrature_V21_Legacy	(){};
		virtual												   ~SKTRAN_SpecsUser_Quadrature_V21_Legacy	(){};
		virtual bool											CreateInternalSpecs( const SKTRAN_SpecsInternal_Quadrature_V21** internalspecs ) const override;

};


/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser_Legacy		2008-1-11*/
/** A class used to define the ray tracing specifications used by the Sasktran
 *	model.  This information includes the location of the observer and sun the
 *	look direction, shell heights and so on
 *
 *	\par Object lifetimes
 *	The class has some rules regarding object lifetime.  The class keeps pointers
 *	to each grid object.  When the dirty flag is set all the references to the
 *	objects are decreased and a new set of objects created.
 *
 *	This means grid objects referenced in the Sasktran engine will not change
 *	if the objects in the grid specifications are changed as long as the Sasktran
 *	code uses the proper reference counting scheme (AddRef/Release).  When the
 *	user calls one of the Get_XXXX methods they will get an object without
 *	a reference count.  The user either has to immediately copy information from
 *	the object or place a reference count on it (with AddRef). All calls to
 *	AddRef() must be matched by a corresponding call to Release().
 **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_SpecificationsUser_Legacy : public SKTRANSO_SpecificationsUser
{
	private:
		bool																m_atmosphericemissions_enabled;
		SKTRAN_SpecsUser_Diffuse_Legacy										m_diffusespecs;
		SKTRAN_SpecsUser_RayTracing_V21										m_raytracingspecs;
		SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height					m_opticalpropspecs;
		SKTRAN_RayTracingRegionManager										m_rayregionmanager;
		SKTRAN_SpecsUser_GroundPoints_V21_Lambertian						m_groundpointspecs;
		SKTRAN_SpecsUser_Quadrature_V21_Legacy								m_quadraturespecs;

	private:
		bool																CheckSpecsAreProperlyDefined() const;


	public:
																			SKTRANSO_SpecificationsUser_Legacy		();
		virtual 														   ~SKTRANSO_SpecificationsUser_Legacy		();

		bool																ConfigureEvenSpacedShells		( double minshellheight_meters,
																											  double shellwidth,
																											  double maxshell);

		bool																ConfigureUserDefinedShells		( const std::vector<double>& shellalts );

		virtual const SKTRAN_SpecsUser_RayTracing_V21*						RayTracingSpecs		() const override { return &m_raytracingspecs;}
		virtual const SKTRAN_SpecsUser_Diffuse_V21*							DiffuseSpecs		() const override { return &m_diffusespecs;}
		virtual const SKTRAN_SpecsUser_GroundPoints_V21*					GroundPointSpecs	() const override { return &m_groundpointspecs;}
		virtual const SKTRAN_SpecsUser_Quadrature_V21*						QuadratureSpecs		() const override { return &m_quadraturespecs;}
		virtual const SKTRAN_SpecsUser_OpticalPropertiesGrid_V21*			OpticalGridSpecs	() const override { return &m_opticalpropspecs;}
		virtual bool														CreateCoordinateSystem ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* usercoords, double groundht, double toaht ) const override;
		virtual bool														UpdateUndefinedParametersFromLinesOfSight( const SKTRAN_LineOfSightArray_V21& linesofsight ) override;
		virtual bool														AtmosphericEmissionsAreEnabled () const override  { return m_atmosphericemissions_enabled;}


		SKTRAN_SpecsUser_Diffuse_Legacy*									DiffuseSpecificationsVar	()				{ return &m_diffusespecs;}
		SKTRAN_SpecsUser_RayTracing_V21*									RayTracingSpecsVar			()				{ return &m_raytracingspecs;}
		SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height*					OpticalTableSpecsVar		()				{ return &m_opticalpropspecs;}
		SKTRAN_SpecsUser_GroundPoints_V21_Lambertian*						GroundPointSpecsVar			()				{ return &m_groundpointspecs;}
		virtual SKTRAN_RayTracingRegionManager*								RayTracingRegionManagerVar	() override		{ return &m_rayregionmanager;}

};

class SKTRAN_GridSpecificationsLegacy_MC_V21 : public SKTRANSO_SpecificationsUser_Legacy
{
	public:
		SKTRAN_GridSpecificationsLegacy_MC_V21		(){ ((SKTRAN_SpecsUser_Diffuse_Legacy*) DiffuseSpecs())->ConfigureForMonteCarlo(true); }
};

#endif
