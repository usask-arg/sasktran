#include "modules/sktran_so/sasktranv21_internals.h"
#include "include/sktran_occ_internals.h"
#include "include/sktran_curvedrays_ace.h"
#include "../sktran_common/sktran_common.h"



/*-----------------------------------------------------------------------------
 *					struct SKOCCULT_LOSFromTangentPoint		 2014- 5- 22*/
/** A small structure to hold definitions of lines of sight for the occultation
 *	engine using tangent altitude and observer altitude.
 **/
/*---------------------------------------------------------------------------*/

struct SKOCCULT_LOSFromTangentPoint
{
	double		TangentAltitude;
	double		ObserverAltitude;
};


/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User		 2014- 4- 24*/
/** The external specifications used by users (and SasktranIF) to configure
 *	the occultation engine. The main purpose of this class is to provide the
 *	user with a configuration class that is approriate and has easy lifetime
 *	management. 
 **/
/*---------------------------------------------------------------------------*/

class SKOCCULT_Specs_User : public SKTRAN_SpecsUser_Base
{
	private:
		std::list<SKOCCULT_LOSFromTangentPoint>								m_tanpoint_los;					// Lines of sight using tangent altitude and observer altitude.
		GEODETIC_INSTANT													m_manualrefpt;					// The manual reference point. Required  if m_tanpoint_los is used
		nxVector															m_manualsun;					// The manual sun location. May be left unset, in which case Sun is taken using mjd of reference point.
		SKTRAN_SpecsUser_RayTracing_V21										m_raytracingspecs;				// Stores the ray tracing grid
		SKTRAN_RayTracingRegionManager										m_rayregionmanager;				// Stores the reference point and coordinate objects
		SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height					m_opticalpropspecs;				// The optical properties grid
//		SKTRAN_SpecsUser_Quadrature_V21_Legacy								m_quadraturespecs;				// May not need this for first version

	private:
		bool																CheckSpecsAreProperlyDefined				() const;

	public:
																			SKOCCULT_Specs_User						();
		virtual 														   ~SKOCCULT_Specs_User						();
		bool																ConfigureEvenSpacedShells					( double minshellheight_meters, double shellwidth, double maxshell);
		bool																ConfigureUserDefinedShells					( const std::vector<double>& shellalts );
		bool																SetReferencePointManually					( double latitude, double longitude, double mjd );
		bool																SetSunManually								( const nxVector& sun );
		bool																AddLineOfSightFromTangentAlt				( double tangentaltitude, double observeraltitude);
		bool																ManualRefPtIsDefined						() const;
		const std::list<SKOCCULT_LOSFromTangentPoint>&						TangentPointLinesOfSight					() const { return m_tanpoint_los;}
		virtual const SKTRAN_SpecsUser_RayTracing_V21*						RayTracingSpecs								() const { return &m_raytracingspecs;}
//		virtual const SKTRAN_SpecsUser_Quadrature_V21*						QuadratureSpecs								() const { return &m_quadraturespecs;}
		virtual const SKTRAN_SpecsUser_OpticalPropertiesGrid_V21*			OpticalGridSpecs							() const { return &m_opticalpropspecs;}
		virtual bool														CreateCoordinateSystem						( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* coords) const;
		virtual bool														UpdateUndefinedParametersFromLinesOfSight	( const SKTRAN_LineOfSightArray_V21& linesofsight );

		SKTRAN_SpecsUser_RayTracing_V21*									RayTracingSpecsVar			()			{ return &m_raytracingspecs;}
		SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height*					OpticalTableSpecsVar		()			{ return &m_opticalpropspecs;}
		SKTRAN_RayTracingRegionManager*										RayTracingRegionManagerVar	()			{ return &m_rayregionmanager;}
};


/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_Internal		 2014- 4- 24*/
/** The internal specifications used for configuration by the occultation
 *	engine. The engine copies the user specs over to the internal specs.
 *	Apart from providing engine configuration the class also has the
 *	following goals (a) allow users to develop an array of configuration
 *	classes derived from SKOCCULT_Specs_User and (b) provide more
 *	sophisticated lifetime management of the configuration objects without placing
 *	excessive demands on the user to keep the objects alive.
 **/
/*---------------------------------------------------------------------------*/

class SKOCCULT_Specs_Internal
{
		std::list<SKOCCULT_LOSFromTangentPoint>					m_tanpoint_los;					// Lines of sight using tangent altitude and observer altitude.
		const SKTRAN_SpecsInternal_RayTracing_V21*				m_raytracingspecs;					// The ray tracing shells
//		const SKTRAN_SpecsInternal_Quadrature_V21*				m_quadraturespecs;					// MAy not need this in first version
		const SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21*	m_opticaltablespecs;				// The optical properties grid
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	m_coordinatesystem;					// Stores the reference point and coordinate object

	private:
		void													ReleaseResources			 ();

	public:
																SKOCCULT_Specs_Internal		 ();
		virtual												   ~SKOCCULT_Specs_Internal		();
		bool													Initialize					 (  const SKOCCULT_Specs_User* userspecifications );
		const std::list<SKOCCULT_LOSFromTangentPoint>&			TangentPointLinesOfSight	 () const { return m_tanpoint_los;}
		const SKTRAN_SpecsInternal_RayTracing_V21*				RayTracingSpecs				 () const { return m_raytracingspecs;}
//		const SKTRAN_SpecsInternal_Quadrature_V21*				QuadratureSpecs				 () const { return m_quadraturespecs;}
		const SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21*	OpticalTableSpecs			 () const { return m_opticaltablespecs;}
		const SKTRAN_CoordinateTransform_V2*					CoordinateSystemPtr			 () const { return m_coordinatesystem.get();}
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	CoordinateSystemObject		 () const { return m_coordinatesystem;}
};



/*-----------------------------------------------------------------------------
 *					SKOCCULT_TableOpticalProperties_Base		 2014- 4- 24*/
/** Base class for the atmospheric optical properties used in Occultation
 *	engines. Occultation engines are different as they want to store
 *	optical properties for multiple wavelengths (eg 2000 wavelengths ).
**/
/*---------------------------------------------------------------------------*/

class SKOCCULT_TableOpticalProperties_Base
{
	public:
															SKOCCULT_TableOpticalProperties_Base	(){};
		virtual 										   ~SKOCCULT_TableOpticalProperties_Base	(){};
		virtual bool										ConfigureOptical						( const std::vector<double>& wavenumber, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) = 0;
		virtual double										TotalExtinctionPerCM					( const HELIODETIC_POINT& point, size_t wavenumber_index  ) const = 0;
		virtual bool										ConfigureGeometry						( const SKOCCULT_Specs_Internal* specs ) = 0;
};

/*-----------------------------------------------------------------------------
 *					SKOCCULT_OpticalProperties1D_HeightWavelength		 2014- 4- 24*/
/** A class used to cache the optical properties of the atmosphere for occultation
 *	engines. This class implements optical properties which are only a function
 *	of altitude. The big difference in this table is that it only stores the extinction
 *	and that it calculates the table for a large array of wavelengths simultaneously.
 *	This gives a significant speed advantage for calculations involving the Voigt function,
 *	which is a big part of the occultation engine domain.
 **/
/*---------------------------------------------------------------------------*/

class SKOCCULT_OpticalProperties1D_HeightWavelength : public SKOCCULT_TableOpticalProperties_Base
{
	private:
		std::vector<double>							m_wavenumber;								//!< The array of wavelengths in this table
		std::vector<double>							m_heights;									//!< The array of heights in this table
		std::vector< std::vector<double> >			m_extinctionpercm;							//!< Extinction per cm as a function of m_extinctionpercm(wavelength)(altitude). Occultation engines only need extinction
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>		m_coords;									//!< The coordinates object. Has coordinate transform and reference point.

	public:
													SKOCCULT_OpticalProperties1D_HeightWavelength	();
		virtual 								   ~SKOCCULT_OpticalProperties1D_HeightWavelength	() override;
		virtual bool								ConfigureOptical								( const std::vector<double>& wavenumber, SKTRAN_AtmosphericOpticalState_V21& opticalstate )override;	
		virtual double								TotalExtinctionPerCM							( const HELIODETIC_POINT& point, size_t wavenumberindex  ) const override;
		virtual bool								ConfigureGeometry								( const SKOCCULT_Specs_Internal* specs )override;
};



/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise		2014-1-24*/
/** Implements a curved ray as a series of points that trace a ray using the
 *	technique outlined in "Ray tracing in a refracting spherically symmetric
 *	atmosphere". Thompson, Pepin and Simon.  J. Opt Soc. Am., 72, 11, 1982.
 **/
/*---------------------------------------------------------------------------*/

class SKOCCULT_RayGeometry_Curved_Piecewise
{
	private:
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	m_coords;
		HELIODETIC_VECTOR										m_observer;
		HELIODETIC_UNITVECTOR									m_look;
		std::vector<HELIODETIC_POINT>							m_path;							//!< The location of each trajectory point
		std::vector<double>										m_curvedpathdistance;			//!< The actual curved distance of the start of a cell from the origin. Accounts for ray curvature.

	public:
															SKOCCULT_RayGeometry_Curved_Piecewise		();
		virtual											   ~SKOCCULT_RayGeometry_Curved_Piecewise		();
		const SKTRAN_CoordinateTransform_V2*				CoordinatesPtr								() const { return m_coords.get();}
		void												ClearRay								();
		const HELIODETIC_POINT&								PathPoint								( const size_t& idx ) const		{ return m_path[idx]; }
		bool												Push_Point								( const HELIODETIC_POINT& pt, double s);						// Updates element in m_quadpoints
		bool												Initialize								( std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords, const HELIODETIC_VECTOR& observer, const HELIODETIC_UNITVECTOR& look);
		bool												ReserveTrajectoryStorage				( size_t n );
		size_t												NumCells								() const	{ return m_curvedpathdistance.size() > 0 ?  m_curvedpathdistance.size()-1 : 0;}
		size_t												NumTrajectoryPoints						() const	{ return m_curvedpathdistance.size();}
		bool												GetCellLength							( size_t cellidx, double*           length	  ) const;
		bool												GetCellMidPoint							( size_t cellidx, HELIODETIC_POINT* pt        ) const;
		bool												GetCellStartLocation					( size_t cellidx, HELIODETIC_POINT* pt        ) const;
		bool												GetCellRayUnitVector					( size_t cellidx, HELIODETIC_UNITVECTOR* look ) const;
		const HELIODETIC_VECTOR&							Observer								() const { return m_observer;}
		const HELIODETIC_UNITVECTOR&						LookVectorAtObserver					() const { return m_look;}
};


/*-----------------------------------------------------------------------------
 *							SKTRAN_RayTracer_Curved_ThomPepSim		2014-01-31*/
/** 
 *	Wrapper class that interfaces the Thompson, Pepin and Simon curved
 *	ray tracing class into the Sasktran V3 paradigm. This curved ray tracer
 *	is not a "fast" ray tracer. It was originally designed in 2014 to service the lines of
 *	sight used in the sktran_occ engine. The code is based upon the refrac.F ray tracing
 *	algorithm used in the ACE-FTS code. I have shown that when the raytracing grid
 *	and atmospheric state are made to match the ACE system the agreement is within
 *	2 microns. 
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayTracer_Curved_ThomPepSim
{
	private:
		SKTRAN_RayRefracted_ThomPepSim				m_curvedraytracer;

	public:
													SKTRAN_RayTracer_Curved_ThomPepSim		( );
		virtual									   ~SKTRAN_RayTracer_Curved_ThomPepSim		( );
		SKTRAN_RayRefracted_ThomPepSim&				CurvedRayTracer							() { return m_curvedraytracer;}
		bool										Initialize								( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>			coords,
																							  std::shared_ptr<const SKTRAN_GridDefRayTracingShells_V21>		raytracinggrid,
																							  skClimatology*								raytracing_atmosphericstate,
																							  const GEODETIC_INSTANT&						riprofile_location,		
																							  double										raytracing_wavenumber);

		bool										TraceRay								( SKOCCULT_RayGeometry_Curved_Piecewise* ray ) const;
};


/*-----------------------------------------------------------------------------
 *					SKOCCULT_TableLinesOfSight		 2014- 4- 27*/
/** **/
/*---------------------------------------------------------------------------*/

class SKOCCULT_TableLinesOfSight
{
	private:
		std::vector< SKOCCULT_RayGeometry_Curved_Piecewise>						m_raytrajectory;
		typedef std::vector< SKOCCULT_RayGeometry_Curved_Piecewise>::iterator	iterator;


	public:
														SKOCCULT_TableLinesOfSight			();
		virtual										   ~SKOCCULT_TableLinesOfSight			();
		void											ClearRayTrajectories				();

		bool											SetLinesOfSight						( const SKTRAN_LineOfSightArray_V21&	observerlinesofsight,
																							  std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	coordinates,
																							  SKTRAN_RayTracingRegionManager*	    raytracingregionmanager);

		bool											AddLinesOfSightFromTangentAltitude	(	const std::list<SKOCCULT_LOSFromTangentPoint>&			tanpoint_los,
																								SKTRAN_RayRefracted_ThomPepSim&							raytracer,
																								std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	coordinates);

		size_t											NumRays								()				const		{ return m_raytrajectory.size();}
		const SKOCCULT_RayGeometry_Curved_Piecewise&	RayAt								(size_t idx)	const		{ return m_raytrajectory.at(idx);}
		SKOCCULT_RayGeometry_Curved_Piecewise&			RayAtVar							(size_t idx)				{ return m_raytrajectory.at(idx);}
};

/*-----------------------------------------------------------------------------
 *					class SKOCCULT_OCC_Engine		 2014- 4- 27*/
/** **/
/*---------------------------------------------------------------------------*/

 class SKOCCULT_OCC_Engine : public SKTRAN_Engine_Base
{
	private:
		double											m_raytracingwavenumber;				// The wavenumber used for ray tracing (default 2400 cm-1)
		double											m_maxpathlength;					// The maximum path length along any section of the ray (default 10,000 meters)
		SKTRAN_RayTracer_Curved_ThomPepSim				m_raytracer;						// The curved ray tracer.
		SKOCCULT_TableLinesOfSight						m_linesofsight;						// The lines of of sight and their curved trajectories
		SKOCCULT_Specs_Internal							m_modelspecifications;				// The various specifications
		SKOCCULT_OpticalProperties1D_HeightWavelength	m_opticalpropertiestable;			// The height grid of optical properties (extinction per cm at various altitudes and wavenumbers.
		bool											m_geometrychanged;


	private:
		bool											TraceLineOfSightRays					(	SKTRAN_AtmosphericOpticalState_V21*		opticalstate);
		bool											CalculateOpticalPropertiesTable			(	double									wavelen,
																									SKTRAN_AtmosphericOpticalState_V21*		opticalstate,
																									bool								    userupdateclimatology );
		GEODETIC_INSTANT								ReferencePoint							() const;
		bool											IntegrateExtinctionAlongRay				( const SKOCCULT_RayGeometry_Curved_Piecewise& ray, const std::vector<double>& wavenumber, std::vector<double>* extinction);

	public:
														SKOCCULT_OCC_Engine						();
		virtual										   ~SKOCCULT_OCC_Engine						();
		bool											SetRayTracingWavenumber					( double wavenumber );
		size_t											NumRays									() const { return m_linesofsight.NumRays();}

		virtual bool									ConfigureModel							(   SKTRAN_SpecsUser_Base&				modelspecifications,
																									const SKTRAN_LineOfSightArray_V21&		linesofsight,
																									size_t									numthreads ) override;


		virtual bool									CalculateMultiWavelengthRadiance		(	std::vector< std::vector<SKTRAN_StokesScalar> >*	losradiance,
																									const std::vector<double>&							wavelen,
																									size_t												numordersofscatter,
																									SKTRAN_AtmosphericOpticalState_V21*					opticalstate,
																									std::vector< std::vector<skRTStokesVector> > *      losvector=nullptr, 
																									bool												updateclimatology = false,
																									SKTRAN_DiagnosticInterface*							diag = NULL);

};
