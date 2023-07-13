//#pragma once


class SKTRAN_RayTracer_Base;

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Base_V3						2013-05-15*/
/** @ingroup rays
 *	@brief A base class for storing the trajectory of a ray through the atmosphere.
 *	The trajectory is normally stored as a series of discrete points along the ray.
 *
 *	The points are drawn by the ray tracers so the ray closely approximates a straight line between the points. The
 *	set of discrete points used to describe the path of the ray are called the 
 *	"trajectory" points although they are also called "quadrature" points
 *	in some parts of the code. The trajectory points break each ray
 *	into discrete sections, called cells, which are used by derived classes to
 *	integrate various parameters such as optical depth and source term integrals.
 *
 *	@par Storage Considerations
 *	Storage of information for the rays can quickly grow, especially
 *	in the successive orders engines. The ray generally has to store information so it can
 *	efficiently provide the engine with:
 *	- The Heliodetic location of each trajectory point (height, radius, sza, solar azimuth, x,y,z)
 *	- The Local effective "straight line" tangent point of each cell (height, radius, x,y,z)
 *	- The Local Ray direction of each cell (x,y,z)
 *	- The Global Tangent point of a ray; this is same as the local tangent point for straight rays but is different for curved rays
 *	- Path length of cell, not necessarily a straight line distance as it may be corrected value for ray curvature between cell endpoints
 *	- Entry point into the atmosphere
 *	- Exit point from atmosphere
 *	- Flag if ray hits the ground
 *	Ideally the engines will interact exclusively with rays
 *	through this class.
 *  The base class provides the basic framework for ray geometry, it contains
 *  basic information common to each ray, such as the observer location and 
 *  the look direction.
 *
 *  Derived classes are responsible for initialization and handling geometry
 *  information from the raytracer.
 */
/*---------------------------------------------------------------------------*/

//
//class SKTRAN_RayGeometry_Base
//{
//	private:
////		HELIODETIC_VECTOR									m_observer;					//!< The location of the obserevr
////		HELIODETIC_UNITVECTOR								m_look;						//!< The observers look direction
////		SKTRAN_RayStorage_Base*								m_trajectorystorage;		//!< The ray storage associated and owned/managed by this instance of this (derived) class
//	
//	protected:
////		virtual void										NotifyDerived_RayInvalid			() {}																//!< Notify derived classes that the ray has changed, typically #MoveObsever was called
////		void												SetTrajectoryStorage				( std::unique_ptr<SKTRAN_RayStorage_Base>	trajectorystorage);		//!< The ray storage associated and owned/managed by this instance of this (derived) class
//
//	private:
//
//	public:	// ---- Methods which should be treated as if protected or private:
////		bool												MoveObserver                        ( const HELIODETIC_VECTOR& observer, const HELIODETIC_UNITVECTOR& look );
////		bool												CalculateBaseLineTangentPointDetails( double groundaltitude, double* Robs, double* Tobs, double* Rt);
//
//	public:
//															SKTRAN_RayGeometry_Base				( SKTRAN_RayStorage_Base* trajectorystorage);
//		virtual											   ~SKTRAN_RayGeometry_Base				();
////		SKTRAN_RayGeometry_Base&							operator=							( SKTRAN_RayGeometry_Base&& moveother );								
////		const HELIODETIC_VECTOR&							GetObserver							() const	{ return m_observer; }		//!< Return the location of the observer in HELIODETIC coordinates.
////		const HELIODETIC_UNITVECTOR&						LookVector							() const	{ return m_look; }			//!< Return the look unit-vector away from the observer in HELIODETIC coordinates.
////		const SKTRAN_CoordinateTransform_V2*				Coordinates							() const	{ return m_trajectorystorage->GetCoordsPtr(); } //!< Return the coordinate object that transforms from heliodetic corodinates to geographic corodinates etc.
////		SKTRAN_RayStorage_Base*								StorageAccessVar                    ()			{ return m_trajectorystorage;} //!< Return the ray's trajectory storage object and allow the caller to have write access.
////		const SKTRAN_RayStorage_Base*						StorageAccess                       () const	{ return m_trajectorystorage;} //!< Return the ray's trajectory storage object with only read access
////		bool												IsDefined							() const	{ return NXFINITE( m_observer.X());} //!<  Return true if the trajectory is defined.
//
//	public:
//
//
////		static bool											GetQuadratureInterpParamsFromPointsAndLook	(	double* r0, double* r1,
////																											double* t0, double* t1, double* rt, 
////																											const HELIODETIC_POINT& startpoint,
////																											const HELIODETIC_POINT& endpoint,
////																											const HELIODETIC_UNITVECTOR& look );
//
//
//};


class SKTRAN_RayFactory_Base;

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayOptical_Base		2013-05-22*/
/** @ingroup rays
 *	Base class for the optically dependant parts of rays.  The optical ray
 *  contains the geometry ray and quadrature information.  
 **/
/*--------------------------------------------------------------------------- */

class SKTRAN_RayOptical_Base
{
	private:
		HELIODETIC_VECTOR								m_observer;					//!< The location of the obserevr
		HELIODETIC_UNITVECTOR							m_look;						//!< The observers look direction
		SKTRAN_RayStorage_Base*							m_trajectorystorage;		//!< The ray storage associated and owned/managed by this instance of this (derived) class
		const SKTRAN_CoordinateTransform_V2*			m_coords;					//!< The corodinate object. This is acquired from the trajectory storage during InitializeStorage
		std::vector< double>							m_opticaldepth;				//!< An array to store the optical depth along the ray.
		double											m_wavelength;

	protected:
		void											InitializeStorage					(SKTRAN_RayStorage_Base* trajectorystorage );
		virtual void									NotifyDerived_RayInvalid			() {} //!< Notify derived classes that the ray has changed, typically #MoveObsever was called

	public:
		bool											MoveObserver                        (const HELIODETIC_VECTOR&    observer, const HELIODETIC_UNITVECTOR& look );
		SKTRAN_RayStorage_Base*							StorageVar							() { return m_trajectorystorage;}

	public:
														SKTRAN_RayOptical_Base				();
		virtual										   ~SKTRAN_RayOptical_Base				();
		const SKTRAN_RayStorage_Base*					Storage								() const								{ return m_trajectorystorage;}
		size_t											GetNumQuadraturePoints				() const								{ return m_trajectorystorage->NumQuadraturePoints();}
		size_t											GetNumCells							() const								{ return m_trajectorystorage->NumCells();}

		std::vector<double>*							OpticalDepthArrayVar				()		 { return &m_opticaldepth;}
		const std::vector<double>&						OpticalDepthArray					() const { return m_opticaldepth;}
		double											TotalOpticalDepth					() const { return m_opticaldepth.back();}
		bool											CalculateBaseLineTangentPointDetails( double groundaltitude, double* Robs, double* Tobs, double* Rt);
		const HELIODETIC_VECTOR&						GetObserver							() const	{ return m_observer; }		//!< Return the location of the observer in HELIODETIC coordinates.
		const HELIODETIC_UNITVECTOR&					LookVector							() const	{ return m_look; }			//!< Return the look unit-vector away from the observer in HELIODETIC coordinates.
		const SKTRAN_CoordinateTransform_V2*			Coordinates							() const	{ return m_coords; } //!< Return the coordinate object that transforms from heliodetic corodinates to geographic corodinates etc.
		bool											IsDefined							() const	{ return NXFINITE( m_observer.X());} //!<  Return true if the trajectory is defined.

		bool											SetWavelength						(double wavelength) { m_wavelength = wavelength; return true;}
		double											GetWavelength						() const	{ return m_wavelength;}

		static bool										GetQuadratureInterpParamsFromPointsAndLook	(	double* r0,						double* r1,
																										double* t0, double* t1, double* rt, 
																										const HELIODETIC_POINT&			startpoint,
																										const HELIODETIC_POINT&			endpoint,
																										const HELIODETIC_UNITVECTOR&	look );

	public:
		virtual bool									TraceRay_NewMethod					() = 0;

	public:














};

