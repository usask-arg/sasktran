
// #define REPLICATE_ACE
/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_ACEFTSProfile		2014-1-27*/
/** Class to store the Refractive index profile used by the ACE-FTS. I have
 *	tried to track down the source reference and it seems to be 
 *	Filippenko (1982, PASP, 94, 715) and the references to Edlen and others
 *	in the document. The reference is older and there are now better equations 
 *	in the literature for refractivity (eg Ciddor 1992).
 *
 *	This class hit one problem in implementation in that we want to get
 *	the refractive index of air when tracing a ray through the atmosphere. This
 *	creates a potential multi-thread issue as the code wants to call
 *	climatologies which are not thread-safe during ray tracing.
 *
 *	The solution performed in this class is to generate the profile when the
 *	profile location is set.
**/
/*---------------------------------------------------------------------------*/

class skRTRefractiveIndex_ACEFTSProfile
{
	private:
		std::shared_ptr<const SKTRAN_GridDefRayTracingShells_V21>	m_raytracingshells;				//!< The ray tracing shells used to generate the refractive index height profile
		skClimatology*								m_atmosphericstate;				//!< The atmospheric state used for ray tracing (needs pressure (Pascal)  and temperature (K))
		double										m_Q;							//!< m_Q used as an intermediate term in refractive index calculkation
		nxSpline									m_P;							//!< Spline used to interpolate log(refractivity) as a function of altitude
		nxSpline									m_T;							//!< Spline used to interpolate log(refractivity) as a function of altitude

	private:
		void					ReleaseResources					();
		bool					INTERP								( double h_meters, double* userR ) const;
		bool					SetWavenumber						( double NUS );										// Update the refractive index profile

	public:
								skRTRefractiveIndex_ACEFTSProfile	();
							   ~skRTRefractiveIndex_ACEFTSProfile	();
		bool					Initialize							( skClimatology* atmospheric_raytracing_state, std::shared_ptr<const SKTRAN_GridDefRayTracingShells_V21>	raytracingshells );
		bool					SetProfileLocationAndWavenumber		( const GEODETIC_INSTANT& point, double NUS);					// Calculate Refractivity height profile
		double					RefractiveIndex						( double h_meters ) const;							// Thread safe interpolation of refractive index table
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_RayRefracted_TrajectoryData		2014-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayRefracted_TrajectoryData
{
	public:
		struct pointentry
		{
			double	r;						// Radial location of this point
			double	angle;
			double	length;

		};
	private:
		std::vector<pointentry>							m_trajectorypoints;							//!< The radius of each point along the trajectory
		typedef std::vector<pointentry>::iterator		iterator;
		
	public:
		typedef std::vector<pointentry>::const_iterator	const_iterator;

	public:
									SKTRAN_RayRefracted_TrajectoryData			(){}
								   ~SKTRAN_RayRefracted_TrajectoryData			(){}
		bool						ReserveSpace								( size_t numelements);
		bool						PushBack									( double r, double angle, double pathlength);
		bool						MirrorPoints								();
		bool						SumAnglesAndLengths							();
		const_iterator				begin										() const { return m_trajectorypoints.begin();}
		const_iterator				end											() const { return m_trajectorypoints.end  ();}
		size_t						size										() const { return m_trajectorypoints.size ();}
		HELIODETIC_VECTOR			ConvertTo3DLocation							( const HELIODETIC_VECTOR& observer, const HELIODETIC_POINT& observerpt, const HELIODETIC_UNITVECTOR& look, const_iterator& iter );
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_RayRefracted_ThomPepSim				2014-1-28*/
/** A ray tracing class that implements a curved ray using the
 *	technique outlined in "Ray tracing in a refracting spherically symmetric
 *	atmosphere". Thompson, Pepin and Simon.  J. Opt Soc. Am., 72, 11, 1982.
 *
 *	This is the technique used by the ACE-FTS analysis and I have compared
 *	this code to the ACE-FTS code. The resultant path lengths agreed to
 *	within 1.5 microns (i.e. the same).
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayRefracted_ThomPepSim
{
	private:
		std::shared_ptr< const SKTRAN_GridDefRayTracingShells_V21>		m_raytracingshells;
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>			m_coordinates;	
		skRTRefractiveIndex_ACEFTSProfile				m_RI;
		GEODETIC_INSTANT								m_point;						//!< Reference Latitude and longitude and time of refractive index porofile
		double											m_zerokmradius;

	private:
		const SKTRAN_CoordinateTransform_V2*			CoordsPtr						() const 	{return m_coordinates.get();}
		double											AltitudeToRadius				( double alt)   const	 { return m_zerokmradius + alt;}
		double											RadiusToAltitude				( double radius) const	 { return radius - m_zerokmradius;}
		bool											ComputeApproximateTangentRadius	( double HOBS, double RSPEC, double NADANG,  double * RT) const ;
		bool											FindGlobalTangentPoint			( double *RT, double* minradius, bool* groundishit) const;
		double											RadialHeightOfCosTheta			( double k, double sintheta	) const;

		bool											TraceRayOutsideAtmosphere		( double				RSPEC,
																						  double                RTOTL,
																						  double                NADANG,
																						  double				RT,
																						  double				minradius,
																						  bool					groundishit,
																						  SKTRAN_RayRefracted_TrajectoryData*	trajectory) const;		//!< Slant Path Geometric Factors

		bool											TraceRayInZenith				(  double									robs,					//!< radius of the observer
																						   SKTRAN_RayRefracted_TrajectoryData*		trajectory) const;		//!< Object to store the ray trajectory parameters

		bool											TraceRayInNadir					(  double				robs,									//!< radius of the observer
																						   SKTRAN_RayRefracted_TrajectoryData*	trajectory) const;		//!< Object to store the ray trajectory parameters

		bool											TraceRayInsideAtmosphere		( double*				HMIN,			//!< MINIMUM HEIGHT FOR THIS GEOMETRY
																						  std::vector<double>*	SPIN,			//!< Slant Path In
																						  std::vector<double>*	SPOUT,			//!< Slant Path out
																						  std::vector<double>*	SP) const;		//!< Slant Path Geometric Factors

		double											IntegrateCurvedPathLengthOld	(  double RT, double REFTAN, double RJ, double RJ1, int    RINC) const;
		bool											IntegrateCurvedPathLength		(  double RT, double REFTAN, double RJ, double RJ1, size_t numsteps, double* L, double* PSI) const;


	public:
														SKTRAN_RayRefracted_ThomPepSim();
													   ~SKTRAN_RayRefracted_ThomPepSim();
		static double									EarthRadius						( double TANLAT );
		bool											Configure						( std::shared_ptr< const SKTRAN_GridDefRayTracingShells_V21> raytracingshells, 
																						  std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords);
		bool											SetRILocationAndWavenumber		( const GEODETIC_INSTANT& point, double NUS );
		bool											SetRIAtmosphericState			( skClimatology* atmosphericstate);
		bool											FindObserverAndLookFromTangentAltitude( double tanalt, double spacecraftalt, const nxVector& sun, nxVector* observer, nxVector* look);
		bool											REFRAC							( double								HOBS,							// Height of observer in meters
																						  double								ZANGLE,							// Zenith angle in degrees of the of look vector which is away from observer
																						  SKTRAN_RayRefracted_TrajectoryData*	trajectory) const;				// returns the trajectory points in the specialized container

};
