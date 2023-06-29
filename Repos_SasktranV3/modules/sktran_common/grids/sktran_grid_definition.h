


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefAngular_V21		2007-11-21*/
/** @ingroup grids
 *	A base class defining a grid of angles **/
/*---------------------------------------------------------------------------*/


class SKTRAN_GridDefAngular_V21 : public SKTRAN_GridDefBase_V2
{
	public:
		virtual				   ~SKTRAN_GridDefAngular_V21()										{};
//		bool					ConfigureGeometry		( const double* angles, size_t numangles )	{ return CopyGridArray( angles, numangles );}
		std::vector<double>&	Angles					()										{ return ArrayVar();}		
		size_t					NumAngles				() const								{ return NumGridPoints(); }
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_GridDefDiffuseCosSZA_V2						2007-11-21*/
/** @ingroup grids
 *	A grid defining the location of diffuse profiles as cosine(SZA). 
 *	The points are stored as cos(SZA) in ascending order (ie -1 to +1). These
 *	values must be in ascending order as we use a binary search algorithm. 
 **/
/*---------------------------------------------------------------------------*/


class SKTRAN_GridDefCosSZA_V21 : public SKTRAN_GridDefAngular_V21
{
	public:
		virtual			   ~SKTRAN_GridDefCosSZA_V21		()							{ NXTRACE(("SKTRAN_GridDefCosSZA_V21::Destructor invoked\n"));};
		double				CosineSZA					( SKTRAN_GridIndex index ) const	{ return SKTRAN_GridDefBase_V2::At(index);}
		double				SZA							( SKTRAN_GridIndex index ) const	{ return acos(SKTRAN_GridDefBase_V2::At(index));}
		double				At							( SKTRAN_GridIndex index ) const	{ NXTRACE(("***** SKTRAN_GridDefCosSZA_V21::At, use SZA or CosineSZA for clarity\n")); return CosineSZA(index);}
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_GridDefDiffuseSLON_V2						2007-11-21*/
/** @ingroup grids
 *	A grid defining the location of diffuse profiles as the solar longitude in radians.
 *	The zero degree solar longitude is the meridian running from the sun (on the Z axis)
 *	and the osculating sphere refernce point. Positive longitude is towards
 *	"solar East". This array must match up with the elements of the correpsonding 
 *	SKTRAN_GridDefDiffuseCosSZA_V2. I.e. each element of both arrays specifies the SZA
 *	and solar longitude. This array cannot guarantee sorted order.
 *
 *	Note that we store the solar longitude in radians. We dont store the sin(SLON) as that
 *	is ambiguous for longitudes in the range 0-360 (or -180 to +180). This is important
 *	in overhead sun conditions (near the pole of the Sun) where we might have diffuse
 *	profiles scatetred over a large range of longitudes.
**/
/*---------------------------------------------------------------------------*/


class SKTRAN_GridDefSLON_V21 : public SKTRAN_GridDefAngular_V21
{
	public:
		virtual			   ~SKTRAN_GridDefSLON_V21		()							{ NXTRACE(("SKTRAN_GridDefSLON_V21::Destructor invoked\n"));};
		double				SinSLON						( SKTRAN_GridIndex index ) const	{ return sin(SKTRAN_GridDefBase_V2::At(index));}
		double				SLON						( SKTRAN_GridIndex index ) const	{ return SKTRAN_GridDefBase_V2::At(index);}
		double				At							( SKTRAN_GridIndex index ) const	{ NXTRACE(("***** SKTRAN_GridDefSLON_V21::At, use SLON or SinSLON for clarity\n")); return SLON(index);}
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_GridDefScatterAngle_V21		2007-11-21*/
/** @ingroup grids
*	A grid defining an array of scattering angles used when generating
*	interpolation tables fro scattering cross-sections. This table stores the
*	cosine of the scattering angle and spans the range -1 to +1. The table is
*	evenly spaced (in cosine space) and I have overloaded
*	FindBoundingIndices so it does the interpolation using direct indexing 
*	which I would expect to be a bit quicker than the generic
*	(upper_bound based) table search used by the base class.
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_GridDefScatterAngle_V21 : public SKTRAN_GridDefAngular_V21
{
	private:
		double				m_invdeltastep;
		double 				m_mincosssa;
		double 				m_maxcosssa;

	public:
							SKTRAN_GridDefScatterAngle_V21()	{ m_invdeltastep = 0.0;}
		virtual			   ~SKTRAN_GridDefScatterAngle_V21()	{}
		virtual bool		FindBoundingIndices		( double x, ENUM_INTERPOLATIONMODE outrange, SKTRAN_GridIndex* lowercell, double* lowerweight, SKTRAN_GridIndex* uppercell, double* upperweight) const;
		virtual bool		CopyGridArray			(  const double* source, size_t numpoints );	//!< Dont call this function, call Configure.
		bool				Configure				( double degrees_resolution, double minssa=0.0, double maxssa=180.0);
//		double				CosineScatter			( SKTRAN_GridIndex index ) const		{ return SKTRAN_GridDefBase_V2::At(index);}
//		double				ScatterAngle			( SKTRAN_GridIndex index ) const 		{ return acos(At(index));}
//		double				At						( SKTRAN_GridIndex index )	const		{ NXTRACE(("***** SKTRAN_GridDefScatterAngle_V21::At, use ScatterAngle or CosineScatter for clarity\n")); return CosineScatter(index);}
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_GridDefRayTracingShells_V21			2008-1-9*/
/**	@ingroup grids
 *	The ray tracing grid is used to define altitudes for tracing ray through
 *	the atmosphere. The intersection of a ray with this set of altitudes defines
 *	a set of infinitely thin shells while the distance between the intersection
 *	points defines a set of cells. Each ray that hits the atmosphere always has
 *	one less cell than shells. The radiative transfer integrals are
 *	evaluated for each cell using numerical quadrature.  The first/lowest altitude
 *	in this grid defines "the ground" and does not have to be located at an
 *	altitude of 0 meters.
 *
 *	\par Round off issues.
 *	Sasktran does attempt to avoid undershooting and over shooting altitude bins by
 *	"microns" due to floting point round off errors. Sasktran will normally round off
 *	ray tracing calculations to the nearest 1 mm. For this reason it is strongly
 *	recommended that altitudes be specified as proper integers only
 *	(formatted as floating point). We cannot see any reason at this time why any
 *	application of Sasktran would need to use fractions of a meter.
 *
 *	\par Default Value
 *	The default version of the tracing shells is a uniform 1 km seperation from
 *	0 to 100 km altitude.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_GridDefRayTracingShells_V21 : public SKTRAN_GridDefBase_V2
{

	public:
									SKTRAN_GridDefRayTracingShells_V21();
		virtual					   ~SKTRAN_GridDefRayTracingShells_V21();
		bool						ConfigureHeights	( const double* shellAlts, size_t numshells );
		bool						ConfigureHeights	( const std::vector<double>& shellAlts );
		size_t						NumShells			() const							{ return NumGridPoints(); };
		const std::vector<double>&	ShellHeight			() const							{ return Array();}
		size_t						NumCells			() const;
		double						LowestShell			() const;
		double						HighestShell		() const;
		bool						IsGroundPoint		( double altitude ) const { return altitude <= LowestShell();}
		double						GroundAltitude		() const { return LowestShell();}

};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_GridDefDiffuseIncomingZenith_V21		2007-12-11*/
/** @ingroup grids
 *	A grid used by the diffuse points in the diffuse table to specify the
 *	zenith angle of rays coming into a point. Zero zenith angle corresponds
 *	to rays coming straight down from vertical. All angles are specified in radians
 *	in the range 0 to pi.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_GridDefDiffuseIncomingZenith_V21 : public SKTRAN_GridDefAngular_V21
{
	private:
		bool				m_isgroundpoint;

	public:
							SKTRAN_GridDefDiffuseIncomingZenith_V21	()	{ m_isgroundpoint = false;};
		virtual			   ~SKTRAN_GridDefDiffuseIncomingZenith_V21	()	{ }
		void				SetIsGroundPoint( bool isground)			{ m_isgroundpoint = isground;}
		bool				IsGroundPoint() const						{ return m_isgroundpoint;} 
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_GridDefDiffuseIncomingAzimuth_V21		2007-11-21*/
/** @ingroup grids
 *	A grid used by the diffuse points in the diffuse table to specify the
 *	azimuth angles of rays coming into to each diffuse point. Zero degrees azimuth
 *	species rays coming in from the sun direction and increases anti-clockwise when
 *	looking down from above. The angles are specified in radians in the
 *	range 0 to two pi.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_GridDefDiffuseIncomingAzimuth_V21 : public SKTRAN_GridDefAngular_V21
{

	public:
							SKTRAN_GridDefDiffuseIncomingAzimuth_V21()	{}
		virtual			   ~SKTRAN_GridDefDiffuseIncomingAzimuth_V21()	{}
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_GridDefDiffuseHeights_V21		2007-12-12*/
/** @ingroup grids
 *	A grid used to specify the altitudes (in meters) of the points in the
 *	the array of diffuse profiles stored in the diffuse profile table.
 *	This grid is often considerbaly coarser than the altitude grid.
 *
 *	The first point is used to define the ground and this should match
 *	the ground height specified in the ray tracing specifications  and
 *	solzar transmission tables.
 *
 *	The internal engine tables will use the first and second points in the array
 *	as ground points. The first point is the ground with only downward
 *	(hemispherical) points and is used for albedo calculations. The second
 *	point is the ground point infinitessimally above the ground which has
 *	both upward and downward points and is used for atmospheric diffuse interpolation.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_GridDefDiffuseHeights_V21 : public SKTRAN_GridDefBase_V2
{

	public:
										SKTRAN_GridDefDiffuseHeights_V21		()										{}
		virtual						   ~SKTRAN_GridDefDiffuseHeights_V21		()										{}
		bool							ConfigureAltitudes					( double* angles, size_t numangles )	{ return CopyGridArray( angles, numangles );}
		const std::vector<double>&		Altitudes							() const								{ return Array();}
		size_t							NumAltitudes						() const								{ return NumGridPoints(); }
		double							GroundAltitude						() const								{ return At(0);}
		bool							IsGroundPoint						( double h ) const						{ return (h <= At(0)); }
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_GridDefOpticalPropertiesRadii_V21		2007-12-12*/
/** @ingroup grids
 *	A grid used to specify the altitude in meters of the optical properties
 *	of the atmosphere. This grid must have the first radius located at the ground
 *	equal to the first radius specified in class #SKTRAN_GridDefRayTracingShells_V21.
 *	This grid is usually offset from the ray tracing grid by about half a shell but is of
 *	equal resolution.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_GridDefOpticalPropertiesRadii_V21 : public SKTRAN_GridDefBase_V2
{

	public:
										SKTRAN_GridDefOpticalPropertiesRadii_V21		()										{}
		virtual						   ~SKTRAN_GridDefOpticalPropertiesRadii_V21		()										{}
		bool							ConfigureAltitudes							( double* angles, size_t numangles )	{ return CopyGridArray( angles, numangles );}
		const std::vector<double>&		Altitudes									() const								{ return Array();}
		size_t							NumAltitudes								() const								{ return NumGridPoints(); }
};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_AlbedoBRDF_V2					2008-7-18*/
/**	@ingroup albedo
 *  A set of classes used to provide some albedo directional outbound 
 *	directional functionality. This implementation is currently constrained to 
 *	variations in the zenith angle.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_AlbedoBRDF_V2 : public nxUnknown
{

	public:
							SKTRAN_AlbedoBRDF_V2() { }
		virtual			   ~SKTRAN_AlbedoBRDF_V2() { }
		virtual double		ReflectanceFunction	( double coszenith) const = 0;
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_AlbedoBRDFOpticalLambertian_V2		2008-7-17*/
/**	@ingroup albedo
 *	Impelments a Lambertian albedo. 
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_AlbedoBRDFLambertian_V2 : public SKTRAN_AlbedoBRDF_V2
{

	public:
		virtual			   ~SKTRAN_AlbedoBRDFLambertian_V2		(){}
		virtual double		ReflectanceFunction	( double coszenith ) const { return 1.0/nxmath::Pi;}
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefWavelengths		2020-03-04*/
 /**	
  *	Wavelength Grid
  **/
/*---------------------------------------------------------------------------*/
class SKTRAN_GridDefWavelength : public SKTRAN_GridDefBase_V2
{

public:
									SKTRAN_GridDefWavelength() {}
	virtual						   ~SKTRAN_GridDefWavelength() {}
	bool							ConfigureWavelengths(double* wavelengths, size_t numwavelengths) { return CopyGridArray(wavelengths, numwavelengths); }
	const std::vector<double>&		Wavelengths() const { return Array(); }
	size_t							NumWavelengths() const { return NumGridPoints(); }
};
