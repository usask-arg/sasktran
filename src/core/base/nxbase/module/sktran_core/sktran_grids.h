class HELIODETIC_VECTOR;
class HELIODETIC_POINT;
class SKTRAN_CoordinateTransform_V2;

#include "geodetic_instant.h"


/*-----------------------------------------------------------------------------
 *					HELIODETIC_UNITVECTOR		2009-1-25*/
/** A class representing a unit vector.  Unit vectors always have unit
 *	magnitude
 **/
/*---------------------------------------------------------------------------*/

class HELIODETIC_UNITVECTOR
{
	private:
		double			m_data[3];

	public:
									HELIODETIC_UNITVECTOR()									{ m_data[0] = -99999; m_data[1] = -99999; m_data[2] = -99999.0;}

		void						SetCoords( double x, double y, double z )				{	m_data[0] = x; m_data[1] = y; m_data[2] = z;
																								#if defined(NXDEBUG)
																									double sum = (m_data[0]*m_data[0] + m_data[1]*m_data[1] + m_data[2]*m_data[2]);
																									NXASSERT((  ((sum > 0.99999) && (sum < 1.00001)) || !NXFINITE(sum)  ));
																								#endif
																							}
		void						Clear		()											{ m_data[0] = 0; m_data[1] = 0; m_data[2] = 0;}
		HELIODETIC_UNITVECTOR&		FromVector	(const HELIODETIC_VECTOR &v); 
		double						X			()									const	{ return m_data[0];}
		double						Y			()									const	{ return m_data[1];}
		double						Z			()								    const	{ return m_data[2];}
		double 						operator&	( const HELIODETIC_UNITVECTOR &v2 ) const	{ return (m_data[0]*v2.m_data[0]) + (m_data[1]*v2.m_data[1]) + (m_data[2]*v2.m_data[2]);}
		HELIODETIC_UNITVECTOR& 		Negate		()											{ m_data[0] = -m_data[0]; m_data[1] = -m_data[1]; m_data[2] = -m_data[2]; return *this;}
//		const double*	DataPtr		() const									{ return m_data;}
};

/*-----------------------------------------------------------------------------
 *					HELIODETIC_VECTOR		2009-1-25*/
/** A vector, typically used to describe a position. The heliodetic
 *	coordinate system places the sun on the Z axis
 **/
/*---------------------------------------------------------------------------*/

class HELIODETIC_VECTOR
{
	private:
		double					m_data[3];

	public:
								HELIODETIC_VECTOR()								{ m_data[0] = -99999; m_data[1] = -99999; m_data[2] = -99999.0;}							
								HELIODETIC_VECTOR(const HELIODETIC_UNITVECTOR& unit, const double& magnitude) {this->SetCoords(unit, magnitude);}
		void					SetCoords	( double x, double y, double z );
		void					SetCoords	( const HELIODETIC_UNITVECTOR& unit, double magnitude );
		void					SetCoords	( const HELIODETIC_UNITVECTOR& unit); 
		double					X			() const							{ return m_data[0];}
		double					Y			() const							{ return m_data[1];}
		double					Z			() const							{ return m_data[2];}
		HELIODETIC_VECTOR&		operator+=	( const HELIODETIC_VECTOR &v2 )		{ m_data[0] += v2.m_data[0];    m_data[1] += v2.m_data[1];    m_data[2] += v2.m_data[2]; return *this;}
		HELIODETIC_VECTOR&		operator*=	( double f )						{ m_data[0] *= f;               m_data[1] *= f;               m_data[2] *= f;            return *this;}
		HELIODETIC_VECTOR&		operator-=	( const HELIODETIC_VECTOR &v2 )		{ m_data[0] -= v2.m_data[0];    m_data[1] -= v2.m_data[1];    m_data[2] -= v2.m_data[2]; return *this;}
		double 					Magnitude	() const							{ return sqrt(m_data[0]*m_data[0] +m_data[1]*m_data[1] + m_data[2]*m_data[2]);}																	//!< Get the magnitude of this vector
//		const double*			DataPtr		() const							{ return m_data;}
		HELIODETIC_VECTOR 		operator*	( double f ) const					  { HELIODETIC_VECTOR  v; v.SetCoords( m_data[0]*f,             m_data[1]*f,            m_data[2]*f ); return v;}
		HELIODETIC_VECTOR 		operator+	( const HELIODETIC_VECTOR &v2 ) const { HELIODETIC_VECTOR  v; v.SetCoords( m_data[0]+v2.m_data[0],  m_data[1]+v2.m_data[1], m_data[2]+v2.m_data[2] ); return v;}
		HELIODETIC_VECTOR 		operator-	( const HELIODETIC_VECTOR &v2 ) const { HELIODETIC_VECTOR  v; v.SetCoords( m_data[0]-v2.m_data[0],  m_data[1]-v2.m_data[1], m_data[2]-v2.m_data[2] ); return v;}
		double 					operator&	( const HELIODETIC_UNITVECTOR &v2 ) const	{ return (X()*v2.X()) + (Y()*v2.Y()) + (Z()*v2.Z());}
		HELIODETIC_UNITVECTOR	UnitVector	() const;
};


/*-----------------------------------------------------------------------------
 *					HELIODETIC_POINT		2009-1-25*/
/** A vector described as a unit vector and a magnitude. Useful
 *	when the unit vector is frequently required.
 **/
/*---------------------------------------------------------------------------*/

class HELIODETIC_POINT
{
	private:
		HELIODETIC_UNITVECTOR	m_direction;			// Unit vector, Sun is along Z axis
		double					m_radius;				// Radius from centre of Osculating sphere in meters
		double					m_heightm;				// Height in meters from surface of osculating sphere
		#if defined(NXDEBUG)
		nxVector				m_geo;					// The vector in geographic coords.
		#endif

	public:
										HELIODETIC_POINT()				{ m_radius = -99999; m_heightm = -99999;}
		void							Clear							();
		void							Initialize						( const HELIODETIC_UNITVECTOR& v, double magnitude, const SKTRAN_CoordinateTransform_V2* coords);
		void							InitializeFromRaw               ( const HELIODETIC_UNITVECTOR& v, double magnitude, double heightm ) { m_direction = v; m_radius = magnitude; m_heightm = heightm; }
		double							CosSZA							() const		{ return m_direction.Z();}
		double							LongitudeX						() const		{ return m_direction.X();}
		double							LongitudeY						() const		{ return m_direction.Y();}
		double							Altitude						() const		{ return m_heightm;}
		double							Radius							() const		{ return m_radius;}
		const HELIODETIC_UNITVECTOR&	LocalZenith						() const		{ return m_direction;}
		double							CosZenithAngle					( const HELIODETIC_UNITVECTOR& look ) const;
		bool							LocalUnitVectors				( HELIODETIC_UNITVECTOR* localunitvectors, size_t numv)  const;
		HELIODETIC_UNITVECTOR			TransformToLocalZenithCoords	( const HELIODETIC_UNITVECTOR& v, const HELIODETIC_UNITVECTOR* localunitvectors ) const;
		HELIODETIC_VECTOR				Vector							() const;
		HELIODETIC_UNITVECTOR           UnitVector                      () const;
		void							FromVector						( const HELIODETIC_VECTOR&	v, const SKTRAN_CoordinateTransform_V2* coords );
};



/*-----------------------------------------------------------------------------
 *					HELIODETIC_BASIS		 2014- 11- 6*/
/** Unit vectors to define an orthogonal coordinate system. The X unit vector 
 *  is always the direction of propagation (opposite to an observer's look 
 *  direction. 
 *  Under the built-in construction, the Y unit vector is made by rotating
 *  X through an angle pi/2 in the Sun direction (direction TO the Sun, not
 *  direction of propagation of sunlight), and Z=cross(X,Y). 

 vector of the ray while the Y and Z
 *	directions are perpendicular to the ray direction.
 **/
/*---------------------------------------------------------------------------*/

struct HELIODETIC_BASIS
{
	public:
		HELIODETIC_UNITVECTOR x;
		HELIODETIC_UNITVECTOR y;
		HELIODETIC_UNITVECTOR z;

	private:
		bool ProduceBasis ( const HELIODETIC_UNITVECTOR& pt, const HELIODETIC_UNITVECTOR& look );

	public:
												HELIODETIC_BASIS	()	{}
												HELIODETIC_BASIS	( const HELIODETIC_POINT&  pt, const HELIODETIC_UNITVECTOR& look )  { ProduceBasis( pt, look);}
												HELIODETIC_BASIS	( const HELIODETIC_VECTOR& pt, const HELIODETIC_UNITVECTOR& look )  { ProduceBasis( pt, look);}
		bool									ProduceBasis		( const HELIODETIC_POINT&  pt, const HELIODETIC_UNITVECTOR& look );
		bool									ProduceBasis		( const HELIODETIC_VECTOR& pt, const HELIODETIC_UNITVECTOR& look );
		const HELIODETIC_UNITVECTOR&			X() const			{ return x;}
		const HELIODETIC_UNITVECTOR&			Y() const			{ return y;}
		const HELIODETIC_UNITVECTOR&			Z() const			{ return z;}
};


/*-----------------------------------------------------------------------------
 *					GEOGRAPHIC_BASIS		 2015- 11-25*/
/** A very simple container to hold the information contained in
 *  HELIODETIC_BASIS, but in geographic coordinates. 
 *  
 **/
/*---------------------------------------------------------------------------*/
struct GEOGRAPHIC_BASIS
{
	public:
		nxVector x, y, z;

	public:
		GEOGRAPHIC_BASIS ( ) {;}
		GEOGRAPHIC_BASIS ( const nxVector& x, const nxVector& y, const nxVector& z ) : x(x), y(y), z(z) {;}
		const nxVector& X() const { return x;}
		const nxVector& Y() const { return y;}
		const nxVector& Z() const { return z;}
};

//class SKTRAN_GridSpecifications_V2;


/*-----------------------------------------------------------------------------
 *					class SKTRAN_CoordinateTransform_V2							2007-11-9*/
/** @ingroup coords
 *	The coordinate system used for the Sasktran calculation. Its main task is to
 *	create an osculating sphere at a reference point suitable for the calculation.
 *	The osculating sphere is a true sphere that has the same North/South
 *	curvature as the Earth at the reference point and has a zenith vector parallel
 *	to the geodetic zenith at the reference point. It significantly simplifies Sasktran
 *	by allowing a spherical approximation to be made with reasonable accuracy at any
 *	point on the Earth. The osculating sphere requires a translational offset from the
 *	true center of the Earth and adjusts the earth's radius to fit.
 *	The code provides methods to translate position vectors	from geographic/geodetic
 *	coordinates into the osculationg sphere system. Directional	vectors do not require
 *	any offset modification. The object also provides mappings between radius and altitude.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_CoordinateTransform_V2 : public nxUnknown
{
	private:
		static size_t				m_numinstances;			//!< The number of instances;
		nxGeodetic					m_truegeoid;			//!< The oblate spheroid used for the real world
		nxGeodetic					m_oscgeoid;				//!< The osculating spheroid used for the spherical approximation
		nxVector					m_centre;				//!< The offset of the center of the osculating earth in "true earth" X,Y,Z meters
		double						m_earthRadius;			//!< The radius of the osculating earth in meters
		nxVector					m_sun;					//!< The unit vector towards the Sun in geographic coordinates
		nxVector					m_heliounitvector[3];	//!< local unit vectors that express "Solar Coordinate" axes x', y' and z' in geographic coordinates;
		nxVector					m_referencepoint_unit;	//!< The unit vector towards the reference point in osculating sphere coordinates
		double						m_latitude;				//!< The reference geodetic latitude for this grid definition (used to define the osculating sphere)
		double						m_longitude;			//!< The reference geodetic longitude for this grid definition (used to define the osculating sphere)
		double						m_mjd;
		double						m_groundaltitude_meters;	//!< Altitude of the ground in meters
		double						m_toaaltitude_meters;		//!< Altitude of the top of the atmosphere in meters.

	private:
		bool						ConfigureGlobalTransform	();		//!< Configure the transform to convert global look directions into local zenith azimuth


	public:
		// ---- Constructors
									SKTRAN_CoordinateTransform_V2	();
		virtual					   ~SKTRAN_CoordinateTransform_V2	();
		static size_t				NumInstances						() { return m_numinstances;}
		bool						ConfigureCoordinates				( double latitude, double longitude, double mjd, const nxVector& sun );
		bool						ConfigureCoordinates				( const nxVector& observer, const nxVector& look, double mjd, const nxVector& sun );
		const nxVector&				GetSunUnitVector					() const							{ return m_sun; };
		double						AltitudeToRadius					( double alt_meters)		const	{ return (alt_meters + m_earthRadius);}
		double						RadiusToAltitude					( double radius_meters)		const	{ return floor((radius_meters - m_earthRadius)*1000.0 + 0.5)/1000.0;} // Round to the nearest micron
		const nxVector&				SunUnit								() const							{ return m_sun;}
		double						ReferencePtLatitude					() const							{ return m_latitude;}
		double						ReferencePtLongitude				() const							{ return m_longitude;}
		const nxVector&				ReferencePointUnitVector			() const							{ return m_referencepoint_unit;}
		double						ReferencePointMJD					() const							{ NXASSERT((m_mjd > 10000)); return m_mjd;}
		HELIODETIC_POINT			ReferencePoint						( double altitude_meters ) const;
		bool						ManuallySetOsculatingSphereRadius	( double radius );

		nxVector					TranslateGeoidToOsculatingSphere	( const nxVector& truecoord ) const	{ return (truecoord -m_centre); };
		nxVector					TranslateOsculatingSphereToGeoid	( const nxVector& osccoord ) const	{ return (osccoord +m_centre); };
		bool						HelioVectorToHelioPoint				( const HELIODETIC_VECTOR& vector, HELIODETIC_POINT* point) const;
		HELIODETIC_VECTOR			GeographicToHelio					( const nxVector& geographic ) const;
		HELIODETIC_UNITVECTOR		GeographicToHelioUnitVector			( const nxVector& geo ) const;
		nxVector					HelioVectorToGeographic				( const HELIODETIC_VECTOR& helio) const;
		nxVector					HelioUnitVectorToGeographic			( const HELIODETIC_UNITVECTOR& helio) const;
		GEODETIC_INSTANT			PointToGeodetic						( const HELIODETIC_POINT& location, double mjd ) const;
		const nxGeodetic&			TrueGeoid							() const 	{ return m_truegeoid;}			//!< The oblate spheroid used for the real world
		const nxGeodetic&			OsculatingGeoid						() const	{ return m_oscgeoid;}				//!< The osculating spheroid used for the spherical approximation
		bool						SetAtmosphereAltitudeBounds			( double groundalt_meters, double toaalt_meters);
		double						GroundAltitude						() const { return m_groundaltitude_meters;} 
		double						TOAAltitude							() const { return m_toaaltitude_meters;} 


};


/*-----------------------------------------------------------------------------
 *					typedef SKTRAN_GridIndex								*/
/** A typedef used to index grids. This type is synonomous with \c size_t.
 *	It was introduced to emphasize that a \c size_t parameter was being used
 *	to index a grid.
**/
/*---------------------------------------------------------------------------*/

typedef size_t SKTRAN_GridIndex;

/*-----------------------------------------------------------------------------
 *					class SKTRAN_GridDefBase_V2						2007-11-21*/
/** This class is the base class for all of the grids that store an array of 
*	values.
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_GridDefBase_V2 : public nxUnknown
{
	public:
		typedef std::vector<double>::iterator		iterator;
		typedef std::vector<double>::const_iterator	const_iterator;
		enum                        ENUM_GRIDSEARCHMODE    { GRIDSEARCH_NONUNIFORM,  GRIDSEARCH_UNIFORM };

	private:
		static size_t				m_numinstances;
//		double*						m_gridvalues;				//!< The values that define the grid (eg shell radii)
//		size_t						m_numgridpoints;			//!< The number of grid points in this grid definition
		std::vector<double>			m_gridvalues;

		ENUM_GRIDSEARCHMODE         m_gridsearchmode;
		double                      m_InvUniformDelta;

	private:
		void						ReleaseResources             ( );
		bool                        SetGridSearchMode_NonUniform ( );
		bool                        SetGridSearchMode_Uniform    ( );
		const_iterator              LowerBound                   ( const_iterator start, const_iterator finish, double value ) const;
		const_iterator              LowerBound_Uniform           ( const_iterator start, const_iterator finish, double value ) const;
		const_iterator              UpperBound                   ( const_iterator start, const_iterator finish, double value ) const;
		const_iterator              UpperBound_Uniform           ( const_iterator start, const_iterator finish, double value ) const;


	public:
		enum						ENUM_INTERPOLATIONMODE { OUTOFBOUND_EXTRAPOLATE, OUTOFBOUND_TRUNCATE, OUTOFBOUND_ZERO, OUTOFBOUND_ERROR};
									SKTRAN_GridDefBase_V2();
		virtual					   ~SKTRAN_GridDefBase_V2();
		static size_t				NumInstances()			{ return m_numinstances;}


		bool						AllocateGridArray		( size_t m_numpoints );
		virtual bool				CopyGridArray			( const double* source, size_t numpoints );
		virtual bool				CopyGridArray			( const std::vector<double>& source);
		bool						DeepCopy				( const SKTRAN_GridDefBase_V2& other	  );		//!< Copies the  grid array from another instance using the virtual function CopyGridArray

		virtual bool				FindBoundingIndices		( double x, ENUM_INTERPOLATIONMODE interpolationmode, SKTRAN_GridIndex* lowercell, double* lowerweight, SKTRAN_GridIndex* uppercell, double* upperweight) const;
		size_t						FindingBoundingIndices	( double x, ENUM_INTERPOLATIONMODE interpolationmode, SKTRAN_GridIndex* index,     double* weight,      size_t maxvertices) const;
		bool						IndexOfPointEqualOrAbove( double value, SKTRAN_GridIndex* index ) const;
		bool						IndexOfPointBelow		( double value, SKTRAN_GridIndex* index ) const;
		bool						IndexOfPointBelowOrEqual( double value, SKTRAN_GridIndex* index ) const;

		double&						AtVar					(SKTRAN_GridIndex index)		{ return m_gridvalues.at(index);} //NXASSERT((  index < m_numgridpoints)); return *(m_gridvalues + index);}
		const double&				At						(SKTRAN_GridIndex index) const	{ return m_gridvalues.at(index); } //NXASSERT((  index < m_numgridpoints)); return *(m_gridvalues + index);}
		iterator					begin					()								{ return m_gridvalues.begin();}
		const_iterator				begin					() const						{ return m_gridvalues.begin();}
		double*						FrontPtr				()								{ return (m_gridvalues.size() > 0) ? &m_gridvalues.front() : NULL;}
		const double*				FrontPtr				() const						{ return (m_gridvalues.size() > 0) ? &m_gridvalues.front() : NULL;}
		double						front					() const						{ return m_gridvalues.front();}
		double						back					() const						{ return m_gridvalues.back();}
		iterator					end						()								{ return m_gridvalues.end();}
		const_iterator				end						() const						{ return m_gridvalues.end();}
		void						erase					()								{ m_gridvalues.clear();}
		size_t						NumGridPoints			() const						{ return m_gridvalues.size();}
		size_t						size					() const						{ return m_gridvalues.size();}
		std::vector<double>&		ArrayVar				()								{ return m_gridvalues;}
		const std::vector<double>&	Array					() const						{ return m_gridvalues;}

		bool                        SetGridSearchMode       ( SKTRAN_GridDefBase_V2::ENUM_GRIDSEARCHMODE mode );
		ENUM_GRIDSEARCHMODE         GetGridSearchMode       ( ) const                       { return m_gridsearchmode; }
};
