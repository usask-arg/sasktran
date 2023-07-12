#if !defined(NXBASE_NXGEODETIC_H)
#define NXBASE_NXGEODETIC_H 1

/** \addtogroup Geodesy */
/**@{*/
/** \name Outline
    This section of the library is involved with different aspects of geodesy.
	The code provides an extensive set of subroutines for addressing the day to day
	problems encountered in geophysical atmospheric sciences. The main sections of the code are:
	\n
	- Geodetic coordinate conversions and various oblate spheroid calculations using class #nxGeodetic
	- Modified Julian Date calculations using class #nxTimeStamp
	- nxSatelliteBase and planetary orbits.
	\n
*/


/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

#define RE_GEOCENTRIC 6371200.0		//!< Radius of geocentric earth in metres.

/*----------------------------------------------------------------------------
 *			class nxGeodetic												 */
/** \par Overview
 *  A class for performing calculations related to an oblate
 *	spheroid Earth. It has code for:
 *	   - Coordinate conversion from geocentric to geodetic and vice-versa.
 *	   - Generation of topocentric (West, South, up) unit vectors,
 *	   - Determining tangent point locations
 *     - intersections of lines of sights withspecific shell heights.
 *	This class, depending upon its use, often appears as a proxy for a vector
 *	specifying the location of a point on the Earth and in this case is like a a replacement
 *	for nxVector. At other times the class is used as a coordinate transformation object that
 *	transforms from X,Y,Z coordinates to lat,long, altitude.
 *
 *	\par Oblate Spheroid
 *	The class represents the Earth using an oblate spheroid and defaults
 *	to the IAU1976 reference geoid. Others are available, like GRS80, MERIT83 and WGS84
 *	but there is very little to choose between the different models.  The important point
 *	is that all of the geoids are within ~100 m of the true (gravitational) shape
 *	of the Earth which is much better than the ~20 km errors if one uses a pure sphere.
 *	I recommend using the default IAU1976 model unless you have a genuine requirement
 *	to use another model.
 *
 *	\par Geocentric Coordinates
 *	This class uses two coordinate systems, a geocentric and a geodetic coordinate system to describe
 *	any point in 3-D space.  Both coordinates systems cover the whole of 3-D space including points
 *	below the surface.  The geocentric system is an orthogonal X,Y,Z cartesian system with the origin
 *	at the centre of the Earth. The Z axis is parallel to the
 *	Earth's spin axis (Z points towards the North Pole). X and Y are in the equatorial plane with X pointing
 *	at the Greenwich meridian. Y is the third axis in a right-handed orthogonal system (points toward ~China).
 *	The class has been coded so that all geocentric coordinates are stored in meters measured from the center
 *	of the Earth.
 *
 *	\par Geodetic Coordinates
 *	Geodetic coordinates are typically the coordinates used on geographic maps and consist of
 *	\e latitude, \e longitude and \e height. This class expresses latitude in degrees, +90 at
 *	the North pole to -90 at the South Pole.  Longitude is measured in degrees with 0 at the
 *	Greenwich meridian and increasing eastwards. Height is measured in meters above the oblate
 *	spheroid surface.  This coordinate system has been tested below the surface (even at the center)
 *	and has been shown to properly cover the whole of 3-D space. It is important to note that the
 *	geodetic latitude is the angle between the local vertical, defined by the tangent to the oblate
 *	spheroid surface and the equatorial plane. The local vertical is generally not parallel to the
 *	geocentric vector from the center of the Earth.
 *
 *	\par Applications for Non-Earth oblate spheroids
 *	It is possible to use the code for oblate spheroids that are much smaller or much greater than the Earth.
 *	However the user should note that some of the functions are approximate in nature and have been designed to
 *	converge when the change in position is less than some fraction of a meter (eg 0.01).  Clearly this convergence
 *	is useless if the oblate spheroid is only microns in size.  Much of the convergence criteria could be easily
 *	coded up and adjusted to follow the oblate spheroid semi-major axis but it has not yet been done.
 *
 *	\par References
 *	- GRS80 (IUGG, 1980), Page K13, Astronomical Almanac 1990
 *	- K.M. Borkowski,1989, Bulletin Geodesique. 63, 50-56
 */
/*---------------------------------------------------------------------------*/

class nxGeodetic
{
	public:
		enum GEOID_MODEL { GEOID_SPHERE=0, IAU1976 = 1, GRS80 = 2, MERIT83=3, WGS84=4};	//!< Enumeration for the different oblate spheroid models

	static const double m_geoid_Re[5];  //= {6378140.0, 6378137.0,  6378137.0, 6378137.0};			//!< Storage for the difefren models Re term
	static const double m_geoid_F [5];  //= {298.257,   298.257222,	298.257,   298.257223563};		//! Storage for the difefrent models flattening term.

   	private:
      	nxVector 	m_location;							//!< Geographic Location of observer in Geocentric coords (X,Y,Z) in metres.
      	double 		m_geodeticlongitude;				//!< The Geodetic/Geocentric Longitude of observer in degrees
      	double 		m_geodeticlatitude;					//!< The Geodetic latitude of the observer
      	double 		m_height;							//!< Height of observer above sea-level in metres.
		double		m_ReA;								//!< Equatorial radius of earth in metres
		double		m_Ref;								//!< Reciprocal Flattening
		bool		m_useexact;							//!< True if we use exact (but slower) geocentric conversion


	private:
		void			init();
		void			ExactGeocentricToGeodetic  ( double r, double z, double* fi,  double* h);
		void			IterateGeocentricToGeodetic( double r, double z );

   	public:
	     				nxGeodetic					();
						nxGeodetic					( GEOID_MODEL id);											//!< Construct using a specific oblate spheroid geoid
	     				nxGeodetic					( double latitude, double longitude,  double Height=0.0  );	//!< Construct and set users location from specified geodetic coordinate s
		void			SelectGeoid					( GEOID_MODEL id );											//!< Select a different pre-defined geoid for calculations (current location is confused)
		void			SetGeoid					( double newReA, double newRef);							//!< Select a different user defined geoid for calculations (current location is confused)
      	void   			FromGeodetic				( double latitude, double longitude,  double Height=0.0 );	//!< Set the current location using the specified geodetic coordinates
      	void   			FromGeocentric				( double latitude, double longitude,  double Height=0.0 );	//!< Set the current location using the specified geocentric coordinates (this is a rarely used function)
      	void   			FromGeocentricVector		( const nxVector &Place );									//!< Set the current location from the specified geocentric X,Y,Z vector. (all in meters).
		nxVector		FromTangentPointLocation	( const nxVector& r, const nxVector& lookv );								//!< Set the current location from the implied tangent point
		bool			FromTangentAltitude			( double h, const nxVector& spacecraftlocation, const nxVector& boresightplane, nxVector* requiredlookvector );	//!< Set the current location from the tangent point at the specified height. Also return the \e look \e vector necessary to do this.
		bool			GetShellHeightLocation		( double H, const nxVector& observerposition, const nxVector& look, nxVector* entrypoint, nxVector* exitpoint, nxVector* lostangentpoint = NULL, double tangenth = 0.0);	//!< Calculate where a look direction enters and exits  a specific shell height.
      	void   			GetGeodeticWestSouthUp		( nxVector *West, nxVector *South, nxVector *Up) const;			//!< Get the topocentric unit vectors at the current location
		void			SetUseExactConversion		( bool useexact)	{ m_useexact = useexact;}				//!< Use the exact conversion between geodetic and geocentric (default is TRUE)
		bool			IsPureSphere				( ) const			{ return m_Ref == 0.0;}			//!< Return nxTRUE if we are using a pure sphere representation
		const nxVector& Location					( ) const			{ return m_location;}			//!< Get the current geocentric location (X,Y,Z in meters)
		double			GeodeticLongitude			( ) const			{ return m_geodeticlongitude;}	//!< Get the current geodetic longitude (0-360) degrees, positive of Greenwich
		double			GeodeticLatitude			( ) const			{ return m_geodeticlatitude;}	//!< Get the current geodetic latitude in degrees
		double			Height						( ) const			{ return m_height;}				//!< Get the current height in meters
		double			SemiMajorAxis				( ) const			{ return m_ReA;}				//!< Get the current semi-major axis of the oblate spheroid
		double			SemiMinorAxis				( ) const			{ return m_ReA*(1-m_Ref);}		//!< Get the current semi-minor axis of the oblate spheroid
		bool			GetOsculatingSpheroid		( double* Radius, nxVector* offset );				//!< Get the sphere that best fits the oblate spheroid at the current point
		void			SetTrueSphere				( double radius );									//!< Set the current geoid to use a pure sphere.
};



/**@}*/
#endif


