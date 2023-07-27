#pragma once

#include <sasktran2/internal_common.h>

namespace sasktran2::coordinates {
#define RE_GEOCENTRIC 6371200.0		//!< Radius of geocentric earth in metres.


    class Geodetic
    {
    public:
        enum GEOID_MODEL { GEOID_SPHERE=0, IAU1976 = 1, GRS80 = 2, MERIT83=3, WGS84=4};	//!< Enumeration for the different oblate spheroid models

        static const double m_geoid_Re[5];  //= {6378140.0, 6378137.0,  6378137.0, 6378137.0};			//!< Storage for the difefren models Re term
        static const double m_geoid_F [5];  //= {298.257,   298.257222,	298.257,   298.257223563};		//! Storage for the difefrent models flattening term.

    private:
        Eigen::Vector3d 	m_location;							//!< Geographic Location of observer in Geocentric coords (X,Y,Z) in metres.
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
        Geodetic					();
        Geodetic					( GEOID_MODEL id);											//!< Construct using a specific oblate spheroid geoid
        Geodetic					( double latitude, double longitude,  double Height=0.0  );	//!< Construct and set users location from specified geodetic coordinate s
        void			SelectGeoid					( GEOID_MODEL id );											//!< Select a different pre-defined geoid for calculations (current location is confused)
        void			SetGeoid					( double newReA, double newRef);							//!< Select a different user defined geoid for calculations (current location is confused)
        void   			FromGeodetic				( double latitude, double longitude,  double Height=0.0 );	//!< Set the current location using the specified geodetic coordinates
        void   			FromGeocentric				( double latitude, double longitude,  double Height=0.0 );	//!< Set the current location using the specified geocentric coordinates (this is a rarely used function)
        void   			FromGeocentricVector		( const Eigen::Vector3d &Place );									//!< Set the current location from the specified geocentric X,Y,Z vector. (all in meters).
        Eigen::Vector3d		FromTangentPointLocation	( const Eigen::Vector3d& r, const Eigen::Vector3d& lookv );								//!< Set the current location from the implied tangent point
        bool			FromTangentAltitude			( double h, const Eigen::Vector3d& spacecraftlocation, const Eigen::Vector3d& boresightplane, Eigen::Vector3d* requiredlookvector );	//!< Set the current location from the tangent point at the specified height. Also return the \e look \e vector necessary to do this.
        bool			GetShellHeightLocation		( double H, const Eigen::Vector3d& observerposition, const Eigen::Vector3d& look, Eigen::Vector3d* entrypoint, Eigen::Vector3d* exitpoint, Eigen::Vector3d* lostangentpoint = NULL, double tangenth = 0.0);	//!< Calculate where a look direction enters and exits  a specific shell height.
        void   			GetGeodeticWestSouthUp		( Eigen::Vector3d *West, Eigen::Vector3d *South, Eigen::Vector3d *Up) const;			//!< Get the topocentric unit vectors at the current location
        void			SetUseExactConversion		( bool useexact)	{ m_useexact = useexact;}				//!< Use the exact conversion between geodetic and geocentric (default is TRUE)
        bool			IsPureSphere				( ) const			{ return m_Ref == 0.0;}			//!< Return nxTRUE if we are using a pure sphere representation
        const Eigen::Vector3d& Location					( ) const			{ return m_location;}			//!< Get the current geocentric location (X,Y,Z in meters)
        double			GeodeticLongitude			( ) const			{ return m_geodeticlongitude;}	//!< Get the current geodetic longitude (0-360) degrees, positive of Greenwich
        double			GeodeticLatitude			( ) const			{ return m_geodeticlatitude;}	//!< Get the current geodetic latitude in degrees
        double			Height						( ) const			{ return m_height;}				//!< Get the current height in meters
        double			SemiMajorAxis				( ) const			{ return m_ReA;}				//!< Get the current semi-major axis of the oblate spheroid
        double			SemiMinorAxis				( ) const			{ return m_ReA*(1-m_Ref);}		//!< Get the current semi-minor axis of the oblate spheroid
        bool			GetOsculatingSpheroid		( double* Radius, Eigen::Vector3d* offset );				//!< Get the sphere that best fits the oblate spheroid at the current point
        void			SetTrueSphere				( double radius );									//!< Set the current geoid to use a pure sphere.
    };

}