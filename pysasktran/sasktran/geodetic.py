import sasktranif.sasktranif as skif
from sasktran.exceptions import wrap_skif_functionfail


class Geodetic(object):
    """
    Class which handles the representation of the Earth geometry.  The object is initialized through one of the
    `from` methods, and then after initialization various parameters are available.

    Parameters
    ----------
    geoid_name: str, optional

    Examples
    --------
    >>> import sasktran as sk
    >>> geodetic = sk.Geodetic()
    >>> geodetic.from_lat_lon_alt(latitude=60, longitude=30, altitude=0)
    >>> print(geodetic)
    ISKGeodetic: IAU 1976
     Latitude: 60.0, Longitude: 30.0, Altitude: 0.0
    >>> print(geodetic.local_up)
    [ 0.4330127  0.25       0.8660254]
    >>> print(geodetic.location)
    [ 2768775.09176071  1598553.04455358  5500479.72571591]
    """
    @wrap_skif_functionfail
    def __init__(self, geoid_name=None):
        self._iskgeodetic = skif.ISKGeodetic()

        if geoid_name is not None:
            self._iskgeodetic.SetProperty(geoid_name, 0)
            self._geoid_name = geoid_name
        else:
            self._geoid_name = 'IAU 1976'

    def __repr__(self):
        return "ISKGeodetic: {}\n Latitude: {}, Longitude: {}, Altitude: {}".format(self._geoid_name, self.latitude,
                                                                                    self.longitude, self.altitude)

    def __getstate__(self):
        state = {k: v for k, v in self.__dict__.items() if k != '_iskgeodetic'}
        return state, self.location

    def __setstate__(self, state):
        self.__dict__.update(state[0])
        self._iskgeodetic = skif.ISKGeodetic()
        if self._geoid_name != 'IAU 1976':
            self._iskgeodetic.SetProperty(self._geoid_name, 0)
        self.from_xyz(state[1])

    @wrap_skif_functionfail
    def from_tangent_altitude(self, altitude, observer, boresight):
        """
        Initialized the Geodetic from a specified tangent altitude, obsever location, and bore sight plane.

        Parameters
        ----------
        altitude : float
            Tangent altitude in meters
        observer : np.ndarray
            Three element array containing the obsever position in geocentric coordinates
        boresight : np.ndarray
            Three element array containing a normalized look vector that is within the bore sight plane.

        Returns
        -------
        np.ndarray
            Three element array containing the normalized look vector to the tangent point.

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> look = geodetic.from_tangent_altitude(15322, [ 3.676013154788849600e+005, 1.009976313640051500e+006,\
                                                          -6.871601202127538600e+006], [0, 0, 1])
        >>> print(look)
        [ 0.28880576  0.7934873   0.535695  ]
        >>> print(geodetic)
        ISKGeodetic: IAU 1976
         Latitude: -57.608525446153806, Longitude: 70.00000000000001, Altitude: 15321.971935943235
        """
        return self._iskgeodetic.SetLocationFromTangentAltitude(altitude, observer, boresight)[1]

    @wrap_skif_functionfail
    def from_tangent_point(self, observer, look_vector):
        """
        Initializes  the Geodetic by calculating the tangent point from an observer position and look vector

        Parameters
        ----------
        observer : np.ndarray
            Three element array containing the observer position in geocentric coordinates
        look_vector : np.ndarray
            Three element array containing a normalized look vector

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_tangent_point([ 3.676013154788849600e+005, 1.009976313640051500e+006,\
                                         -6.871601202127538600e+006], [ 2.884568631765662100e-001,\
                                          7.925287180643269000e-001,  5.372996083468238900e-001])
        >>> print(geodetic)
        ISKGeodetic: IAU 1976
         Latitude: -57.499726734289936, Longitude: 70.0, Altitude: 10000.000072686584
        """
        self._iskgeodetic.SetLocationFromTangentPoint(observer, look_vector)

    @wrap_skif_functionfail
    def from_lat_lon_alt(self, latitude, longitude, altitude):
        """
        Initializes the Geodetic based on a specifiec latitude, longitude, and altitude.

        Parameters
        ----------
        latitude : float
            Latitude in degrees (-90 to 90)
        longitude : float
            Longitude in degrees (0 to 360 or -180 to 180)
        altitude : float
            Altitude above the geoid in metres

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_lat_lon_alt(latitude=-15, longitude=-20, altitude=7342)
        >>> print(geodetic)
        ISKGeodetic: IAU 1976
         Latitude: -15.0, Longitude: 340.0, Altitude: 7342.0
        """
        self._iskgeodetic.SetLocationLatLonAlt(latitude, longitude, altitude)

    @wrap_skif_functionfail
    def from_xyz(self, location):
        """
        Initializes the Geodetic from a geocentric location

        Parameters
        ----------
        location : np.ndarray
            Three element vector containing a location in geocentric coordinates

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic)
        ISKGeodetic: IAU 1976
         Latitude: -15.000000000000009, Longitude: 340.00000000000006, Altitude: 7341.999999995809
        """
        self._iskgeodetic.SetLocationXYZ(location)

    @property
    @wrap_skif_functionfail
    def altitude(self):
        """
        float: Altitude in meters

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic.altitude)
        7341.999999995809
        """
        return self._iskgeodetic.GetAlt()

    @property
    @wrap_skif_functionfail
    def latitude(self):
        """
        float: Latitude in degrees in the range (-90 to 90)

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic.latitude)
        -15.000000000000009
        """
        return self._iskgeodetic.GetLatitude()

    @property
    @wrap_skif_functionfail
    def local_south(self):
        """
        np.ndarray : Three element vector unit vector tangent to the Earth pointing south

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic.local_south)
        [-0.24321035  0.08852133 -0.96592583]
        """
        return self._iskgeodetic.GetLocalSouth()

    @property
    @wrap_skif_functionfail
    def local_up(self):
        """
        np.ndarray : Three element unit vector normal to the surface of the Earth pointing up

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic.local_up)
        [ 0.90767337 -0.33036609 -0.25881905]
        """
        return self._iskgeodetic.GetLocalUp()

    @property
    @wrap_skif_functionfail
    def local_west(self):
        """
        np.ndarray : Three element unit vector tangent to the surface of the Earth pointing west

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic.local_west)
        [-0.34202014 -0.93969262  0.        ]
        """
        return self._iskgeodetic.GetLocalWest()

    @property
    @wrap_skif_functionfail
    def longitude(self):
        """
        float : Longitude in degrees in the range (0 to 360)

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic.longitude)
        340.00000000000006
        """
        return self._iskgeodetic.GetLongitude()

    @property
    @wrap_skif_functionfail
    def location(self):
        """
        np.ndarray : Three element vector representing the position in geocentric coordinates

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_lat_lon_alt(latitude=-15, longitude=-20, altitude=7342)
        >>> print(geodetic.location)
        [ 5797230.47518212 -2110019.3341472  -1642001.16317228]
        """
        return self._iskgeodetic.GetLocationXYZ()

    def altitude_intercepts(self, altitude, observer, look_vector):
        """
        Calculate the two intersections of a line of sight and an altitude.

        Parameters
        ----------
        altitude : float
            Altitude in meters.
        observer : np.ndarray
            Three element array containing the obsever position in geocentric coordinates.
        look_vector : np.ndarray
            Three element array containing a normalized look vector.

        Returns
        -------
        np.ndarray
            Three element array containing the first (entering) intercept in geocentric coordinates.
        np.ndarray
            Three element array containing the second (exiting) intercept in geocentric coordinates.

        Examples
        --------
        >>> import sasktran as sk
        >>> import numpy as np
        >>> geodetic = sk.Geodetic()
        >>> look = geodetic.from_tangent_altitude(15322, [3.676013154788849600e+005, 1.009976313640051500e+006, \
                                                      -6.871601202127538600e+006], [0, 0, 1])
        >>> obs = geodetic.location
        >>> intercept1, intercept2 = geodetic.altitude_intercepts(16000, obs, look)
        >>> print(np.array_str(intercept1, precision=3))
        [ 1147300.6    3152182.491 -5425366.163]
        >>> print(np.array_str(intercept2, precision=3))
        [ 1201097.087  3299987.124 -5325581.071]
        """

        try:
            ok, int1, int2 = self._iskgeodetic.GetAltitudeIntercepts(altitude, observer, look_vector)
            return int1, int2
        except:
            return None, None

    @property
    @wrap_skif_functionfail
    def osculating_spheroid_center(self):
        """
        np.ndarray : Three element vector representing the center of the osculating spheroid, relative to the center
        of the geoid, in meters.

        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic.osculating_spheroid_center)
        [ 36183.74814049 -13169.80728732    735.81389216]
        """
        return self._iskgeodetic.GetOsculatingSpheroidCenter()

    @property
    @wrap_skif_functionfail
    def osculating_spheroid_radius(self):
        """
        np.ndarray : Radius of the osculating spheroid in meters.
        Examples
        --------
        >>> import sasktran as sk
        >>> geodetic = sk.Geodetic()
        >>> geodetic.from_xyz([ 5797230.47518212, -2110019.3341472, -1642001.16317228])
        >>> print(geodetic.osculating_spheroid_radius)
        6339703.29902925
        """
        return self._iskgeodetic.GetOsculatingSpheroidRadius()
