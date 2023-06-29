import sasktran as sk
import numpy as np
import logging
from sasktran.util import ListWithCallback
from sasktran.geodetic import Geodetic


class Geometry(object):
    """
    Class which represents the Geometry for a radiative transfer calculation.  Currently this is made up of
    a set of :py:class:`sasktran.LineOfSight` which define the viewing geometry, and a :py:attr:`sasktran.Geometry.sun`
    property which controls the solar position.  In the future this class may be extended to also control options
    such as the top of the atmosphere, and the 'ground' altitude.

    Note that it is not always necessary to explicitly set the sun position.  If no sun position is given then it
    is automatically calculated from the average `mjd` given in the lines of sight.

    This class can be thought of as the geometry interface to the radiative transfer model.  For simple
    problems, we can construct this class and modify the lines of sight directly, as is done in the example below.
    For more complicated problems it may be desirable to create a new class that inherits from this class, as is
    done in the :class:`VerticalImage` class.

    Examples
    --------
    >>> import sasktran as sk
    >>> geometry = sk.Geometry()
    >>> los1 = sk.LineOfSight(mjd=54832.5, observer=[3.6760131547888e+005, 1.0099763136400e+006, -6.871601202127e+006],\
                              look_vector=[2.884568631765662e-001, 7.925287180643269e-001,  5.372996083468238e-001]),
    >>> los2 = sk.LineOfSight(mjd=54832.5, observer=[3.6928085406796e+005, 1.0145908079886e+006, -6.870844156040e+006],\
                              look_vector=[2.88456863176566e-001, 7.925287180643269e-001,  5.372996083468238e-001])
    >>> geometry.lines_of_sight = [los1, los2]
    >>> print(geometry)
    Geometry object containing 2 lines of sight
    >>> geometry.sun = [0, 0, 1]
    >>> print(geometry)
    Geometry object containing 2 lines of sight
    Sun position forced to [0, 0, 1]
    """
    def __init__(self):
        self._los = ListWithCallback()
        self._los.add_modified_callback(self._set_has_changed)
        self._has_changed = False
        self._sun = None
        self._ref_point = None

    def __repr__(self):
        ret = "Geometry object containing {} lines of sight".format(len(self._los))
        if self._sun is not None:
            ret += '\nSun position forced to {}'.format(self._sun)
        return ret

    @property
    def lines_of_sight(self):
        """
        list : A list of :py:class:`sasktran.LineOfSight` in which to calculate radiances.

        Examples
        --------
        >>> import sasktran as sk
        >>> geometry = sk.Geometry()
        >>> arg_office = sk.Geodetic()  # line of sight terminator
        >>> arg_office.from_lat_lon_alt(52.131638, -106.633873, 0)
        >>> observer = sk.Geodetic()   # line of sight observer
        >>> observer.from_lat_lon_alt(0, -100, 35786000)
        >>> look = (arg_office.location - observer.location)
        >>> look = look / np.linalg.norm(look)  # don't forget to normalize the look vector
        >>> # in this example we'll only add one line of sight (done below) but any number could be added to this list
        >>> geometry.lines_of_sight = [ \
             sk.LineOfSight(look_vector=look, observer=observer.location, mjd=57906), \
        ]
        """
        return self._los

    @lines_of_sight.setter
    def lines_of_sight(self, value):
        self._los = ListWithCallback(value)
        self._has_changed = True

    @property
    def sun(self):
        """
        np.ndarray((3,)) : Unit vector pointing to the sun in geodetic coordinates. See
        :py:attr:`sasktran.Geodetic.location`.

        Examples
        --------
        >>> import sasktran as sk
        >>> geometry = sk.Geometry()
        >>> sun = sk.Geodetic()
        >>> sun.from_lat_lon_alt(22.26666667, -48.98333333, 0) # altitude (argument 3) doesnt matter here
        >>> geometry.sun = sun.local_up # a geocentric unit vector pointing at the sun
        """
        return self._sun

    @sun.setter
    def sun(self, value):
        self._sun = value
        self._has_changed = True

    @property
    def reference_point(self):
        """
        np.ndarray((4,)) : The reference point as a four element array. Format: [latitude, longitude,
                           surface_altitude, mjd]. Note that longitude must be a value in the range [0, 360).

        Notes
        -----
        Reference points semantics are engine specific. Refer to your engines specific documentation for details.
        """
        return self._ref_point

    @reference_point.setter
    def reference_point(self, value):
        self._ref_point = value
        self._has_changed = True

    @property
    def has_changed(self):
        return self._has_changed

    def _set_has_changed(self):
        self._has_changed = True

    def _added_to_engine(self):
        self._has_changed = False


class VerticalImage(Geometry):
    """
    A specialized :class:`Geometry` class which has various convenience methods to construct lines of sight
    corresponding to a vertical image of the atmosphere based on input solar zenith angles.

    Examples
    --------
    >>> from sasktran.geometry import VerticalImage
    >>> geometry = VerticalImage()
    >>> geometry.from_sza_saa(60, 60, 0, 0, [10, 20, 30, 40], 54372, 0)
    >>> print(geometry)
    Geometry object containing 4 lines of sight
    Sun position forced to [ 0.5        0.75       0.4330127]
    """
    def __init__(self):
        super().__init__()

        self.reftp = None
        self.reflook = None
        self.refmjd = None

    def from_sza_saa(self, sza, saa, lat, lon, tanalts_km, mjd, locallook, satalt_km=600, refalt_km=20):
        """
        Constructs the geometry from a given solar zenith angle and solar azimuth angle

        Parameters
        ----------
        sza : float
            Solar zenith angle in degrees of the tangent point at refalt_km.
        saa : float
            Solar azimuth angle in degrees of the tangent point at refalt_km.
        lat : float
            Latitude in degrees (-90 to 90)
        lon : float
            Longitude in degrees
        tanalts_km : np.ndarray
            Tangent altitudes in km
        mjd : float
            Modified Julian Date for the measurements
        locallook : float
            Angle in degrees defining the look direction on the surface of the Earth.  0 corresponds to true north,
            90 to East and so on.
        satalt_km : float, optional
            Altitude of the observers in km.  Default is 600 km.
        refalt_km : float, optional
            Altitude that the given solar zenith angle and solar azimuth angle are valid at.  Since the geometry
            represents a vertical image, the solar angles change for every line of sight.

        """
        geo = sk.Geodetic()

        geo.from_lat_lon_alt(lat, lon, refalt_km * 1000)

        self.reftp = geo.location
        west = geo.local_west
        south = geo.local_south
        up = geo.local_up

        sunh = -np.cos(np.deg2rad(saa)) * south - np.sin(np.deg2rad(saa)) * west
        self._sun = np.cos(np.deg2rad(sza)) * up + np.sin(np.deg2rad(sza)) * sunh

        self.reflook = -np.cos(np.deg2rad(locallook)) * south - np.sin(np.deg2rad(locallook)) * west

        a = (np.linalg.norm(self.reftp) + satalt_km * 1000) ** 2

        s = -np.dot(self.reftp, -self.reflook) + np.sqrt(
            np.dot(self.reftp, -self.reflook) ** 2 - np.linalg.norm(self.reftp) ** 2 + a)

        obs = self.reftp - s * self.reflook
        self.refobs = self.reftp - s * self.reflook
        self._los = []
        for alt in tanalts_km:
            look = geo.from_tangent_altitude(alt * 1000.0, obs, self.reflook)

            self._los.append(sk.LineOfSight(mjd, obs, look))
        self.refmjd = mjd
        self._calc_real_parameters()

    def from_sza_ssa(self, sza, ssa, lat, lon, tanalts_km, mjd, locallook, satalt_km=600, refalt_km=20):
        sza_r = np.deg2rad(sza)
        ssa_r = np.deg2rad(ssa)

        if np.cos(ssa_r) / np.sin(sza_r) > 1:
            logging.warning('Unphysical geometry, forcing scattering angle to be the maximum allowed, %f',
                            np.rad2deg(np.arccos(np.sin(sza_r))))
            saa = locallook
        else:
            saa = np.rad2deg(np.arccos(np.cos(ssa_r) / np.sin(sza_r))) + locallook

        self.from_sza_saa(sza, saa, lat, lon, tanalts_km, mjd, locallook, satalt_km, refalt_km)

    def _calc_real_parameters(self):
        geo = sk.Geodetic()
        geo.from_tangent_point(self.refobs, self.reflook)
        self.reftp = geo.location

        self.lon = geo.longitude
        self.lat = geo.latitude

        # SSA is sun dot look
        self.ssa = np.rad2deg(np.arccos(np.dot(self._sun, self.reflook)))
        # SZA is sun dot tp
        self.sza = np.rad2deg(np.arccos(np.dot(self._sun, self.reftp / np.linalg.norm(self.reftp))))
        # TODO: calc saa
        self.saa = None


class NadirGeometry(Geometry):
    def __init__(self):
        super().__init__()
        self._int_lat = None
        self._int_lon = None
        self._oza = None
        self._oaa = None
        self._sza = None
        self._saa = None
        self._osc_adjust = False

    def from_lat_lon(self, mjd, observer: Geodetic, lats, lons, elevations=0):
        """
        Construct lines of sight from arrays of points (latitude, longitude[, elevation]).

        Notes
        -----
        Arguments for this function are associated sequences (ie. the first observation is at (lats[0], lons[0], the
        second is at (lats[1], lons[1]) and so on). This means that if you have a geoegraphic grid of sample points,
        these arrays should be flattened before being passed to this function.

        Arguments are broadcast together using numpy's broadcast rules.

        Examples
        --------
        >>> geometry = sk.NadirGeometry()
        >>> # make the look vector from TEMPO to the ARG office
        >>> tempo = sk.Geodetic()
        >>> tempo.from_lat_lon_alt(0, -100, 35786000)
        >>> geometry.from_lat_lon(mjd=57906, observer=tempo, lats=52.131638, lons=-106.633873, elevations=0)

        Parameters
        ----------
        mjd : float, np.ndarray
            Modified Julian Data of observations. If ``mjd`` is a ``float`` then this mjd will be used for each line of
            sight. Else mjd's are associated element-wise.
        observer : sasktran.Geodetic or array of sasktran.Geodetic
            Observer locations.
        lats : float, np.ndarray
            Latitudes of observations. If this argument is a float then this latitude will be used for each line of
            sight, otherwise latitudes will be flattened and associated element wise.
        lons : float, np.ndarray
            Same as ``lats`` but for longitudes.
        elevations : float, np.ndarray
            Used to specify the ground elevation (in meters) for each sample point.
        """
        # broadcast arguments to sequence of sample points
        mjd = np.atleast_1d(mjd).astype(float)
        obs = np.atleast_1d(observer)
        lats = np.atleast_1d(lats).astype(float)
        lons = np.atleast_1d(lons).astype(float)
        elev = np.atleast_1d(elevations).astype(float)

        shape = np.broadcast(mjd, obs, lats, lons, elev).shape
        mjd = np.broadcast_to(mjd, shape).flatten()
        obs = np.broadcast_to(obs, shape).flatten()
        lats = np.broadcast_to(lats, shape).flatten()
        lons = np.broadcast_to(lons, shape).flatten()
        elev = np.broadcast_to(elev, shape).flatten()

        # calculate lines of sight
        self.lines_of_sight = []
        for latitude, longitude, observer, mj_date, elevation in zip(lats, lons, obs, mjd, elev):
            terminator = sk.Geodetic()
            terminator.from_lat_lon_alt(latitude, longitude, elevation)

            look = terminator.location - observer.location
            look = look / np.linalg.norm(look)

            self.lines_of_sight.append(
                sk.LineOfSight(mjd=mj_date, observer=observer.location, look_vector=look)
            )

        self._calc_intersection()
        self._calc_los_angles()
        self._calc_sun_angles()
        self._has_changed = True

    def from_fov(self, mjd: float, observer: Geodetic,
                 center: Geodetic, ntheta: int=1, nphi: int=1,
                 dtheta: float=0.0, dphi: float=0.0,
                 bearing: float=None):
        """
        Adds lines of sight by creating a rectangular grid of equally
        spaced angles surrounding a center line of sight. The center line of
        sight is  defined by the observer position and the center position.

        Parameters
        ----------
        observer : sk.Geodetic
            Position of the observer.
        center : sk.Geodetic
            The center of the bundle of lines of sight.
        ntheta : int
            The number of rows in the angle grid.
        nphi : int
            The number of columns in the angle grid.
        dtheta: float
            The spacing between columns in the angle grid.
        dphi : float
            The spacing between rows in the angle grid.
        bearing : float
            The bearing (degrees, CW from north) of the center column's
            intersection with the local horizontal at center. Defaults to None,
            which places the column in the plane containing the observer and
            the local vertical at center, or in the north/south direction
            if the los is vertical.
        """

        # get the center look vector
        look_center = center.location - observer.location
        look_center = look_center / np.linalg.norm(look_center)

        # get the axis perpendicular to the center column plane
        if bearing is None:
            if abs(np.dot(look_center, center.local_up) + 1) < 1e-6:
                # if the local up and the los are parallel (ie if the los is vertical)
                col_axis = center.local_west
            else:
                col_axis = np.cross(look_center, center.local_up)
        else:
            cosb = np.cos(np.deg2rad(bearing))
            sinb = np.sin(np.deg2rad(bearing))
            local_horizontal = -(cosb * center.local_south + sinb * center.local_west)
            col_axis = np.cross(look_center, local_horizontal)

        # get the LOSs of the center column
        theta_range = (ntheta - 1) * dtheta
        theta = np.linspace(-theta_range / 2, theta_range / 2, ntheta)
        col = np.zeros((ntheta, 3))
        if ntheta == 1:
            col[0, :] = look_center
        else:
            for i in range(ntheta):
                R = self._rotation_matrix(col_axis, theta[i])
                col[i, :] = np.dot(R, look_center)

        # fan the center column out into a grid
        phi_range = (nphi - 1) * dphi
        phi = np.linspace(-phi_range / 2, phi_range / 2, nphi)
        row_axis = np.cross(look_center, col_axis)
        looks = np.zeros((ntheta, nphi, 3))
        if nphi == 1:
            looks[:, 0, :] = col
        else:
            for j in range(nphi):
                R = self._rotation_matrix(row_axis, phi[j])
                looks[:, j, :] = np.dot(col, R)

        self._los = []
        for i in range(ntheta):
            for j in range(nphi):
                self._los.append(sk.LineOfSight(mjd, observer.location, looks[i, j, :]))

        self._calc_intersection()
        self._calc_los_angles()
        self._calc_sun_angles()
        self._has_changed = True

    def from_zeniths_and_azimuths(self, solar_zenith: float, solar_azimuth: float,
                                  observer_mjds, observer_zeniths, observer_azimuths, observer_altitudes=600e3,
                                  reference_point=None):
        """
        Configures the geometry from arrays of zenith and azimuth angles. Only one solar zenith angle can be specified.

        Notes
        -----
        Observation parameters are forced to np.arrays and then broadcasted to a common shape. This means that any
        observation parameter can be a scalar (value used for all entries) or an array (each entry is individually
        specified).

        A common application is to specify a meshgrid of zenith and azimuth angles. In this case, simply flatten the
        arrays before passing them as np.arrays.

        Parameters
        ----------
        solar_zenith : float
            Solar zenith angle in degrees.
        solar_azimuth : float
            Solar azimuth angle (relative to North) in degrees.
        observer_mjds : float, np.array
            Observation MJDs.
        observer_zeniths : float, np.array
            Observation zenith angles in degrees.
        observer_azimuths : float, np.array
            Observation azimuth angles (relative to North) in degrees.
        observer_altitudes : float, np.array
            Observation surface altitude in meters.
        reference_point : np.array
            The reference point as a four element array (lat, lon, alt, mjd).
        """
        # fix solar_zenith of 0
        dither = 1e-2
        if np.abs(solar_zenith) < dither:
            solar_zenith = dither

        # update reference point if necessary
        if reference_point is None and self.reference_point is None:
            raise ValueError('The reference point has not been specified.')
        elif reference_point is not None:
            self.reference_point = reference_point

        # get the reference point
        reference_point = sk.Geodetic()
        reference_point.from_lat_lon_alt(*self.reference_point[:3])

        north = -reference_point.local_south
        east = -reference_point.local_west
        up = reference_point.local_up

        def horizontal(azimuth):
            return north * np.cos(np.deg2rad(azimuth)) + east * np.sin(np.deg2rad(azimuth))

        def unit_vector(zenith, azimuth):
            return up * np.cos(np.deg2rad(zenith)) + horizontal(azimuth) * np.sin(np.deg2rad(zenith))

        # set sun
        sun = unit_vector(solar_zenith, solar_azimuth)

        # make type consistent (i.e. at least 1d np.array)
        oza_deg = np.atleast_1d(observer_zeniths).astype(float).flatten()
        oaa_deg = np.atleast_1d(observer_azimuths).astype(float).flatten()
        o_mjd = np.atleast_1d(observer_mjds).astype(float).flatten()

        # fix parallel azimuths
        oaa_deg[np.where(np.abs(oaa_deg - solar_azimuth) < dither)] += dither

        # get broadcast size
        shape = np.broadcast(o_mjd, oza_deg, oaa_deg, observer_altitudes).shape
        o_mjd = np.broadcast_to(o_mjd, shape)
        oza_deg = np.broadcast_to(oza_deg, shape)
        oaa_deg = np.broadcast_to(oaa_deg, shape)
        observer_altitudes = np.broadcast_to(observer_altitudes, shape)

        # calculate lines of sight
        self.lines_of_sight = []
        for mjd, zenith, azimuth, obs_alt in zip(o_mjd, oza_deg, oaa_deg, observer_altitudes):
            observer = reference_point.location + obs_alt * unit_vector(zenith, azimuth)
            look = -unit_vector(zenith, azimuth)
            self.lines_of_sight.append(sk.LineOfSight(mjd, observer, look))
        self.sun = sun

    def from_zeniths_and_azimuth_difference(self, solar_zenith, observer_zeniths, azimuth_differences,
                                            solar_azimuth=180.0, mjd=58197.666, observer_alt=600e3,
                                            reference_point=(52.1332, 253.33, 0, 58197.666)):
        """
        Configures the geometry from arrays of solar zenith angles, observer zenith angles, and azimuth angle
        differences. If only one solar zenith angle is provided, from_zeniths_and_azimuths is called.

        If multiple solar zenith angles are provided, a north-facing (azimuth=0) sun is positioned over the reference
        point with the average solar zenith angle, ground locations to the south are chosen to match the given solar
        zenith angles, and observers and look vectors are calculated from these ground locations. Since ground
        intersection points will be away from the reference point, calling transform_los_to_preserve_ground_angles()
        is recommended.

        Notes
        -----
        Observation parameters are forced to np.arrays and then broadcasted to a common shape. This means that any
        observation parameter can be a scalar (value used for all entries) or an array (each entry is individually
        specified).

        A common application is to specify a meshgrid of zenith and azimuth angles. In this case, simply flatten the
        arrays before passing them as np.arrays.

        Parameters
        ----------
        solar_zenith : float
            Solar zenith angle in degrees.
        solar_azimuth : float
            Solar azimuth angle (relative to North) in degrees.
        observer_mjds : float, np.array
            Observation MJDs.
        observer_zeniths : float, np.array
            Observation zenith angles in degrees.
        observer_azimuths : float, np.array
            Observation azimuth angles (relative to North) in degrees.
        observer_altitudes : float, np.array
            Observation surface altitude in meters.
        reference_point : np.array
            The reference point as a four element array (lat, lon, alt, mjd).
        """

        # make type consistent (i.e. at least 1d np.array)
        sza_deg = np.atleast_1d(solar_zenith).astype(float).flatten()
        oza_deg = np.atleast_1d(observer_zeniths).astype(float).flatten()
        daz_deg = np.atleast_1d(azimuth_differences).astype(float).flatten()
        o_mjd = np.atleast_1d(mjd).astype(float).flatten()

        if any(sza_deg < 0.) or any(sza_deg > 180.) or any(oza_deg < 0.) or any(oza_deg > 180.):
            raise ValueError('Zenith angles must be within [0, 180] degrees')

        if sza_deg.size == 1:
            oaa = solar_azimuth + daz_deg
            return self.from_zeniths_and_azimuths(solar_zenith, solar_azimuth, mjd, observer_zeniths, oaa,
                                                  observer_alt, reference_point)

        # note: solar_azimuth is ignored for multiple sza

        # get broadcast size
        shape = np.broadcast(o_mjd, sza_deg, oza_deg, daz_deg, observer_alt).shape
        o_mjd = np.broadcast_to(o_mjd, shape)
        sza_deg = np.broadcast_to(sza_deg, shape)
        oza_deg = np.broadcast_to(oza_deg, shape)
        daz_deg = np.broadcast_to(daz_deg, shape)
        observer_alt = np.broadcast_to(observer_alt, shape)

        # update reference point if necessary
        if reference_point is None and self.reference_point is None:
            raise ValueError('The reference point has not been specified.')
        elif reference_point is not None:
            self.reference_point = reference_point

        # get the reference point
        reference_point = sk.Geodetic()
        reference_point.from_lat_lon_alt(*self.reference_point[:3])

        north = -reference_point.local_south
        east = -reference_point.local_west
        up = reference_point.local_up

        def horizontal(azimuth):
            return north * np.cos(np.deg2rad(azimuth)) + east * np.sin(np.deg2rad(azimuth))

        def unit_vector(zenith, azimuth):
            return up * np.cos(np.deg2rad(zenith)) + horizontal(azimuth) * np.sin(np.deg2rad(zenith))

        def make_geodetic(lat, lon):
            lat = (lat + 180.0) % 360.0 - 180.0  # shift to [-180, 180) (shouldn't be necessary)
            if lat > 90.0:  # wrap around the north pole
                lat = 180.0 - lat
                lon += 180.0
            if lat < -90.0:  # wrap around the south pole
                lat = -180.0 - lat
                lon += 180.0
            geo = sk.Geodetic()
            geo.from_lat_lon_alt(lat, lon, 0.0)
            return geo

        # define the sun so that the sza at the reference point is the average sza, and saa = 0.0
        sun_lat = reference_point.latitude + np.average(sza_deg)  # could be > 90.0 (caught by make_geodetic)
        sun_lon = reference_point.longitude
        sun_geo = make_geodetic(sun_lat, sun_lon)

        lines_of_sight = []
        for mjd, sza, oza, daz, alt in zip(o_mjd, sza_deg, oza_deg, daz_deg, observer_alt):
            # set ground intersections south of the sun
            gnd_lat = sun_lat - sza  # could be > 90.0 or < -90.0 (caught by make_geodetic)
            gnd_lon = sun_lon
            gnd_geo = make_geodetic(gnd_lat, gnd_lon)
            up, north, east = gnd_geo.local_up, -gnd_geo.local_south, -gnd_geo.local_west
            # if the ground intersection was wrapped across a pole, then saa = 180.0
            look = -unit_vector(oza, daz + gnd_geo.longitude - gnd_lon)
            obs = gnd_geo.altitude_intercepts(alt, gnd_geo.location, look)[0]
            lines_of_sight.append(sk.LineOfSight(mjd, obs, look))

        self.lines_of_sight = lines_of_sight
        self.sun = sun_geo.local_up

    def transform_los_to_preserve_ground_angles(self):
        """
        Adjusts all lines of sight that intersect the ground so that the solar/observer zenith/azimuth angles will
        be preserved on the osculating sphere. Currently the reference point must be manually set to use this method.
        If the ground intersection of a line of sight coincides with the reference point, no change will occur.

        Notes
        -----
        Locations on the geoid and on the osculating sphere with the same lat/lon coordinates will have the same
        local unit directions, and thus will have the same solar/observer angles if the sun vector and the look
        vector are kept constant.

        Thus the only operation required is to translate the observer by the difference between the geoid/osculating
        sphere locations that share lat/lon coordinates with the ground intersection.
        """

        if self.reference_point is None:
            raise ValueError('The reference point must be set to call transform_los_to_preserve_ground_angles')

        if len(self.lines_of_sight) == 0:
            raise ValueError('Lines of sight must be defined before calling this method')

        if self._osc_adjust:
            raise ValueError('The transform has already been applied; reset the lines_of_sight property to transform '
                             'it again')

        geo = sk.Geodetic()
        geo.from_lat_lon_alt(*self.reference_point[:3])
        osc_radius = geo.osculating_spheroid_radius
        osc_center = geo.osculating_spheroid_center
        lines_of_sight = []
        for los, lon, lat in zip(self._los, self._int_lon, self._int_lat):
            if not any(np.isnan([lat, lon])):  # nadir
                geo.from_lat_lon_alt(lat, lon, 0.0)
                obs = los.observer + (osc_center + osc_radius * geo.local_up) - geo.location
            else:  # limb
                obs = los.observer
            lines_of_sight.append(sk.LineOfSight(los.mjd, obs, los.look_vector))
        self.lines_of_sight = lines_of_sight
        self._osc_adjust = True

    def _calc_intersection(self):
        """
        Calculate the intersection of every line of sight with the earth.
        Fills with nans where no intersection occurs.
        """
        self._int_lon = []
        self._int_lat = []
        geo = sk.Geodetic()
        for i, los in enumerate(self._los):
            i1, i2 = geo.altitude_intercepts(0.0, los.observer, los.look_vector)
            if i1 is None:
                self._int_lon.append(np.nan)
                self._int_lat.append(np.nan)
            else:
                geo.from_xyz(i1)
                self._int_lon.append(geo.longitude)
                self._int_lat.append(geo.latitude)

    def _calc_los_angles(self):
        """
        Calculates the zenith and azimuth angles of the observer from the
        ground intersection of every line of sight with a ground intersection
        """
        self._oza = []
        self._oaa = []

        for los, lon, lat in zip(self._los, self._int_lon, self._int_lat):
            if any(np.isnan([lat, lon])):
                self._oza.append(np.nan)
                self._oaa.append(np.nan)
            else:
                geo = sk.Geodetic()
                geo.from_lat_lon_alt(lat, lon, 0.0)
                self._oza.append(self._zenith(-los.look_vector, geo.local_up))
                self._oaa.append(self._azimuth(-los.look_vector, geo.local_south, geo.local_west))

    def _calc_sun_angles(self):
        """
        Calculates the zenith and azimuth angles of the sun from the ground
        intersection of every line of sight with a ground intersection.
        :return:
        """
        if self._sun is not None:
            self._sza = []
            self._saa = []

            for los, lon, lat in zip(self._los, self._int_lon, self._int_lat):
                if any(np.isnan([lat, lon])):
                    self._sza.append(np.nan)
                    self._saa.append(np.nan)
                else:
                    geo = sk.Geodetic()
                    geo.from_lat_lon_alt(lat, lon, 0.0)
                    self._sza.append(self._zenith(self._sun, geo.local_up))
                    self._saa.append(self._azimuth(self._sun, geo.local_south, geo.local_west))

    def _zenith(self, vec, up):
        dot = np.max([-1.0, np.min([1.0, np.dot(vec, up)])])
        return np.rad2deg(np.arccos(dot))

    def _azimuth(self, vec, south, west):
        return np.rad2deg(np.arctan2(np.dot(vec, -west), np.dot(vec, -south)))

    def _rotation_matrix(self, axis, angle):
        """
        Create a 3D rotation matrix using the Euler-Rodrigues formula

        :param axis: axis of rotation
        :param angle: angle to rotate around the axis (radians)
        :return: 3x3 rotation matrix
        """

        unit_axis = axis / np.linalg.norm(axis)

        half_angle_rad = angle / 2.0
        sin_half_angle_rad = np.sin(half_angle_rad)
        cos_half_angle_rad = np.cos(half_angle_rad)
        a = cos_half_angle_rad
        b = unit_axis[0] * sin_half_angle_rad
        c = unit_axis[1] * sin_half_angle_rad
        d = unit_axis[2] * sin_half_angle_rad

        R = np.zeros((3, 3))
        R[0, 0] = a * a + b * b - c * c - d * d
        R[0, 1] = 2 * (b * c - a * d)
        R[0, 2] = 2 * (b * d + a * c)
        R[1, 0] = 2 * (b * c + a * d)
        R[1, 1] = a * a + c * c - b * b - d * d
        R[1, 2] = 2 * (c * d - a * b)
        R[2, 0] = 2 * (b * d - a * c)
        R[2, 1] = 2 * (c * d + a * b)
        R[2, 2] = a * a + d * d - b * b - c * c

        return R

    @property
    def sun(self):
        return self._sun

    @sun.setter
    def sun(self, value):
        self._sun = value
        self._calc_sun_angles()
        self._has_changed = True

    @property
    def lines_of_sight(self):
        return self._los

    @lines_of_sight.setter
    def lines_of_sight(self, value):
        self._osc_adjust = False  # clean slate allows another osculating sphere adjustment
        self._los = ListWithCallback(value)
        self._calc_intersection()
        self._calc_los_angles()
        self._calc_sun_angles()
        self._has_changed = True

    @property
    def latitude(self):
        return self._int_lat

    @property
    def longitude(self):
        return self._int_lon

    @property
    def solar_zenith(self):
        return self._sza

    @property
    def solar_azimuth(self):
        return self._saa

    @property
    def observer_zenith(self):
        return self._oza

    @property
    def observer_azimuth(self):
        return self._oaa
