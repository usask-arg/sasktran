import numpy as np
import sasktran as sk


class LineOfSight(object):
    """
    Class which represents a single line of sight in SASKTRAN.  A single line of sight is defined by the observer
    position, a timestamp, and a unit look vector.

    Parameters
    ----------
    mjd : float
        Modified julian date for the measurement

    observer : np.array
        Three element array indicating the observer position for the measurement.

    look_vector : np.array
        Three element array which is the unit look vector away from the instrument

    Examples
    --------
    >>> from sasktran import LineOfSight
    >>> los = LineOfSight(mjd=54832.5, observer=[3.6760131547888e+005, 1.0099763136400e+006, -6.871601202127e+006],\
                          look_vector=[2.884568631765662e-001, 7.925287180643269e-001,  5.372996083468238e-001])
    >>> print(los)
    Observer: [367601.31547888, 1009976.31364, -6871601.202127]
    Look: [0.2884568631765662, 0.7925287180643269, 0.5372996083468238]
    MJD: 54832.5
    >>> print(los.mjd)
    54832.5
    >>> print(los.observer)
    [367601.31547888, 1009976.31364, -6871601.202127]
    >>> print(los.look_vector)
    [0.2884568631765662, 0.7925287180643269, 0.5372996083468238]
    """
    def __init__(self, mjd: float, observer: np.array, look_vector: np.array):
        self._mjd = mjd
        self._observer = observer
        self._look_vector = look_vector

    def __repr__(self):
        ret = "Observer: {}\nLook: {}\nMJD: {}".format(self._observer, self._look_vector, self._mjd)

        return ret

    @property
    def mjd(self):
        return self._mjd

    @property
    def observer(self):
        return self._observer

    @property
    def look_vector(self):
        return self._look_vector

    def ground_intersection(self, altitude: float=0.0):
        """
        Returns an sk.Geodetic object containing the location where the line of sight intersects the earth at the given
        altitude. Returns None if the intersection does not exist.

        Parameters
        ----------
        altitude : float
            The altitude of the desired intersection in meters. Default 0.

        Examples
        --------
        >>> from sasktran import LineOfSight
        >>> los = LineOfSight(mjd=54832.5, observer=[3.6760131547888e+005, 1.0099763136400e+006, -6.871601202127e+006],\
                              look_vector=[2.878657667526608e-001, 7.909046939869273e-001, 5.400028382900848e-001])
        >>> print(los.ground_intersection(altitude=1000.0))
        ISKGeodetic: IAU 1976
         Latitude: -57.4997267221534, Longitude: 69.99999999999979, Altitude: 999.9999956705142
        """
        geo = sk.Geodetic()
        int1, in2 = geo.altitude_intercepts(altitude, self.observer, self.look_vector)
        np.set_printoptions()
        if int1 is None:
            return
        else:
            geo.from_xyz(int1)
            return geo

    def tangent_location(self):
        """
        Returns an sk.Geodetic object containing the location where the line of sight is tangent to the surface of the
        earth. Returns None if the line of sight intersects the earth.

        Examples
        --------
        >>> from sasktran import LineOfSight
        >>> los = LineOfSight(mjd=54832.5, observer=[3.6760131547888e+005, 1.0099763136400e+006, -6.871601202127e+006],\
                              look_vector=[2.884568631765662e-001, 7.925287180643269e-001,  5.372996083468238e-001])
        >>> print(los.tangent_location())
        ISKGeodetic: IAU 1976
         Latitude: -57.49972673428996, Longitude: 69.99999999999979, Altitude: 10000.000072206138
        """
        geo = sk.Geodetic()
        geo.from_tangent_point(self.observer, self.look_vector)
        if geo.altitude < 0.0:
            return
        else:
            return geo