from typing import List, Dict, Any, Tuple
import numpy as np
from .geometry import Geometry
from .atmosphere import Atmosphere
from .engine import Engine


# ------------------------------------------------------------------------------
#           class EngineSO(Engine):
# ------------------------------------------------------------------------------
class EngineSO(Engine):
    """
    The ``SO`` engine solves the scalar radiative transfer equation using the successive orders technique. It was the
    first engine developed for the sasktran framework. The **SO** engine is deprecated and is no longer under development.
    It has been superceded by the :class:`~.EngineHR` engine which also uses a successive orders algorithm but has significant
    accuracy improvements as well as additional important features such as weighting function generation and adaptive ray
    splitting.

    The theoretical basis for the **SO** technique can be found in [soBourassa2008]_. The ``SO`` engine, which was developed for
    ozone retrievals from odin-osiris measurements, is optimized for computational speed at the expense of memory.
    The model solves the radiative transfer over a geographic region surrounding a location known as the *Reference Point*;
    this is typically the location where the solution is required (eg location of limb scan etc.) The model
    approximates the Earth as a sphere with the radius matching the North-South curvature of the true Earth
    at the *Reference Point*.

    The model works in two stages. The first stage is a *geometry* stage where all of the rays used by the model are traced
    through the atmosphere and split into manageable segments throughout the atmosphere. Many geometrical factors and calculations generated
    during the ray-tracing are internally cached and tabulated for later use by the second stage. The second stage peforms all the
    *wavelength specific* calculations.

    The primary advantage of the ``SO`` engine is that the ray tracing is performed
    once and is shared between subsequent wavelength calculations. This saving can be significant within atmospheric retrievals where
    the optical properties of the atmosphere are modified multiple times as part of the fitting process.

    On the other hand this strength is also a weakness as (i) the ray tracing is not adaptive and cannot split rays when the
    optical depth is too large and (ii) the internal memory footprint can get large if many diffuse profiles are required across
    the region of interest. The model works well for scenarios that avoids these two problem conditions, for example
    uv-visible calculations with the Sun well above the horizon and small horizontal gradients in atmospheric properties. The
    :class:`~.EngineHR` addresses these issues and is the recommended engine for successive orders calculations.

    References
    ----------
    .. [soBourassa2008] | A.E. Bourassa, D.A. Degenstein, E.J. Llewellyn, SASKTRAN: A spherical geometry radiative transfer code for efficient estimation of limb scattered sunlight, Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 109, Issue 1, January 2008, Pages 52-73, ISSN 0022-4073, http://dx.doi.org/10.1016/j.jqsrt.2007.07.007.

    """

    # ------------------------------------------------------------------------------
    #           __init__
    # ------------------------------------------------------------------------------

    def __init__(self,
                 geometry: Geometry = None,
                 atmosphere: Atmosphere = None,
                 wavelengths: List[float] = None,
                 options: Dict[str, Any] = None):
        super().__init__('SO', geometry, atmosphere, wavelengths, options)

    # ------------------------------------------------------------------------------
    #           reference_point_target_altitude
    # ------------------------------------------------------------------------------
    @property
    def userdefined_reference_point_target_altitude(self) -> float:
        """
        The target altitude to be used when determining the reference point. The target altitude is used in zenith
        and nadir observations to find the location where a ray intersects this altitude. This altitude is then used to
        find an average reference point location. The target altitude is used with the target range variable in limb
        viewing geometries to apply an increased weighting to lines of sight tangential in the vicinity of the target
        altitude. This encourages the reference point to be closer to the lines of sight which are tangential in the
        region of the reference point.

        Parameters
        ----------
        heightm : float
            The target height in meters.
        """
        return self.options.get('setreferencepoint_targetaltitude', np.nan)

    @userdefined_reference_point_target_altitude.setter
    def userdefined_reference_point_target_altitude(self, heightm: float):
        self.options['setreferencepoint_targetaltitude'] = heightm
        
    # ------------------------------------------------------------------------------
    #               userdefined_reference_point_target_range:
    # ------------------------------------------------------------------------------
    @property
    def userdefined_reference_point_target_range(self) -> float:
        """
        The reference point target range parameter. This variable is only used when calculating the reference point
        from limb viewing geometries. Specifies the altitude range above and below the target altitude for enhanced
        weighting of limb viewing lines of sight, default 15000 meters.

        Parameters
        ----------
        range : float
            The height range value in meters.
        """
        return self.options.get('setreferencepoint_targetrange', np.nan)
        
    @userdefined_reference_point_target_range.setter
    def userdefined_reference_point_target_range(self, range: float):
        self.options['setreferencepoint_targetrange'] = range

    # ------------------------------------------------------------------------------
    #           userdefined_ground_altitude
    # ------------------------------------------------------------------------------
    @property
    def userdefined_ground_altitude(self) -> float:
        """
        Sets the altitude in meters of the ground shell above the oblate spheroid/sea level.

        Parameters
        ----------
        heightm : float
            The height of the ground surface above sea-level in meters
        """
        return self.options.get('setgroundaltitude', np.nan)
    
    @userdefined_ground_altitude.setter
    def userdefined_ground_altitude(self, heightm: float):
        self.options['setgroundaltitude'] = heightm

    # ------------------------------------------------------------------------------
    #           userdefined_upper_bound_altitude
    # ------------------------------------------------------------------------------
    @property
    def userdefined_upper_bound_altitude(self) -> float:
        """
        The maximum altitude of the atmosphere in meters used when considering the reference point

        Parameters
        ----------
        heightm : float
            The maximum altitude in meters
        """
        return self.options.get('setupperboundaltitude', None)

    @userdefined_upper_bound_altitude.setter
    def userdefined_upper_bound_altitude(self, heightm: float):
        self.options['setupperboundaltitude'] = heightm

    # ------------------------------------------------------------------------------
    #           userdefined_lower_bound_altitude
    # ------------------------------------------------------------------------------
    @property
    def userdefined_lower_bound_altitude(self) -> float:
        """
        The lower altitude of the atmosphere in meters used when considering the reference point

        Parameters
        ----------
        heightm : float
            The minimum altitude in meters
        """
        return self.options.get('setlowerboundaltitude', None)

    @userdefined_lower_bound_altitude.setter
    def userdefined_lower_bound_altitude(self, heightm: float):
        self.options['setlowerboundaltitude'] = heightm

    # ------------------------------------------------------------------------------
    #           userdefined_reference_point
    # ------------------------------------------------------------------------------
    @property
    def userdefined_reference_point(self) -> np.ndarray:
        """
        Sets the reference point to the given location.

        Parameters
        ----------
        location : array[4]
             three element array specifing [latitude, longitude, height, mjd]. Note that
             the height field is ignored.

        """
        return self.options.get('setreferencepoint', None)

    @userdefined_reference_point.setter
    def userdefined_reference_point(self, location: Tuple[float, float, float, float]):
        self.options['setreferencepoint'] = np.array(location)

    # ------------------------------------------------------------------------------
    #           userdefined_sun
    # ------------------------------------------------------------------------------
    @property
    def userdefined_sun(self) -> np.ndarray:
        """
        Manually sets the sun position. Use with caution as many line-of-sight helper functions will also
        set the sun's position inside the internal engine.

        Parameters
        ----------
        sun : array[3]
            The three element array of the sun's unit vector. Expressed in global geographic coordinates and is in the
            directions from the Earth toward the sun.
        """
        return self.options.get('setsun', None)

    @userdefined_sun.setter
    def userdefined_sun(self, sun: Tuple[float, float, float]):
        self.options['setsun'] = np.array(sun)

    # ------------------------------------------------------------------------------
    #           num_orders_of_scatter
    # ------------------------------------------------------------------------------
    @property
    def num_orders_of_scatter(self) -> int:
        """
        Sets the number of orders of scatter used in the successive orders approximation. This number should be set high
        enough to ensure the solution has converged. We have found that values as high as 50 are often required at UV-VIS
        wavelengths close to 350 nm to ensure convergence

        Parameters
        ----------
        n : int
             The number of successive orders to be calculated.
        """
        return int(self.options.get('numordersofscatter', 10))

    @num_orders_of_scatter.setter
    def num_orders_of_scatter(self, n: int):
        self.options['numordersofscatter'] = n

    # ------------------------------------------------------------------------------
    #           num_diffuse_profiles
    # ------------------------------------------------------------------------------

    @property
    def num_diffuse_profiles(self) -> int:
        """
        Sets the number of diffuse profiles used in the calculation. The default is 1, implying there is no spatial variation
        variation in the diffuse signal across the calculation region. We recommend setting odd values for this number so the
        middle profile is centered on the reference point. Typical values are the odd numbers up to ~11 although smaller
        values (3 or 5) probably capture the variation of the diffuse signal across the scene for most situations.  This
        parameter is not intended to accurately capture the variation in the diffuse profile across the terminator

        Parameters
        ----------
        n : int
             The number of diffuse profiles.
        """
        return int(self.options.get('numdiffuseprofiles', 1))

    @num_diffuse_profiles.setter
    def num_diffuse_profiles(self, n: int):
        self.options['numdiffuseprofiles'] = n

    # ------------------------------------------------------------------------------
    #           userdefined_shells
    # ------------------------------------------------------------------------------
    @property
    def userdefined_shells(self) -> np.ndarray:
        """
        Configures the model to use the altitude shells defined by the user. The
        shell altitudes must be in ascending order and are specified in meters above sea level.
        By default, shell heights are evenly spaced at 1000 m intervals from 0 to 100,000 m.
        The diffuse points, by default, are placed in the middle of the shells and the optical
        properties are placed both in the middle and on the boundary of each shell.

        Parameters
        ----------
        heights : array[n]
            The array of shell altitudes in meters above sea level.
        """
        return self.options.get('configureuserdefinedshells', None)

    @userdefined_shells.setter
    def userdefined_shells(self, heights: List[float]):
        self.options['configureuserdefinedshells'] = np.array(heights)

    # ------------------------------------------------------------------------------
    #           userdefined_diffuse_heights
    # ------------------------------------------------------------------------------
    @property
    def userdefined_diffuse_heights(self) -> np.ndarray:
        """
        Manually overides the default placement of diffuse points. This is an advanced option and most users do not need to
        call this property as the default values are usually sufficient.

        Parameters
        ----------
        heights : array[n]
            The array of diffuse altitudes specified in meters above sea level. Diffuse points are placed at each altitude
            on each diffuse profile.
        """
        return self.options.get('manualdiffuseheights', None)

    @userdefined_diffuse_heights.setter
    def userdefined_diffuse_heights(self, heights: List[float]):
        self.options['manualdiffuseheights'] = np.array(heights)

    # ------------------------------------------------------------------------------
    #           userdefined_diffuse_incoming_resolution
    # ------------------------------------------------------------------------------
    @property
    def userdefined_diffuse_incoming_resolution(self) -> np.ndarray:
        """
        Specifies four parameters that specify the distribution of azimuth and zenith angles across the incoming unit sphere of
        each diffuse point. This distribution of points is used to calculate the incoming signal which is then scattered to the
        outbound unit sphere of each diffuse point.The algorithm breaks the zenith angles of the incoming unit sphere of each
        diffuse point into 3 regions,

        * Upward incoming signal from the lower (ground) regions.
        * Horizontal incoming signalfrom the angles close to the limb
        * Downward signal from altitudes above the diffuse point.

        This is an alternative method to using property **userdefined_diffuse_incoming_zenith_angles** and **userdefined_diffuse_incoming_azimuth_angles**
        set the incoming zenith and azimuth angles. Note that the default values are set in these other functions.

        Parameters
        -----------
        parameters : array[4] parameters:
            An array of 4 integers that configure the distribution of incoming points on the incoming unit sphere of each diffuse
            point,

            * [0] Ground resolution: number of zenith angles between 180 and 100 degrees on the incoming diffuse unit sphere.
            * [1] Horizon resolution: number of zenith angles between 100 and 80 degrees on the incoming diffuse unit sphere.
            * [2] Atmosphere resolution: number of zenith angles between 0 and 80 degrees on the incoming diffuse unit sphere.
            * [3] Number of azimuths: the number of horizontal azimuths around the incoming diffuse unit sphere.
        """
        return self.options.get('diffuseincomingresolution', None)

    @userdefined_diffuse_incoming_resolution.setter
    def userdefined_diffuse_incoming_resolution(self, parameters: Tuple[float, float, float, float]):
        self.options['diffuseincomingresolution'] = np.array(parameters)

    # ------------------------------------------------------------------------------
    #           userdefined_diffuse_incoming_zenith_angles
    # ------------------------------------------------------------------------------
    @property
    def userdefined_diffuse_incoming_zenith_angles(self) -> np.ndarray:
        """
        Configures the incoming zenith angles for the incoming unit sphere of each diffuse point. This is an array of ascending
        zenith angles in degrees. Each element of the array specifies the zenith angle at which the solid angle of the
        corresponding segment finishes. I.E. The first value should not be 0.0 as that implies the solid angle segment starts
        and ends at 0.0. Likewise the last value should be 180.0 as the last segment will finish at 180.0 zenith angle.

        Parameters
        ----------
        zenithangles : array[n]
            The array of zenith angles in degrees in ascending order. The default value is, [15, 30, 40, 50, 60, 70, 75, 80,
            85, 87, 89, 90, 91, 93, 95, 100, 105, 110, 120, 130, 140, 150, 165, 180].
        """
        return self.options.get('configureincomingzenithangles', None)

    @userdefined_diffuse_incoming_zenith_angles.setter
    def userdefined_diffuse_incoming_zenith_angles(self, zenithangles: List[float]):
        self.options['configureincomingzenithangles'] = np.array(zenithangles)

    # ------------------------------------------------------------------------------
    #           userdefined_diffuse_incoming_azimuth_angles
    # ------------------------------------------------------------------------------
    @property
    def userdefined_diffuse_incoming_azimuth_angles(self) -> np.ndarray:
        """
        Configures the incoming azimuth angles for the incoming unit sphere of each diffuse point.

        Parameters
        ----------
        azimuthangles: array[n]
            This is an array of ascending azimuth angles in degrees. The default value is [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
        """
        return self.options.get('configureincomingazimuthangles', None)

    @userdefined_diffuse_incoming_azimuth_angles.setter
    def userdefined_diffuse_incoming_azimuth_angles(self, azimuthangles: List[float]):
        self.options['configureincomingazimuthangles'] = np.array(azimuthangles)