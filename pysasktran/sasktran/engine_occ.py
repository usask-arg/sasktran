from typing import List, Dict, Any, Tuple
import numpy as np
from .geometry import Geometry
from .atmosphere import Atmosphere
from .engine import Engine


class EngineOCC(Engine):
    """"
    The ``OCC`` engine was built as a radiative transfer model for solar occultation observations from a satellite. It was
    originally built as a plug-in subsititute engine for the Fortran code used in the `ACE-FTS <http://www.ace.uwaterloo.ca/instruments_acefts.php>`_
    analysis and is well-suited for atmospheric transmission calculations in near infra-red and infra-red micro-windows. The algorithm
    traces curved rays using a spherically symmetric version of Snells law,  Thompson et al. 1982. The refractive index of the
    atmosphere follows the dry-air formula published by Filippenko 1982 which only requires a single profile of pressure
    and temperature; we note that the atmospheric refractive index does not include effects due to moist air which
    may be significant at tropospheric altitudes.

    The code for the OCC engine was frequently compared to the original Fortran code during development to ensure
    they are identical. The engine calculates transmission through the atmosphere, rather than radiance, and has
    the following features:

        * curved ray tracing in a spherically symmetric atmosphere. The path length of the
          curved ray is the same as the ACE-FTS code to within 1 micron.
        * The code provides efficient calculation of Voigt profiles across micro-windows in the near infra-red and infra-red.

    References
    ----------
        **Dennis A. Thompson**, Theodore J. Pepin, and Frances W. Simon. Ray tracing in a refracting spherically symmetric atmosphere.
        J. Opt. Soc. Am., 72(11):1498–1501, Nov 1982. URL: `http://www.osapublishing.org/abstract.cfm?URI=josa-72-11-1498 <http://www.osapublishing.org/abstract.cfm?URI=josa-72-11-1498>`_,
        doi:10.1364/JOSA.72.001498.

        **A. V. Filippenko**. The importance of atmospheric differential refraction in spectrophotometry. Publications of the Astronomical Society of the Pacific , 94:715–721, August 1982. doi:10.1086/131052.

    Example
    --------
    Here is an example that calculates the optical depth at 5000 wavenumbers across a CO2 band at 934 to 936 cm-1. This example is almost the same as ACE-FTS scan ace.sr20314 except the sun is no longer
    in the correct location. The test example creates a user defined profile of CO2 number density used for CO2 absorption and the MSIS climatology for Rayleigh extinction. Lines of sights are specified
    as an array of tangent altitudes and the satellite is placed at 600 km. Optical depth through the atmosphere is calculated for each tangent altitude across the 4000 wavenumbers::

        import numpy as np
        import matplotlib.pyplot as plt
        import sasktran as sk

        def test_occ_engine():

            co2profile =  np.array( [    0.000, 9.5620350469e+15,  1000.000, 8.5604676285e+15,  2000.000, 7.7062120091e+15,  3000.000, 6.9531991470e+15,  4000.000, 6.2702731320e+15,  5000.000, 5.6375862919e+15,
                            6000.000, 5.0651291274e+15,  7000.000, 4.4975838604e+15,  8000.000, 3.9468136861e+15,  9000.000, 3.4348048814e+15, 10000.000, 2.9871067830e+15, 11000.000, 2.5656416175e+15,
                           12000.000, 2.1874053365e+15, 13000.000, 1.8533816021e+15, 14000.000, 1.6023327829e+15, 15000.000, 1.3568375796e+15, 16000.000, 1.1279532788e+15, 17000.000, 9.7672446573e+14,
                           18000.000, 8.4173283897e+14, 19000.000, 7.1576699275e+14, 20000.000, 6.2070908062e+14, 21000.000, 5.2364297410e+14, 22000.000, 4.3181248841e+14, 23000.000, 3.7860567983e+14,
                           24000.000, 3.2428122678e+14, 25000.000, 2.7110791383e+14, 26000.000, 2.3526785090e+14, 27000.000, 2.0344493146e+14, 28000.000, 1.7304110039e+14, 29000.000, 1.4714113133e+14,
                           30000.000, 1.2544180466e+14, 31000.000, 1.0738346125e+14, 32000.000, 9.2442937053e+13, 33000.000, 8.0342242281e+13, 34000.000, 6.9455591820e+13, 35000.000, 5.9095214441e+13,
                           36000.000, 5.0374561563e+13, 37000.000, 4.3515754800e+13, 38000.000, 3.7794009046e+13, 39000.000, 3.2874895083e+13, 40000.000, 2.8685628465e+13, 41000.000, 2.4978923024e+13,
                           42000.000, 2.1682117851e+13, 43000.000, 1.8864809592e+13, 44000.000, 1.6431826141e+13, 45000.000, 1.4348899126e+13, 46000.000, 1.2595260698e+13, 47000.000, 1.1093125765e+13,
                           48000.000, 9.8376261311e+12, 49000.000, 8.8026864921e+12, 50000.000, 7.8993464447e+12, 51000.000, 7.0038829664e+12, 52000.000, 6.0771348455e+12, 53000.000, 5.2887296427e+12,
                           54000.000, 4.6787494256e+12, 55000.000, 4.1667051367e+12, 56000.000, 3.6751620506e+12, 57000.000, 3.1811011797e+12, 58000.000, 2.7604364326e+12, 59000.000, 2.4249492298e+12,
                           60000.000, 2.1420175118e+12, 61000.000, 1.8772791073e+12, 62000.000, 1.6195294613e+12, 63000.000, 1.3994285676e+12, 64000.000, 1.2229247260e+12, 65000.000, 1.0734951007e+12,
                           66000.000, 9.3270881894e+11, 67000.000, 7.9345730980e+11, 68000.000, 6.7795327304e+11, 69000.000, 5.9174431127e+11, 70000.000, 5.2173619614e+11, 71000.000, 4.5523334147e+11,
                           72000.000, 3.8840635314e+11, 73000.000, 3.3304529951e+11, 74000.000, 2.9045416707e+11, 75000.000, 2.5517516779e+11, 76000.000, 2.2127024526e+11, 77000.000, 1.8582366434e+11,
                           78000.000, 1.5596546276e+11, 79000.000, 1.3362547386e+11, 80000.000, 1.1541990113e+11, 81000.000, 9.8756976417e+10, 82000.000, 8.2629944315e+10, 83000.000, 6.8563739750e+10,
                           84000.000, 5.6814363571e+10, 85000.000, 4.6797966799e+10, 86000.000, 3.8795906044e+10, 87000.000, 3.2908654369e+10, 88000.000, 2.7811184596e+10, 89000.000, 2.2974282383e+10,
                           90000.000, 1.8716304570e+10, 91000.000, 1.5254396937e+10, 92000.000, 1.2548308770e+10, 93000.000, 1.0295593615e+10, 94000.000, 8.3338827301e+09, 95000.000, 6.6488536883e+09,
                           96000.000, 5.2936443303e+09, 97000.000, 4.2242029799e+09, 98000.000, 3.3594428424e+09, 99000.000, 2.6511281727e+09]).reshape( [100,2])
            altitudes = co2profile[:,0]
            values    = co2profile[:,1]
            co2numberdensity = sk.ClimatologyUserDefined(altitudes, {'SKCLIMATOLOGY_CO2_CM3': values})

            atmosphere = sk.Atmosphere()
            atmosphere['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
            atmosphere['co2']      = sk.Species(sk.HITRANChemical('CO2'), co2numberdensity )

            tanalts = np.array( [ 95.542934418655, 93.030998230911, 90.518486023880, 87.999366761185, 85.485855103470, 82.971916199661, 80.457603455521, 77.942962647415, 75.421955109573, 72.906806946732,
                         70.391479493118, 67.869934082962, 65.354396820999, 62.838825226761, 60.317182541824, 57.801700592972, 55.286336899734, 52.765050888992, 50.250070572830, 47.735359192825,
                         45.220966340042, 42.682148825007, 40.254586233282, 37.901745439957, 35.689252976524, 33.619203107470, 31.633878541417, 29.706157206720, 27.941217916525, 26.315136637345,
                         24.759740931714, 23.273057386890, 21.701357220703, 20.435203333687, 19.296175927224, 18.238125008002, 17.137857798933, 15.665431416870, 14.623809766528, 13.581115284387,
                         12.793781944543, 11.901170623281, 10.978181776555, 10.1851695349872, 9.4383271471788, 8.7424541473265, 8.0540969039894, 7.5483134223615, 7.0824804787830, 6.7903857771487,
                         6.3015475934096] )

            mjd = 54242.26386852                    # MJD("2007-05-22 06:19:58.24") Use ACE-FTS scan ace.sr20314 as an example.
            lat = 68.91
            lng = -79.65
            sza = 60.0
            saa = 157.5
            rayazi = 0.0
            geometry = sk.VerticalImage()
            geometry.from_sza_saa(sza, saa, lat, lng, tanalts, mjd, rayazi, satalt_km=600, refalt_km=20)

            wavenum = np.arange( 934.0, 936.0, 0.0005)
            wavelengths = 1.0E7/wavenum
            engine = sk.EngineOCC(geometry=geometry, atmosphere=atmosphere, wavelengths=wavelengths)

            extinction = engine.calculate_radiance()
            plt.plot( wavenum, np.log10(extinction[:,25]))
            plt.xlabel('Wavenumber cm-1')
            plt.ylabel('Log10(Optical Depth)')
            plt.title('CO2 Optical Depth. Tan Alt={:7.3f} km'.format(tanalts[25]))
            plt.show()
    """
    # ------------------------------------------------------------------------------
    #           __init__
    # ------------------------------------------------------------------------------

    def __init__(self,
                 geometry: Geometry = None,
                 atmosphere: Atmosphere = None,
                 wavelengths: List[float] = None,
                 options: Dict[str, Any] = None):
        super().__init__('OCC', geometry, atmosphere, wavelengths, options)

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
    #           SetReferencePoint_TargetAltitude
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
    #           ray_tracing_wavenumber
    # ------------------------------------------------------------------------------
    @property
    def ray_tracing_wavenumber(self) -> float:
        """
        Sets the wavenumber in cm-1 used for tracing curved rays through the atmosphere. Rays are only traced once in
        each model and the trajectories are shared between all wavenumber extinction calculations.

        Parameters
        ----------
        wavenumber : float
            The wavenumber used for tracing rays.

        Returns
        -------

        """
        return self.options.get('setraytracingwavenumber', np.nan)

    @ray_tracing_wavenumber.setter
    def ray_tracing_wavenumber(self, wavenumber: float):
        self.options['setraytracingwavenumber'] = wavenumber

    # ------------------------------------------------------------------------------
    #           AddLineOfSightFromTangentAlt
    # ------------------------------------------------------------------------------
    @property
    def userdefined_lines_of_sight_from_tangent_alt(self) -> np.ndarray:
        """
        Adds a line of sight generated from the tangent height and an observer height.

        Parameters
        ----------
        tangentheights : Array[2*n]
            A 2*n element array specifying ``n`` lines of sight, the first element of each line of sight entry is the
            target tangent height and the 2nd element is the observer height. Both are in meters.
        """
        return self.options.get('addlinesofsightfromtangentalt', None)

    @userdefined_lines_of_sight_from_tangent_alt.setter
    def userdefined_lines_of_sight_from_tangent_alt(self, tangentheights: List[float]):
        self.options['addlinesofsightfromtangentalt'] = np.array(tangentheights)

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
        return self.options.get('setraytracingshells', None)

    @userdefined_shells.setter
    def userdefined_shells(self, heights: List[float]):
        self.options['setraytracingshells'] = np.array(heights)
