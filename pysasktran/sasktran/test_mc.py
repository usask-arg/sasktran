import unittest
import sasktran as sk
import numpy as np
import sys


class TestEngineMC(unittest.TestCase):

    def _default_config(self):
        atmosphere = sk.Atmosphere()
        atmosphere['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere.brdf = 1.0
        mjd = 54832.0
        lat = -57.5
        lng = 70.0
        sza = 60.0
        saa = 157.5
        rayazi = 0.0
        geometry = sk.VerticalImage()
        tanalts_km = np.arange(10.0, 18.0, 2.0)
        geometry.from_sza_saa(sza, saa, lat, lng, tanalts_km, mjd, rayazi, satalt_km=600, refalt_km=0)
        # in the equivalent c++ test, SetRaysFromTangentHeightArray has const look, variable obs,
        # but here from_sza_saa has const obs, variable look,
        # therefore override to match c++ test exactly
        geometry.lines_of_sight = [
            sk.LineOfSight(
                mjd=mjd,
                observer=np.array([3.6608172984716902e+05, 1.0058012864299840e+06, -6.8744281347604860e+06]),
                look_vector=np.array([2.8845686317656621e-01, 7.9252871806432690e-01, 5.3729960834682389e-01])
            ),
            sk.LineOfSight(
                mjd=mjd,
                observer=np.array([3.6776456985338684e+05, 1.0104248513476220e+06, -6.8736649392269310e+06]),
                look_vector=np.array([2.8845686317656621e-01, 7.9252871806432690e-01, 5.3729960834682389e-01])
            ),
            sk.LineOfSight(
                mjd=mjd,
                observer=np.array([3.6944996016853803e+05, 1.0150554231814672e+06, -6.8728969933124129e+06]),
                look_vector=np.array([2.8845686317656621e-01, 7.9252871806432690e-01, 5.3729960834682389e-01])),
            sk.LineOfSight(
                mjd=mjd,
                observer=np.array([3.7113791328355065e+05, 1.0196930362500623e+06, -6.8721242737504719e+06]),
                look_vector=np.array([2.8845686317656621e-01, 7.9252871806432690e-01, 5.3729960834682389e-01])
            )
        ]
        wavelengths = [600.0, 340.0]
        engine = sk.EngineMC(geometry=geometry, atmosphere=atmosphere, wavelengths=wavelengths)
        engine.max_photons_per_los = 500                # 500 photons (radiance samples per scatter order)
        engine.target_std = 0.0                         # trace all 500 photon paths; don't stop at a target precision
        engine.min_relative_path_weight = 0.0           # don't truncate even if contributions are small
        engine.num_ray_tracing_shells = 101             # 100 atmospheric layers
        engine.ray_tracer_type = 0                      # straight rays intersecting spherical shells only
        engine.solar_table_type = 0                     # calculate single scatter source terms on the fly; no cache
        engine.optical_property_integration_type = 0    # no adaptive integration to split large optical depths
        engine.optical_table_type = 0                   # 1-dimensional atmosphere
        engine.debug_mode = 1234                        # turn off multithreading, fix the rng seed
        engine.scatter_position_resolution = 50.0       # solve scatter positions to within 50 m
        engine.min_fraction_higher_order = 1.0          # no optimization: always trace photons to the maximum order
        engine.num_orders_of_scatter = 50               # trace 50 orders of scatter
        return engine

    def test_mc(self):
        engine = self._default_config()

        engine_output = engine.calculate_radiance()
        hardcodeMC = np.array([
            [0.055773125345084, 0.042552458325431, 0.037556817927853, 0.029011461653996],
            [0.171188856097297, 0.166793306703896, 0.168159914831352, 0.154864050048638]
        ])
        diff = np.abs((engine_output.radiance - hardcodeMC) / hardcodeMC)
        maxdiff = np.max(diff)
        self.assertAlmostEqual(maxdiff, 0.0, places=5)

    def test_mc_simultaneous(self):
        engine = self._default_config()
        engine.simultaneous_wavelength = True

        engine_output = engine.calculate_radiance()
        hardcodeMC = np.array([
            [0.055773125345084, 0.042552458325431, 0.037556817927853, 0.029011461653996],
            [0.174631688252108, 0.207222859479211, 0.159144998846106, 0.139795410371997]
        ])
        diff = np.abs((engine_output.radiance - hardcodeMC) / hardcodeMC)
        maxdiff = np.max(diff)
        self.assertAlmostEqual(maxdiff, 0.0, places=5)

    # @unittest.skipIf(sys.platform == 'linux', "Temporarily skipping raman on linux")
    def test_mc_inelastic(self):
        engine = self._default_config()
        engine.simultaneous_wavelength = True
        engine.atmosphere['rayleigh'] = sk.Species(sk.InelasticRayleigh(), sk.MSIS90())
        engine.optical_property_wavelengths = [330., 340., 350., 590., 600., 610.]
        engine.max_raman_orders = [3, 2]  # calculate paths with 1 Raman event up to 5 total scatters, and 2 up to 2
        engine.min_fraction_higher_raman_order = 0.15

        engine_output = engine.calculate_radiance()
        hardcodeMC = np.array([
            [0.056824427667601, 0.051363737943832, 0.031991126813500, 0.029536049150431],
            [0.161535158610863, 0.164237525384473, 0.138946123097123, 0.135347155516999]
        ])
        diff = np.abs((engine_output.radiance - hardcodeMC) / hardcodeMC)
        maxdiff = np.max(diff)
        self.assertAlmostEqual(maxdiff, 0.0, places=5)

    # @unittest.skipIf(sys.platform == 'linux', "temporarily skipping raman on linux")
    def test_mc_raman(self):
        engine = self._default_config()
        engine.simultaneous_wavelength = True
        engine.atmosphere['rayleigh'] = sk.Species(sk.InelasticRayleigh(), sk.MSIS90())
        engine.optical_property_wavelengths = [330., 340., 350., 590., 600., 610.]
        engine.secondary_output = 2
        engine.max_raman_orders = [3, 2]
        engine.min_fraction_higher_raman_order = [0.05] * 3

        engine_output = engine.calculate_radiance()
        hardcodeMC = np.array([
            [6.40423107195561e-05, 9.00503515726087e-05, 7.17288311565054e-05, 6.05103163474209e-05],
            [-3.56401502763240e-05, -9.05373136546410e-05, -3.17115580209637e-05, -5.91555494448382e-05]
        ])
        diff = np.abs((engine_output.ring_spectrum - hardcodeMC) / hardcodeMC)
        maxdiff = np.max(diff)
        self.assertAlmostEqual(maxdiff, 0.0, places=5)

    def test_mc_amf(self):
        mjd = 58330.75
        observer = sk.Geodetic()
        observer.from_lat_lon_alt(0., 260., 35.786e6)
        lats, lons = 50., 260.
        geometry = sk.NadirGeometry()
        geometry.from_lat_lon(mjd, observer, lats, lons)

        engine = self._default_config()
        engine.geometry = geometry
        engine.wavelengths = [440.]
        engine.air_mass_factor = 1
        engine.air_mass_factor_shells = np.linspace(0, 5e4, 11)
        engine_output = engine.calculate_radiance()
        hardcodeMC = np.array([
            3.6577441, 3.57388592, 3.34411161, 3.16308768, 3.08838479,
            3.05418787, 3.03382158, 3.02590744, 3.0218465, 3.01734845
        ])
        diff = np.abs((engine_output.air_mass_factor - hardcodeMC) / hardcodeMC)
        maxdiff = np.max(diff)
        self.assertAlmostEqual(maxdiff, 0.0, places=5)