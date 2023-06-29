import unittest
import sasktran as sk
import numpy as np


class TestHITRAN(unittest.TestCase):
    def setUp(self):
        pass

    @unittest.skip
    def test_simple_access(self):
        opt_prop = sk.HITRANChemical('H2O')

        opt_prop.calculate_cross_sections(sk.MSIS90(), 0, 0, 10000, 54372, 1e7 / np.array([1353, 1354]))

    @unittest.skip
    def test_cached_access(self):
        opt_prop = sk.HITRANChemical('O2')
        opt_prop_nocache = sk.HITRANChemical('O2', use_cache=False)

        atmo = sk.Atmosphere()
        atmo['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
        atmo['O2'] = sk.Species(opt_prop, sk.MSIS90(include_o2=True), 'SKCLIMATOLOGY_O2_CM3')

        atmo_nocache = sk.Atmosphere()
        atmo_nocache['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
        atmo_nocache['O2'] = sk.Species(opt_prop_nocache, sk.MSIS90(include_o2=True), 'SKCLIMATOLOGY_O2_CM3')

        tanalts_km = np.arange(10, 11, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        engine = sk.EngineHR(atmosphere=atmo, geometry=geometry)
        engine_nocache = sk.EngineHR(atmosphere=atmo, geometry=geometry)

        engine.wavelengths = [762.914, 762.915]
        engine_nocache.wavelengths = [762.914, 762.915]

        rad_cache = engine.calculate_radiance()
        rad_nocache = engine_nocache.calculate_radiance()

        np.testing.assert_array_equal(rad_cache, rad_nocache)

        engine.wavelengths = [763.425, 763.426]
        engine_nocache.wavelengths = [763.425, 763.426]

        rad_cache = engine.calculate_radiance()
        rad_nocache = engine_nocache.calculate_radiance()

        np.testing.assert_array_equal(rad_cache, rad_nocache)
