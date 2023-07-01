import unittest
import sasktran as sk
import numpy as np


def _netcdf4_available():
    try:
        import netCDF4
        return True
    except ImportError:
        return False


class TestSpecies(unittest.TestCase):
    def setUp(self):
        self.test_altitudes = np.linspace(500, 99500, 100)

    def _basic_access_tests(self, species, clim: sk.Climatology):
        clim.get_parameter(species, 0, 0, self.test_altitudes, 54372)

    def test_aerosol(self):

        values = np.ones_like(self.test_altitudes)
        climatology = sk.SpeciesAerosol(self.test_altitudes, {'SKCLIMATOLOGY_AEROSOL_CM3': values},
                                        {'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': values * 0.08,
                                         'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': values * 1.6}, 'H2SO4')

        climatology.optical_property.calculate_cross_sections(sk.MSIS90(), 0, 0, 1000, 54372, 750)
        self._basic_access_tests('SKCLIMATOLOGY_AEROSOL_CM3', climatology)
        self._basic_access_tests('SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS', climatology)
        self._basic_access_tests('SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH', climatology)

    @unittest.skip
    def test_glossac(self):
        geom0 = sk.VerticalImage()
        geom0.from_sza_saa(0, 0, 0, 0, tanalts_km=[10, 20], mjd=54372, locallook=180)

        atmo = sk.Atmosphere()
        atmo['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
        atmo['O3'] = sk.Species(sk.O3DBM(), sk.Labow())
        atmo['aer'] = sk.SpeciesAerosolGloSSAC()

        atmo.brdf = sk.Lambertian(0.4)

        eng = sk.EngineHR(geom0, atmo)
        eng.num_orders_of_scatter = 1
        eng.wavelengths = [350]
        eng.calculate_radiance()

    @unittest.skip
    def test_baum(self):
        baum = sk.SpeciesBaumIceCloud(particlesize_microns=50, cloud_top_m=1000, cloud_width_fwhm_m=200,
                                      vertical_optical_depth=0.01, vertical_optical_depth_wavel_nm=750)

        numden = baum.climatology.get_parameter('icecloud', 0, 0, baum.altitude_grid, 54732)
        xs = baum.optical_property.calculate_cross_sections(sk.MSIS90(), 0, 0, 0, 54372, 750).total

        self.assertAlmostEqual(0.01, np.sum(numden * xs) * (10 * 100))  # numden * xs will be ext in /cm
