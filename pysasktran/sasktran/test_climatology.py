import unittest
import sasktran as sk
import numpy as np


class TestClimatology(unittest.TestCase):
    def setUp(self):
        self.test_altitudes = np.linspace(500, 99500, 100)

    def _basic_access_tests(self, species, clim: sk.Climatology):
        clim.get_parameter(species, 0, 0, self.test_altitudes, 54372)

    def test_msis90(self):
        climatology = sk.MSIS90()

        self._basic_access_tests('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', climatology)

        test_altitudes = [10000, 30000, 50000]

        air_dens = climatology.get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 0, 0, test_altitudes, 54372)
        pressure = climatology.get_parameter('SKCLIMATOLOGY_PRESSURE_PA', 0, 0, test_altitudes, 54372)
        temperature = climatology.get_parameter('SKCLIMATOLOGY_TEMPERATURE_K', 0, 0, test_altitudes, 54372)

        np.testing.assert_array_almost_equal_nulp(air_dens, np.array([8.65527009584726e+18,
                                                                      3.7683737787462624e+17,
                                                                      2.31561579580648e+16]), nulp=100)

        np.testing.assert_array_almost_equal_nulp(pressure, np.array([28431.822158652583,
                                                                      1199.0148598042165,
                                                                      85.68265287987921]), nulp=100)

        np.testing.assert_array_almost_equal_nulp(temperature, np.array([237.92516401048934,
                                                                         230.45540624725945,
                                                                         268.0048846289609]), nulp=100)

    def test_ecmwf(self):
        climatology = sk.ECMWF()

        try:
            self._basic_access_tests('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', climatology)
        except sk.SasktranError:
            unittest.skip('ECMWF Directory could not be found, check registry settings')

    def test_labow(self):
        climatology = sk.Labow()

        self._basic_access_tests('SKCLIMATOLOGY_O3_CM3', climatology)

    @unittest.skip
    def test_pratmo(self):
        climatology = sk.Pratmo()

        self._basic_access_tests('SKCLIMATOLOGY_NO2_CM3', climatology)

    def test_lognormal_aerosol(self):
        values = np.ones_like(self.test_altitudes)
        climatology = sk.ClimatologyUserDefined(self.test_altitudes,
                                                {'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': values * 0.08,
                                                 'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': values * 1.6})

        self._basic_access_tests('SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS', climatology)
        self._basic_access_tests('SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH', climatology)

    def test_userdefined_reconfigure(self):
        climatology = sk.ClimatologyUserDefined(self.test_altitudes, {'SKCLIMATOLOGY_O3_CM3': np.ones_like(self.test_altitudes)})

        self.assertAlmostEqual(1, climatology.get_parameter('SKCLIMATOLOGY_O3_CM3', 0, 0, 1000, 54372))

        climatology['SKCLIMATOLOGY_O3_CM3'] *= 2

        self.assertAlmostEqual(2, climatology.get_parameter('SKCLIMATOLOGY_O3_CM3', 0, 0, 1000, 54372))

        climatology['SKCLIMATOLOGY_O3_CM3'][0] += 10

        self.assertAlmostEqual(12, climatology.get_parameter('SKCLIMATOLOGY_O3_CM3', 0, 0, self.test_altitudes[0], 54372)[0])

    def test_userdefined_arbitrary_guid(self):
        climatology = sk.ClimatologyUserDefined(self.test_altitudes, {'ozone': np.ones_like(self.test_altitudes)})

        self.assertAlmostEqual(1, climatology.get_parameter('ozone', 0, 0, 1000, 54372))

    def test_userdefined_twodim(self):
        # Do some geometry calculations to set up the climatology
        geo = sk.Geodetic()

        geo.from_lat_lon_alt(10, 10, 0)
        reference = geo.location

        geo.from_lat_lon_alt(10, 20, 0)
        temp = geo.location

        normal = np.cross(temp, reference)
        normal /= np.linalg.norm(normal)
        reference /= np.linalg.norm(reference)

        angle_grid = np.linspace(-10, 10, 21)
        alt_grid = np.linspace(0, 100000, 101)

        values = np.ones((len(angle_grid), len(alt_grid)))

        values[:10, :] = 10

        climatology = sk.ClimatologyUserDefined2D(angle_grid, alt_grid, {'ozone': values}, reference, normal)

        # We have made a plane climatology where the angular dimension is almost entirely in longitude, values less
        # than 10 longitude should be 10, values greater than 10 longitude should be 1

        val_1 = climatology.get_parameter('ozone', 10, 8, 10000, 54372)
        val_2 = climatology.get_parameter('ozone', 10, 12, 10000, 54372)

        self.assertAlmostEqual(val_1, 10)
        self.assertAlmostEqual(val_2, 1)
