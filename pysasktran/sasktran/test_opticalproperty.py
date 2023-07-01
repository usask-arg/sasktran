import unittest
import sasktran as sk
import numpy as np


class TestOpticalProperty(unittest.TestCase):
    def setUp(self):
        self.test_altitudes = np.linspace(500, 99500, 100)

    def _basic_cross_section_tests(self, optprop: sk.OpticalProperty, wavelength=750.0):
        optprop.calculate_cross_sections(sk.MSIS90(), 0, 0, 1000, 54372, wavelength)

    def test_rayleigh(self):
        rayleigh = sk.Rayleigh()
        self._basic_cross_section_tests(rayleigh)

    def test_mieaerosol(self):
        climatology = sk.ClimatologyUserDefined(self.test_altitudes,
                                                {'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS':
                                                    np.ones(self.test_altitudes.shape) * 0.08,
                                                 'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH':
                                                    np.ones(self.test_altitudes.shape) * 1.6})

        mie_aerosol = sk.MieAerosol(climatology, 'H2SO4')
        self._basic_cross_section_tests(mie_aerosol)

    def test_o3_dbm(self):
        o3_dbm = sk.O3DBM()
        self._basic_cross_section_tests(o3_dbm)

    def test_no2_vandaele(self):
        no2_van = sk.NO2Vandaele1998()
        self._basic_cross_section_tests(no2_van)

    @unittest.skip
    def test_baum(self):
        baum = sk.BaumIceCrystal(50)
        self._basic_cross_section_tests(baum)

    def test_convolved(self):
        hires_opt = sk.O3DBM()
        convolved_opt = sk.OpticalPropertyConvolved(hires_opt, np.linspace(100, 1000, 10), np.linspace(0.5, 10, 10))

        self._basic_cross_section_tests(convolved_opt)

        # Test with scalar psf
        convolved_opt = sk.OpticalPropertyConvolved(hires_opt, np.linspace(100, 1000, 10), 2)
        self._basic_cross_section_tests(convolved_opt)

        # Test with defined output spacing
        convolved_opt = sk.OpticalPropertyConvolved(hires_opt, np.linspace(100, 1000, 10), 2, output_spacing=1)
        self._basic_cross_section_tests(convolved_opt)

    def test_userdefinedpressure(self):
        pressures = np.array([0, 100000])
        temperatures = np.array([[100, 150, 300], [100, 150, 300]])
        wavel = np.array([0, 10000])
        vmr = np.array([1.0])

        xs = np.ones((len(vmr), len(pressures), temperatures.shape[1], len(wavel)))

        opt_prop = sk.UserDefinedAbsorptionPressure(vmr, pressures, temperatures, wavel, xs)

        opt_prop.calculate_cross_sections(sk.MSIS90(), 0, 0, 20000, 54372, [750])
