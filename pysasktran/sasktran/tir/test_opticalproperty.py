import unittest
import sasktran as sk
import numpy as np
from sasktran.tir.opticalproperty import HITRANChemicalTIR


class TestOpticalPropertyTIR(unittest.TestCase):
    def _basic_cross_section_tests(self, optprop: sk.OpticalProperty, wavelength=np.array([10000, 10001])):
        optprop.calculate_cross_sections(sk.MSIS90(), 0, 0, 1000, 54372, wavelength)

    @unittest.skip
    def test_hitran_tir(self):
        hitran_co2 = HITRANChemicalTIR('CO2')

        self._basic_cross_section_tests(hitran_co2)

        hitran_co2 = HITRANChemicalTIR('CO2', micro_window_margin=100)
        test_wavelengths = np.array([10000, 10001])
        atmospheric_state = sk.MSIS90()
        xs = hitran_co2.calculate_cross_sections(atmospheric_state, latitude=0, longitude=0, altitude=20000, mjd=54372,
                                                 wavelengths=test_wavelengths)

        np.testing.assert_array_almost_equal_nulp(xs.scattering, np.array([0.0, 0.0]), nulp=100)
        np.testing.assert_array_almost_equal_nulp(xs.absorption, np.array([2.0644486002134603e-28,
                                                                           3.8329826741978217e-28]), nulp=100)
        np.testing.assert_array_almost_equal_nulp(xs.total, np.array([2.0644486002134603e-28,
                                                                      3.8329826741978217e-28]), nulp=100)
