import unittest
import sasktran as sk
import numpy as np
from sasktran.tir.climatology import ClimatologySpeciesCustomPT, ClimatologySpecies, ClimatologyAtmosphericState,\
    ClimatologyFull


class TestClimatologyTIR(unittest.TestCase):
    def setUp(self):
        self.test_altitudes = np.linspace(500, 99500, 100)

    def _basic_access_tests(self, species, clim: sk.Climatology):
        clim.get_parameter(species, 0, 0, self.test_altitudes, 54372)

    def test_atmospheric_state(self):
        climatology = ClimatologyAtmosphericState()

        self._basic_access_tests('SKCLIMATOLOGY_TEMPERATURE_K', climatology)

        test_altitudes = [10000, 30000, 50000]

        pressure = climatology.get_parameter('SKCLIMATOLOGY_PRESSURE_PA', 0, 0, test_altitudes, 54372)
        temperature = climatology.get_parameter('SKCLIMATOLOGY_TEMPERATURE_K', 0, 0, test_altitudes, 54372)

        np.testing.assert_array_almost_equal_nulp(pressure, np.array([26500.0, 1197.0, 79.78]), nulp=100)
        np.testing.assert_array_almost_equal_nulp(temperature, np.array([223.3, 226.5, 270.7]), nulp=100)

    def test_climatology_species(self):
        climatology = ClimatologySpecies('CO2')

        self._basic_access_tests('SKCLIMATOLOGY_CO2_CM3', climatology)

        test_altitudes = [10000, 30000, 50000]

        co2_dens = climatology.get_parameter('SKCLIMATOLOGY_CO2_CM3', 0, 0, test_altitudes, 54372)

        np.testing.assert_array_almost_equal_nulp(co2_dens, np.array([2836532788301347.5,
                                                                      126315487132049.84,
                                                                      7044277132027.684]), nulp=100)

    def test_climatology_full(self):
        climatology = ClimatologyFull()

        self._basic_access_tests('SKCLIMATOLOGY_CO2_CM3', climatology)

        test_altitudes = [10000, 30000, 50000]

        pressure = climatology.get_parameter('SKCLIMATOLOGY_PRESSURE_PA', 0, 0, test_altitudes, 54372)
        temperature = climatology.get_parameter('SKCLIMATOLOGY_TEMPERATURE_K', 0, 0, test_altitudes, 54372)
        co2_dens = climatology.get_parameter('SKCLIMATOLOGY_CO2_CM3', 0, 0, test_altitudes, 54372)

        np.testing.assert_array_almost_equal_nulp(pressure, np.array([26500.0, 1197.0, 79.78]), nulp=100)
        np.testing.assert_array_almost_equal_nulp(temperature, np.array([223.3, 226.5, 270.7]), nulp=100)
        np.testing.assert_array_almost_equal_nulp(co2_dens, np.array([2836532788301347.5,
                                                                      126315487132049.84,
                                                                      7044277132027.684]), nulp=100)

    def test_climatology_custom_pt(self):
        atmospheric_state = sk.MSIS90()

        climatology = ClimatologySpeciesCustomPT('CO2', atmospheric_state, latitude=0, longitude=0, mjd=54372)

        self._basic_access_tests('SKCLIMATOLOGY_CO2_CM3', climatology)

        test_altitudes = [10000, 30000, 50000]

        co2_dens = climatology.get_parameter('SKCLIMATOLOGY_CO2_CM3', 0, 0, test_altitudes, 54372)

        np.testing.assert_array_almost_equal_nulp(co2_dens, np.array([2856239131629593.0,
                                                                      124356334698626.53,
                                                                      7641532126161.376]), nulp=100)

    def test_fascode_dataset(self):
        dataset = 'fascode'
        climatology_ids = ['mls', 'mlw', 'sas', 'saw', 'std', 'tro']
        species_major = ['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2']
        species_minor = ['NO', 'SO2', 'NO2', 'NH3', 'HNO3', 'OH', 'HF', 'HCl', 'HBr', 'HI', 'ClO', 'OCS', 'H2CO',
                         'HOCl', 'N2', 'HCN', 'CH3Cl', 'H2O2', 'C2H2', 'C2H6', 'PH3', 'COF2', 'SF6', 'H2S', 'CFCl3',
                         'CF2Cl2', 'CClF3', 'CF4', 'CHCl2F', 'CHClF2', 'C2Cl3F3', 'C2Cl2F4', 'C2ClF5', 'CCl4',
                         'ClONO2', 'N2O5', 'HNO4']

        for clim in climatology_ids:
            for species in species_major + species_minor:
                climatology = ClimatologySpecies(species, dataset=dataset, climatology=clim)
                self._basic_access_tests('SKCLIMATOLOGY_' + species + '_CM3', climatology)

    def test_mipas_1998_dataset(self):
        dataset = 'mipas_1998'
        climatology_ids = ['day_imk', 'ngt_imk', 'sum_imk', 'win_imk']
        species_major = ['N2', 'O2', 'O3P', 'CO2', 'O3', 'H2O', 'CH4', 'N2O', 'HNO3', 'CO', 'NO2', 'N2O5', 'ClO',
                         'HOCl', 'ClONO2', 'NO']
        species_minor = ['HCN', 'H2O2', 'F12', 'F14', 'F22', 'COF2', 'OCS', 'NH3', 'SO2', 'CFCl3', 'C2H2', 'C2H6',
                         'CCl4', 'SF6', 'HNO4', 'CH3Cl', 'CClF3', 'CHCl2F', 'C2Cl3F3', 'C2Cl2F4']

        for clim in climatology_ids:
            for species in species_major + species_minor:
                climatology = ClimatologySpecies(species, dataset=dataset, climatology=clim)
                self._basic_access_tests('SKCLIMATOLOGY_' + species + '_CM3', climatology)

    def test_mipas_2001_dataset(self):
        dataset = 'mipas_2001'
        climatology_ids = ['day', 'equ', 'ngt', 'sum', 'win']
        species_major = ['N2', 'O2', 'CO2', 'O3', 'H2O', 'CH4', 'N2O', 'HNO3', 'CO', 'NO2', 'N2O5', 'ClO', 'HOCl',
                         'ClONO2', 'NO', 'HNO4', 'HCN', 'NH3', 'F11', 'F12', 'F14', 'F22', 'CCl4', 'COF2', 'H2O2',
                         'C2H2', 'C2H6', 'OCS', 'SO2', 'SF6']
        species_minor = ['CClF3', 'CHCl2F', 'C2Cl3F3', 'C2Cl2F4', 'C2ClF5', 'CH3Cl', 'H2S']

        for clim in climatology_ids:
            for species in species_major + species_minor:
                climatology = ClimatologySpecies(species, dataset=dataset, climatology=clim)
                self._basic_access_tests('SKCLIMATOLOGY_' + species + '_CM3', climatology)
