import unittest
import sasktran as sk
import numpy as np


class TestEmission(unittest.TestCase):

    def test_userdefined(self):
        alts = np.arange(500, 100500, 1000)
        wavel = np.linspace(300, 350, 50)
        ver = np.ones((len(alts), len(wavel)))
        emission = sk.EmissionTable(alts, wavel, ver)
        atmospheric_state = sk.MSIS90()
        radiance = emission.isotropic_emission(0, 0, 22500, 54372, [350])
        self.assertAlmostEqual(1, radiance)

    def test_hitran_emissions(self):
        msis = sk.MSIS90()
        constant = sk.ClimatologyUserDefined(np.array([0.0, 100000.0]), {'SKEMISSION_PHOTOCHEMICAL_O2': [1.0E6, 1.0E6]})
        emission = sk.HITRANPhotoChemical_O2_SingletDelta(msis, (constant, 'SKEMISSION_PHOTOCHEMICAL_O2'), self_broadening_climatology=(msis, 'SKCLIMATOLOGY_O2_CM3'))
        wavelen = np.arange(1223.0, 1321.0, 0.0001)
        signal = emission.isotropic_emission(52.0, -106.0, 80000.0, 57005.8, wavelen, False)

    def test_with_engine(self):
        from sasktran.geometry import VerticalImage

        # First recreate our geometry and atmosphere classes
        geometry = VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=[10, 20, 30, 40], mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        alts = np.arange(500, 100500, 1000)
        wavel = np.linspace(300, 350, 50)
        ver = np.ones((len(alts), len(wavel)))
        emission = sk.EmissionTable(alts, wavel, ver)

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

        atmosphere.emissions['test'] = emission

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)

        engine.options['useemissions'] = 1

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [340, 600]

        # And do the calculation
        radiance = engine.calculate_radiance()
