import unittest
import sasktran as sk
import numpy as np
from sys import platform


class TestEngineSO(unittest.TestCase):
    @unittest.skip
    def test_so(self):
        if platform == 'linux' or platform == 'linux2':
            self.skipTest('Skipping test on Linux')
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
        geometry.from_sza_saa(sza, saa, lat, lng, tanalts_km, mjd, rayazi, satalt_km=600, refalt_km=20)
        wavelengths = [600.0, 340.0]
        engine = sk.EngineSO(geometry=geometry, atmosphere=atmosphere, wavelengths=wavelengths)
        engine.num_orders_of_scatter = 10
        radiance = engine.calculate_radiance()
        hardcodedSO = np.array([[0.05592585348539569, 0.045392502931164834, 0.03633963097994816, 0.02884492287048933],      # Hard coded values from 2020-09-16.
                               [0.15607917070783842, 0.15527077617343513, 0.15424333937565035, 0.15258545307242732]])
        diff = np.abs((radiance - hardcodedSO) / radiance)
        maxdiff = np.max(diff)
        self.assertAlmostEqual(maxdiff, 0.0, places=5)
