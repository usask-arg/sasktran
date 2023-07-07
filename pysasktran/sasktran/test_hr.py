import unittest
import sasktran as sk
import numpy as np


class TestEngineHR(unittest.TestCase):
    def test_format_xarray(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

        atmosphere.wf_species = 'ozone'

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [340, 600]

        # And do the calculation
        radiance = engine.calculate_radiance('xarray')
        correct_check1 = [0.0583454507287009, 0.0575745091193133, 0.0304976942923949, 0.0100584374020211]
        correct_check2 = [0.0312688589955404, 0.0090109076095719, 0.0027612316622947, 0.0009875801095094]

        np.testing.assert_allclose(correct_check1, radiance.isel(wavelength=0, los=[0, 10, 20, 30])['radiance'].values, rtol=1e-6)
        np.testing.assert_allclose(correct_check2, radiance.isel(wavelength=1, los=[0, 10, 20, 30])['radiance'].values, rtol=1e-6)

    @unittest.skip
    def test_with_baum(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
        atmosphere['cloud'] = sk.SpeciesBaumIceCloud(50, 15000, 1000, 0.01, 750)

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1
        
        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [750]

        # And do the calculation
        radiance = engine.calculate_radiance('xarray')
        correct_check = np.array([0.03388033, 0.00774329, 0.00169492, 0.00042484])

        np.testing.assert_allclose(correct_check, radiance.isel(wavelength=0, los=[0, 10, 20, 30])['radiance'].values, rtol=1e-5)

    def test_format_xarray_scalar(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

        atmosphere.wf_species = 'ozone'

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = 340

        # And do the calculation
        radiance = engine.calculate_radiance('xarray')

    def test_polarized(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

        atmosphere.wf_species = 'ozone'

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [340, 600]

        # And do the calculation
        radiance = engine.calculate_radiance('xarray', full_stokes_vector=True)

    def test_disable_los_scattering(self):
        # First recreate our geometry and atmosphere classes
        geometry = sk.NadirGeometry()
        geometry.from_zeniths_and_azimuths(60, 30, 54372, 0, 0, reference_point=(30, 0, 0, 54372))

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1
        engine.disable_los_single_scattering = True

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [340, 600]

        # And do the calculation
        radiance = engine.calculate_radiance('xarray')

    def test_prefill_singlescatter(self):
        szas = [0, 25, 50, 75, 85, 88]
        saas = [0, 90, 180]

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

        wavel = [340, 600]

        for sza in szas:
            for saa in saas:
                geo = sk.VerticalImage()
                geo.from_sza_saa(sza, saa, 0, 0, [10, 20, 30, 40, 50], 54372, 0)

                engine_normal = sk.EngineHR(geometry=geo, atmosphere=atmosphere, wavelengths=wavel)
                engine_normal.num_orders_of_scatter = 1
                engine_normal.grid_spacing = 1000
                engine_normal.use_solartable_for_singlescatter = True

                engine_prefill = sk.EngineHR(geometry=geo, atmosphere=atmosphere, wavelengths=wavel)
                engine_prefill.use_solartable_for_singlescatter = True
                engine_prefill.prefill_solartransmission_table = True
                engine_prefill.num_orders_of_scatter = 1
                engine_prefill.grid_spacing = 1000
                # engine_prefill.options['numsolartablesza'] = 3600

                radiance_normal = engine_normal.calculate_radiance()
                radiance_prefill = engine_prefill.calculate_radiance()

                self.failUnless(np.allclose(radiance_normal, radiance_prefill, rtol=2e-3, atol=0), "Prefill table radiance significantly different than normal radiance for SZA={}, SAA={}".format(sza, saa))

    def test_hitran_hrssapprox(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
        alts = np.arange(0, 100001, 1000)
        co2 = sk.MSIS90().get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 0, 0, alts, 54372) * 400e-6

        atmosphere['co2'] = sk.Species(sk.HITRANChemical('CO2'), sk.ClimatologyUserDefined(alts, {'co2': co2}))

        engine = sk.EngineHRSSApprox(geometry, atmosphere, [1600, 1601, 1602], ms_wavelengths=[1600, 1602])

        # And do the calculation
        radiance = engine.calculate_radiance('xarray')
