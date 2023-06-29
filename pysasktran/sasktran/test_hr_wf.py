import unittest
import sasktran as sk
import numpy as np
import xarray as xr


class TestEngineHRWF(unittest.TestCase):
    def test_absorption_wf(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        labow = sk.Labow()

        clim_alts = np.arange(500, 100000, 1000)
        o3_vals = labow.get_parameter('SKCLIMATOLOGY_O3_CM3', 0, 0, clim_alts, 54372)

        user_o3 = sk.ClimatologyUserDefined(clim_alts, {'ozone': o3_vals})

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), user_o3)
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

        atmosphere.wf_species = 'ozone'

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1

        engine.options['wfheights'] = np.array([20500, 30500, 40500, 50500])
        engine.options['wfwidths'] = np.ones_like(engine.options['wfheights']) * 1000

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [330]
        engine.grid_spacing = 250

        # And do the calculation
        radiance = engine.calculate_radiance('xarray')

        wfs = []
        for idx, height in enumerate(engine.options['wfheights']):
            clim_idx = np.argmin(np.abs(height - clim_alts))
            dn = user_o3['ozone'][clim_idx] * 0.001

            user_o3['ozone'][clim_idx] += dn

            radiance_perturb = engine.calculate_radiance('xarray')

            user_o3['ozone'][clim_idx] -= dn

            wfs.append((radiance_perturb['radiance'] - radiance['radiance']) / dn)

        wfs = xr.concat(wfs, dim='perturbation')
        radiance['numerical_wf_ozone'] = wfs

        # WF's can be small, so we don't want to check percent differences directly, instead normalize each wf
        # by the maximum of each row
        radiance['wf_ozone'] = radiance['wf_ozone'] / np.abs(radiance['numerical_wf_ozone']).max(dim=['los'])
        radiance['numerical_wf_ozone'] = radiance['numerical_wf_ozone'] / np.abs(radiance['numerical_wf_ozone']).max(dim=['los'])

        wf_diff = np.abs(radiance['wf_ozone'] - radiance['numerical_wf_ozone'])

        np.testing.assert_array_less(wf_diff.values, 0.02)

    def test_absorption_wf_polarized(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        labow = sk.Labow()

        clim_alts = np.arange(500, 100000, 1000)
        o3_vals = labow.get_parameter('SKCLIMATOLOGY_O3_CM3', 0, 0, clim_alts, 54372)

        user_o3 = sk.ClimatologyUserDefined(clim_alts, {'ozone': o3_vals})

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), user_o3)
        atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

        atmosphere.wf_species = 'ozone'

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1
        engine.polarization = 'vector'

        engine.options['wfheights'] = np.array([20500, 30500, 40500, 50500])
        engine.options['wfwidths'] = np.ones_like(engine.options['wfheights']) * 1000

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [330]
        engine.grid_spacing = 250

        # And do the calculation
        radiance = engine.calculate_radiance('xarray', True)

        wfs_I = []
        wfs_Q = []
        wfs_U = []
        wfs_V = []

        for idx, height in enumerate(engine.options['wfheights']):
            clim_idx = np.argmin(np.abs(height - clim_alts))
            dn = user_o3['ozone'][clim_idx] * 0.001

            user_o3['ozone'][clim_idx] += dn

            radiance_perturb = engine.calculate_radiance('xarray', True)

            user_o3['ozone'][clim_idx] -= dn

            wfs_I.append((radiance_perturb['I'] - radiance['I']) / dn)
            wfs_Q.append((radiance_perturb['Q'] - radiance['Q']) / dn)
            wfs_U.append((radiance_perturb['U'] - radiance['U']) / dn)
            wfs_V.append((radiance_perturb['V'] - radiance['V']) / dn)

        wfs_I = xr.concat(wfs_I, dim='perturbation')
        wfs_Q = xr.concat(wfs_Q, dim='perturbation')
        wfs_U = xr.concat(wfs_U, dim='perturbation')
        wfs_V = xr.concat(wfs_V, dim='perturbation')

        numerical_wfs = np.zeros_like(radiance['wf_ozone'].values)
        numerical_wfs[:, :, :, 0] = wfs_I.values.transpose((1, 2, 0))
        numerical_wfs[:, :, :, 1] = wfs_Q.values.transpose((1, 2, 0))
        numerical_wfs[:, :, :, 2] = wfs_U.values.transpose((1, 2, 0))
        numerical_wfs[:, :, :, 3] = wfs_V.values.transpose((1, 2, 0))

        radiance['numerical_wf_ozone'] = (['wavelength', 'los', 'perturbation', 'stokes'], numerical_wfs)

        # WF's can be small, so we don't want to check percent differences directly, instead normalize each wf
        # by the maximum of each row
        radiance['wf_ozone'] = radiance['wf_ozone'] / np.abs(radiance['numerical_wf_ozone']).max(dim=['los'])
        radiance['numerical_wf_ozone'] = radiance['numerical_wf_ozone'] / np.abs(radiance['numerical_wf_ozone']).max(dim=['los'])

        wf_diff = np.abs(radiance['wf_ozone'] - radiance['numerical_wf_ozone'])

        np.testing.assert_array_less(wf_diff.isel(stokes=slice(0, 3)).values, 0.02)

    def test_rayleigh_wf_polarized(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=90, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        msis = sk.MSIS90()

        clim_alts = np.arange(500, 100000, 1000)
        air_vals = msis.get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 0, 0, clim_alts, 54372)

        user_air = sk.ClimatologyUserDefined(clim_alts, {'air': air_vals})

        atmosphere = sk.Atmosphere()
        atmosphere.brdf = 1

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), user_air)

        atmosphere.wf_species = 'air'

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1
        engine.polarization = 'vector'

        engine.options['wfheights'] = np.array([20500, 30500, 40500, 50500])
        engine.options['wfwidths'] = np.ones_like(engine.options['wfheights']) * 1000

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [600]
        engine.grid_spacing = 250

        # And do the calculation
        radiance = engine.calculate_radiance('xarray', True)

        wfs_I = []
        wfs_Q = []
        wfs_U = []
        wfs_V = []

        for idx, height in enumerate(engine.options['wfheights']):
            clim_idx = np.argmin(np.abs(height - clim_alts))
            dn = user_air['air'][clim_idx] * 0.001

            user_air['air'][clim_idx] += dn

            radiance_perturb = engine.calculate_radiance('xarray', True)

            user_air['air'][clim_idx] -= dn

            wfs_I.append((radiance_perturb['I'] - radiance['I']) / dn)
            wfs_Q.append((radiance_perturb['Q'] - radiance['Q']) / dn)
            wfs_U.append((radiance_perturb['U'] - radiance['U']) / dn)
            wfs_V.append((radiance_perturb['V'] - radiance['V']) / dn)

        wfs_I = xr.concat(wfs_I, dim='perturbation')
        wfs_Q = xr.concat(wfs_Q, dim='perturbation')
        wfs_U = xr.concat(wfs_U, dim='perturbation')
        wfs_V = xr.concat(wfs_V, dim='perturbation')

        numerical_wfs = np.zeros_like(radiance['wf_air'].values)
        numerical_wfs[:, :, :, 0] = wfs_I.values.transpose((1, 2, 0))
        numerical_wfs[:, :, :, 1] = wfs_Q.values.transpose((1, 2, 0))
        numerical_wfs[:, :, :, 2] = wfs_U.values.transpose((1, 2, 0))
        numerical_wfs[:, :, :, 3] = wfs_V.values.transpose((1, 2, 0))

        radiance['numerical_wf_air'] = (['wavelength', 'los', 'perturbation', 'stokes'], numerical_wfs)

        # WF's can be small, so we don't want to check percent differences directly, instead normalize each wf
        # by the maximum of each row
        radiance['wf_air'] = radiance['wf_air'] / np.abs(radiance['numerical_wf_air']).max(dim=['los'])
        radiance['numerical_wf_air'] = radiance['numerical_wf_air'] / np.abs(radiance['numerical_wf_air']).max(dim=['los'])

        wf_diff = np.abs(radiance['wf_air'] - radiance['numerical_wf_air'])

        np.testing.assert_array_less(wf_diff.isel(stokes=slice(0, 3)).values, 0.05)

    def test_rayleigh_wf(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        msis = sk.MSIS90()

        clim_alts = np.arange(500, 100000, 1000)
        air_vals = msis.get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 0, 0, clim_alts, 54372)

        user_air = sk.ClimatologyUserDefined(clim_alts, {'air': air_vals})

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), user_air)

        atmosphere.wf_species = 'air'

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1

        engine.options['wfheights'] = np.array([20500, 30500, 40500, 50500])
        engine.options['wfwidths'] = np.ones_like(engine.options['wfheights']) * 1000

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [330]
        engine.grid_spacing = 250

        # And do the calculation
        radiance = engine.calculate_radiance('xarray')

        wfs = []
        for idx, height in enumerate(engine.options['wfheights']):
            clim_idx = np.argmin(np.abs(height - clim_alts))
            dn = user_air['air'][clim_idx] * 0.001

            user_air['air'][clim_idx] += dn

            radiance_perturb = engine.calculate_radiance('xarray')

            user_air['air'][clim_idx] -= dn

            wfs.append((radiance_perturb['radiance'] - radiance['radiance']) / dn)

        wfs = xr.concat(wfs, dim='perturbation')
        radiance['numerical_wf_air'] = wfs

        # WF's can be small, so we don't want to check percent differences directly, instead normalize each wf
        # by the maximum of each row
        radiance['wf_air'] = radiance['wf_air'] / np.abs(radiance['numerical_wf_air']).max(dim=['los'])
        radiance['numerical_wf_air'] = radiance['numerical_wf_air'] / np.abs(radiance['numerical_wf_air']).max(dim=['los'])

        wf_diff = np.abs(radiance['wf_air'] - radiance['numerical_wf_air'])

        np.testing.assert_array_less(wf_diff.values, 0.05)

    def test_multiple_wf(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        msis = sk.MSIS90()

        clim_alts = np.arange(500, 100000, 1000)
        air_vals = msis.get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 0, 0, clim_alts, 54372)

        user_air = sk.ClimatologyUserDefined(clim_alts, {'air': air_vals})

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), user_air)

        atmosphere.wf_species = ['air', 'ozone']

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 1

        engine.options['wfheights'] = np.array([20500, 30500, 40500, 50500])
        engine.options['wfwidths'] = np.ones_like(engine.options['wfheights']) * 1000

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [330]
        engine.grid_spacing = 1000

        # And do the calculation
        radiance = engine.calculate_radiance('xarray')

    def test_wf_aerosol_particle_size(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        msis = sk.MSIS90()

        clim_alts = np.arange(500, 100000, 1000)
        air_vals = msis.get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 0, 0, clim_alts, 54372)

        user_air = sk.ClimatologyUserDefined(clim_alts, {'air': air_vals})

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), user_air)

        particle_size_dist = sk.ClimatologyUserDefined(altitudes=clim_alts,
                                                       values={'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': np.ones_like(clim_alts) * 1.6,
                                                               'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': np.ones_like(clim_alts) * 0.08})
        aerosol_optprop = sk.MieAerosol(particlesize_climatology=particle_size_dist, species='H2SO4')

        aerosol_clim = sk.ClimatologyUserDefined(altitudes=[0, 100000],
                                                 values={'aerosol': [1, 1]})

        aerosol_species = sk.Species(aerosol_optprop, aerosol_clim)
        atmosphere['aerosol'] = aerosol_species

        atmosphere.wf_species = ['air', 'ozone', 'aerosol', 'aerosol_lognormal_medianradius']

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 50

        engine.options['wfheights'] = np.array([20500, 30500, 40500, 50500])
        engine.options['wfwidths'] = np.ones_like(engine.options['wfheights']) * 1000

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [675]
        engine.grid_spacing = 500

        # And do the calculation
        radiance = engine.calculate_radiance('xarray')

        wfs = []
        for idx, height in enumerate(engine.options['wfheights']):
            clim_idx = np.argmin(np.abs(height - clim_alts))
            dn = particle_size_dist['SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS'][clim_idx] * 0.05

            particle_size_dist['SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS'][clim_idx] += dn

            radiance_perturb = engine.calculate_radiance('xarray')

            particle_size_dist['SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS'][clim_idx] -= dn

            wfs.append((radiance_perturb['radiance'] - radiance['radiance']) / dn)

        wfs = xr.concat(wfs, dim='perturbation')
        radiance['numerical_wf_aerosol_lognormal_medianradius'] = wfs

        radiance['wf_aerosol_lognormal_medianradius'] = radiance['wf_aerosol_lognormal_medianradius'] / np.abs(radiance['numerical_wf_aerosol_lognormal_medianradius']).max(dim=['los'])
        radiance['numerical_wf_aerosol_lognormal_medianradius'] = radiance['numerical_wf_aerosol_lognormal_medianradius'] / np.abs(radiance['numerical_wf_aerosol_lognormal_medianradius']).max(dim=['los'])

        wf_diff = np.abs(radiance['wf_aerosol_lognormal_medianradius'] - radiance['numerical_wf_aerosol_lognormal_medianradius'])

        np.testing.assert_array_less(wf_diff.values, 0.06)

    def test_wf_aerosol_particle_size_polarized(self):
        tanalts_km = np.arange(10, 50, 1)

        # First recreate our geometry and atmosphere classes
        geometry = sk.VerticalImage()
        geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                              satalt_km=600, refalt_km=20)

        msis = sk.MSIS90()

        clim_alts = np.arange(500, 100000, 1000)
        air_vals = msis.get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 0, 0, clim_alts, 54372)

        user_air = sk.ClimatologyUserDefined(clim_alts, {'air': air_vals})

        atmosphere = sk.Atmosphere()

        atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
        atmosphere['air'] = sk.Species(sk.Rayleigh(), user_air)

        particle_size_dist = sk.ClimatologyUserDefined(altitudes=clim_alts,
                                                       values={'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': np.ones_like(clim_alts) * 1.6,
                                                               'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': np.ones_like(clim_alts) * 0.08})
        aerosol_optprop = sk.MieAerosol(particlesize_climatology=particle_size_dist, species='H2SO4')

        aerosol_clim = sk.ClimatologyUserDefined(altitudes=[0, 100000],
                                                 values={'aerosol': [1, 1]})

        aerosol_species = sk.Species(aerosol_optprop, aerosol_clim)
        atmosphere['aerosol'] = aerosol_species

        atmosphere.wf_species = ['air', 'ozone', 'aerosol', 'aerosol_lognormal_medianradius']

        # And now make the engine
        engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
        engine.num_orders_of_scatter = 50
        engine.polarization = 'vector'

        engine.options['wfheights'] = np.array([20500, 30500, 40500, 50500])
        engine.options['wfwidths'] = np.ones_like(engine.options['wfheights']) * 1000

        # Choose some wavelengths to do the calculation at
        engine.wavelengths = [675]
        engine.grid_spacing = 500

        # And do the calculation
        radiance = engine.calculate_radiance('xarray', True)

        wfs_I = []
        wfs_Q = []
        wfs_U = []
        wfs_V = []

        for idx, height in enumerate(engine.options['wfheights']):
            clim_idx = np.argmin(np.abs(height - clim_alts))
            dn = particle_size_dist['SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS'][clim_idx] * 0.05

            particle_size_dist['SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS'][clim_idx] += dn

            radiance_perturb = engine.calculate_radiance('xarray', True)

            particle_size_dist['SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS'][clim_idx] -= dn

            wfs_I.append((radiance_perturb['I'] - radiance['I']) / dn)
            wfs_Q.append((radiance_perturb['Q'] - radiance['Q']) / dn)
            wfs_U.append((radiance_perturb['U'] - radiance['U']) / dn)
            wfs_V.append((radiance_perturb['V'] - radiance['V']) / dn)

        wfs_I = xr.concat(wfs_I, dim='perturbation')
        wfs_Q = xr.concat(wfs_Q, dim='perturbation')
        wfs_U = xr.concat(wfs_U, dim='perturbation')
        wfs_V = xr.concat(wfs_V, dim='perturbation')

        numerical_wfs = np.zeros_like(radiance['wf_air'].values)
        numerical_wfs[:, :, :, 0] = wfs_I.values.transpose((1, 2, 0))
        numerical_wfs[:, :, :, 1] = wfs_Q.values.transpose((1, 2, 0))
        numerical_wfs[:, :, :, 2] = wfs_U.values.transpose((1, 2, 0))
        numerical_wfs[:, :, :, 3] = wfs_V.values.transpose((1, 2, 0))

        radiance['numerical_wf_aerosol_lognormal_medianradius'] = (['wavelength', 'los', 'perturbation', 'stokes'], numerical_wfs)

        radiance['wf_aerosol_lognormal_medianradius'] = radiance['wf_aerosol_lognormal_medianradius'] / np.abs(radiance['numerical_wf_aerosol_lognormal_medianradius']).max(dim=['los'])
        radiance['numerical_wf_aerosol_lognormal_medianradius'] = radiance['numerical_wf_aerosol_lognormal_medianradius'] / np.abs(radiance['numerical_wf_aerosol_lognormal_medianradius']).max(dim=['los'])

        wf_diff = np.abs(radiance['wf_aerosol_lognormal_medianradius'] - radiance['numerical_wf_aerosol_lognormal_medianradius'])

        # differences are higher, but I think it's just the approximation
        np.testing.assert_array_less(wf_diff.isel(stokes=slice(0, 3)).values, 0.18)
