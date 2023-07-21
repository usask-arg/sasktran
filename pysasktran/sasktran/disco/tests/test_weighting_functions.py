import sasktran as sk

from sasktran.disco.tests.common import default_atmosphere, default_geometry, standard_test_wavelengths
import sasktran.disco.interface as do
import numpy as np


NUM_WF_CALC = 100


def numerical_wf(engine: do.EngineDO, climatology: sk.ClimatologyUserDefined, species_name: str, rel_eps=0.001):
    altitudes = climatology.altitudes

    start_rad = engine.calculate_radiance(output_format='xarray')

    numerical_wf = np.zeros((len(start_rad.wavelength), len(start_rad.los), len(altitudes)))

    for idx in range(NUM_WF_CALC):
        dx = climatology[species_name][idx] * rel_eps
        climatology[species_name][idx] += dx

        pert_rad_above = engine.calculate_radiance(output_format='xarray')
        climatology[species_name][idx] -= 2*dx
        pert_rad_below = engine.calculate_radiance(output_format='xarray')
        climatology[species_name][idx] += dx

        numerical_wf[:, :, idx] = ((pert_rad_above['radiance'].values - pert_rad_below['radiance'].values) / (2*dx))

    start_rad['wf_numerical_' + species_name] = (['wavelength', 'los', 'perturbation'], numerical_wf)

    return start_rad


def test_absorption_wf():
    """
    Verifies that the weighting function for ozone absorption is equal to a numerical calculation to 7 decimal places
    """
    ozone_wf_wavelengths = [310, 330, 600]

    atmo = default_atmosphere()
    geo = default_geometry()

    engine = do.EngineDO(atmosphere=atmo, geometry=geo, wavelengths=ozone_wf_wavelengths)
    engine.num_layers = 50
    engine.options['forcenumazimuthterms'] = 4
    engine.num_streams = 4

    altitudes = atmo['o3'].climatology.altitudes

    engine.options['wfaltitudes'] = altitudes
    engine.options['wfwidths'] = np.gradient(altitudes)
    atmo.wf_species = 'o3'

    rad = numerical_wf(engine, atmo['o3'].climatology, 'o3', rel_eps=0.01)

    rad = rad.isel(perturbation=slice(0, NUM_WF_CALC))

    p_diff = (rad['wf_o3'] - rad['wf_numerical_o3']) / rad['wf_numerical_o3'] * 100

    np.testing.assert_almost_equal(p_diff.values, 0.0, decimal=4)


def test_scattering_wf():
    """
    Tests that the weighting function for Rayleigh scattering is equal to the numerical one up to 7 decimal places
    """
    air_wavelengths = [330]

    atmo = default_atmosphere()
    geo = default_geometry()

    engine = do.EngineDO(atmosphere=atmo, geometry=geo, wavelengths=air_wavelengths)
    engine.num_layers = 50
    engine.options['forcenumazimuthterms'] = 4
    engine.num_streams = 4

    altitudes = atmo['air'].climatology.altitudes

    engine.options['wfaltitudes'] = altitudes
    engine.options['wfwidths'] = np.gradient(altitudes)
    atmo.wf_species = 'air'

    rad = numerical_wf(engine, atmo['air'].climatology, 'air', rel_eps=0.01)

    rad = rad.isel(perturbation=slice(0, NUM_WF_CALC))

    p_diff = (rad['wf_air'] - rad['wf_numerical_air']) / rad['wf_numerical_air'] * 100

    np.testing.assert_almost_equal(p_diff.values, 0.0, decimal=5)


def test_multiple_species_wf():
    """
    Tests that the weighting function works for multiple species at the same time
    """
    air_wavelengths = [330]

    atmo = default_atmosphere()

    atmo['air2'] = sk.Species(sk.SimpleRayleigh(), sk.MSIS90())

    geo = default_geometry()

    engine = do.EngineDO(atmosphere=atmo, geometry=geo, wavelengths=air_wavelengths)
    engine.num_layers = 20
    engine.options['forcenumazimuthterms'] = 1
    engine.num_streams = 4

    altitudes = atmo['air'].climatology.altitudes

    engine.options['wfaltitudes'] = altitudes
    engine.options['wfwidths'] = np.gradient(altitudes)
    atmo.wf_species = ['air']

    rad = numerical_wf(engine, atmo['air'].climatology, 'air', rel_eps=0.01)

    rad = rad.isel(perturbation=slice(0, NUM_WF_CALC))

    p_diff = (rad['wf_air'] - rad['wf_numerical_air']) / rad['wf_numerical_air'] * 100

    # This is currently bugged since we assume the legendre coefficients of the species do not change in altitude
    # but they slightly do for Rayleigh scattering, but it is still close
    np.testing.assert_almost_equal(p_diff.values, 0.0, decimal=1)


def test_brdf_wf():
    """

    """
    brdf_wavelengths = np.arange(300, 800, 1)

    atmo = default_atmosphere()
    geo = default_geometry()

    atmo.brdf = 0.3

    engine = do.EngineDO(atmosphere=atmo, geometry=geo, wavelengths=brdf_wavelengths)
    engine.num_layers = 50
    engine.options['forcenumazimuthterms'] = 4
    engine.num_streams = 4

    altitudes = atmo['air'].climatology.altitudes

    engine.options['wfaltitudes'] = altitudes
    engine.options['wfwidths'] = np.gradient(altitudes)
    atmo.wf_species = 'brdf'

    rad = engine.calculate_radiance(output_format='xarray')

    dx = 0.001
    atmo.brdf.albedo += dx
    rad_u = engine.calculate_radiance(output_format='xarray')
    atmo.brdf.albedo -= 2*dx
    rad_l = engine.calculate_radiance(output_format='xarray')

    wf = (rad_u['radiance'] - rad_l['radiance']) / (2*dx)

    p_diff = (rad['wf_brdf'] - wf) / wf * 100

    np.testing.assert_almost_equal(p_diff.values, 0.0, decimal=4)


def test_inside_atmosphere_wf():
    """
    """
    ozone_wf_wavelengths = [310]

    atmo = default_atmosphere()
    geo = default_geometry(obs_alt=30000)

    engine = do.EngineDO(atmosphere=atmo, geometry=geo, wavelengths=ozone_wf_wavelengths)
    engine.num_layers = 50
    engine.options['forcenumazimuthterms'] = 4
    engine.num_streams = 4

    altitudes = atmo['o3'].climatology.altitudes

    engine.options['wfaltitudes'] = altitudes
    engine.options['wfwidths'] = np.gradient(altitudes)
    engine.options['usepseudospherical'] = 0
    atmo.wf_species = 'o3'

    rad = numerical_wf(engine, atmo['o3'].climatology, 'o3', rel_eps=0.01)

    rad = rad.isel(perturbation=slice(0, NUM_WF_CALC))

    p_diff = (rad['wf_o3'] - rad['wf_numerical_o3']) / rad['wf_numerical_o3'] * 100

    np.testing.assert_almost_equal(p_diff.values, 0.0, decimal=4)
