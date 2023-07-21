import sasktran as sk

from sasktran.disco.tests.common import default_atmosphere, default_geometry, standard_test_wavelengths
import sasktran.disco.interface as do
import numpy as np


NUM_WF_CALC = 100


def numerical_wf(engine: do.EngineDO, climatology: sk.ClimatologyUserDefined, species_name: str, rel_eps=0.001):
    altitudes = climatology.altitudes

    start_rad = engine.calculate_radiance(output_format='xarray')

    numerical_wf = np.zeros((len(start_rad.wavelength), len(start_rad.los), len(altitudes), 3))

    for idx in range(NUM_WF_CALC):
        dx = climatology[species_name][idx] * rel_eps
        climatology[species_name][idx] += dx

        pert_rad_above = engine.calculate_radiance(output_format='xarray')
        climatology[species_name][idx] -= 2*dx
        pert_rad_below = engine.calculate_radiance(output_format='xarray')
        climatology[species_name][idx] += dx

        numerical_wf[:, :, idx, :] = ((pert_rad_above['radiance'].values - pert_rad_below['radiance'].values) / (2*dx))

    start_rad['wf_numerical_' + species_name] = (['wavelength', 'los', 'perturbation', 'stokes'], numerical_wf)

    return start_rad


def test_absorption_wf():
    """
    Verifies that the weighting function for ozone absorption is equal to a numerical calculation to 7 decimal places
    """
    ozone_wf_wavelengths = [310, 330, 600]

    atmo = default_atmosphere()
    geo = default_geometry()

    engine = do.EngineDO(atmosphere=atmo, geometry=geo, wavelengths=ozone_wf_wavelengths)
    engine.num_layers = 20
    engine.options['forcenumazimuthterms'] = 4
    engine.num_stokes = 3
    engine.num_streams = 4

    altitudes = atmo['o3'].climatology.altitudes

    engine.options['wfaltitudes'] = altitudes
    engine.options['wfwidths'] = np.gradient(altitudes)
    atmo.wf_species = 'o3'

    rad = numerical_wf(engine, atmo['o3'].climatology, 'o3', rel_eps=0.01)

    rad = rad.isel(perturbation=slice(0, NUM_WF_CALC))

    p_diff = (rad['wf_o3'] - rad['wf_numerical_o3']) / rad['wf_numerical_o3'] * 100

    np.testing.assert_almost_equal(p_diff.values[:, :, :, 0], 0.0, decimal=4)

    # Use a little lower precision here since the derivatives on Q/U are quite sensitive to numerical differences
    np.testing.assert_almost_equal(p_diff.values[:, :, :, 1], 0.0, decimal=3)
    np.testing.assert_almost_equal(p_diff.values[:, :, :, 2], 0.0, decimal=3)


def test_scattering_wf():
    """
    Tests that the weighting function for Rayleigh scattering is equal to the numerical one up to 7 decimal places
    """
    air_wavelengths = [330]

    atmo = default_atmosphere()
    geo = default_geometry()

    engine = do.EngineDO(atmosphere=atmo, geometry=geo, wavelengths=air_wavelengths)
    engine.num_layers = 20
    engine.options['forcenumazimuthterms'] = 4
    engine.num_streams = 4
    engine.num_stokes = 3

    altitudes = atmo['air'].climatology.altitudes

    engine.options['wfaltitudes'] = altitudes
    engine.options['wfwidths'] = np.gradient(altitudes)
    atmo.wf_species = 'air'

    rad = numerical_wf(engine, atmo['air'].climatology, 'air', rel_eps=0.01)

    rad = rad.isel(perturbation=slice(0, NUM_WF_CALC))

    p_diff = (rad['wf_air'] - rad['wf_numerical_air']) / rad['wf_numerical_air'] * 100

    np.testing.assert_almost_equal(p_diff.values[:, :, :, 0], 0.0, decimal=5)

    # Q and U are trickier since there are conditions where dQ is ~ 0
    np.testing.assert_array_almost_equal_nulp(rad['wf_air'].values[:, :, :, 1], rad['wf_numerical_air'].values[:, :, :, 1], nulp=1e13)
    np.testing.assert_array_almost_equal_nulp(rad['wf_air'].values[:, :, :, 2], rad['wf_numerical_air'].values[:, :, :, 2], nulp=1e13)


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
    engine.num_stokes = 3

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

    # Check I normally
    np.testing.assert_almost_equal(p_diff.values[:, :, 0], 0.0, decimal=4)

    # Q/U can have zeros or are really small, check with ULP
    Q_analytic = rad['wf_brdf'].values[:, :, 1].flatten()
    Q_numerical = wf.values[:, :, 1].flatten()

    U_analytic = rad['wf_brdf'].values[:, :, 2].flatten()
    U_numerical = wf.values[:, :, 2].flatten()

    Q_good = abs(Q_numerical) > 1e-10

    np.testing.assert_array_almost_equal_nulp(Q_analytic[Q_good], Q_numerical[Q_good], nulp=1e10)

    U_good = abs(U_numerical) > 1e-10

    np.testing.assert_array_almost_equal_nulp(U_analytic[U_good], U_numerical[U_good], nulp=1e10)