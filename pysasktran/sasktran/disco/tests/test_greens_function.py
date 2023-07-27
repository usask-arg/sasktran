import sasktran as sk

from sasktran.disco.tests.common import default_atmosphere, default_geometry, standard_test_wavelengths
import sasktran.disco.interface as do
import numpy as np


NUM_WF_CALC = 100


def numerical_wf(engine: do.EngineDO, climatology: sk.ClimatologyUserDefined, species_name: str, rel_eps=0.001):
    altitudes = climatology.altitudes

    start_rad = engine.calculate_radiance(output_format='xarray')
    engine.atmosphere.wf_species = None

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


def test_greens_function_equal():
    atmo = default_atmosphere()
    geo = default_geometry()

    atmo.wf_species = ['o3', 'air']

    wavelengths = standard_test_wavelengths()
    for num_layer in [2, 10, 50, 100]:
        for num_stream in [2, 8, 16]:
            engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
            engine.num_layers = num_layer
            engine.num_streams = num_stream

            rad_standard = engine.calculate_radiance(output_format='xarray')

            engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
            engine.num_layers = num_layer
            engine.num_streams = num_stream
            engine.options['usegreensfunction'] = 1

            rad_greens = engine.calculate_radiance(output_format='xarray')

            np.testing.assert_almost_equal(rad_standard.radiance.values, rad_greens.radiance.values)
            np.testing.assert_almost_equal(rad_standard.wf_o3.values, rad_greens.wf_o3.values)

            np.testing.assert_almost_equal(rad_standard.radiance.values, rad_greens.radiance.values)
            np.testing.assert_almost_equal(rad_standard.wf_air.values, rad_greens.wf_air.values)


def test_greens_function_losspherical_wf():
    """
    Verifies that the weighting function for ozone absorption is equal to a numerical calculation to 7 decimal places
    """
    ozone_wf_wavelengths = [310, 330, 600]

    atmo = default_atmosphere()
    geo = default_geometry()

    engine = do.EngineDO(atmosphere=atmo, geometry=geo, wavelengths=ozone_wf_wavelengths)
    engine.num_layers = 10
    engine.options['forcenumazimuthterms'] = 4
    engine.num_streams = 4

    altitudes = atmo['o3'].climatology.altitudes

    engine.options['wfaltitudes'] = altitudes
    engine.options['wfwidths'] = np.gradient(altitudes)
    engine.options['uselosspherical'] = 1
    engine.options['usegreensfunction'] = 1

    atmo.wf_species = 'o3'

    rad = numerical_wf(engine, atmo['o3'].climatology, 'o3', rel_eps=0.001)

    rad = rad.isel(perturbation=slice(0, NUM_WF_CALC))

    p_diff = (rad['wf_o3'] - rad['wf_numerical_o3']) / rad['wf_numerical_o3'] * 100

    np.testing.assert_almost_equal(p_diff.values, 0.0, decimal=4)