import sasktran as sk

from sasktran.co.tests.common import default_atmosphere, default_geometry, standard_test_wavelengths
import numpy as np


NUM_WF_CALC = 100


def numerical_wf(engine: sk.EngineCO, climatology: sk.ClimatologyUserDefined, species_name: str, rel_eps=0.001):
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

    atmo = default_atmosphere(altitude_spacing=1000)
    geo = default_geometry()

    engine = sk.EngineCO(atmosphere=atmo, geometry=geo, wavelengths=ozone_wf_wavelengths,
                         options={'numssmoments': 2,
                                  'numdostreams': 2})

    altitudes = atmo['o3'].climatology.altitudes
    atmo.wf_species = 'o3'

    rad = numerical_wf(engine, atmo['o3'].climatology, 'o3', rel_eps=0.01)

    rad = rad.isel(perturbation=slice(0, NUM_WF_CALC))

    p_diff = (rad['wf_o3'] - rad['wf_numerical_o3']) / rad['wf_numerical_o3'] * 100

    np.testing.assert_almost_equal(p_diff.values, 0.0, decimal=4)


def test_scattering_wf():
    """
    Verifies that the weighting function for ozone absorption is equal to a numerical calculation to 7 decimal places
    """
    ozone_wf_wavelengths = [310, 350, 600]

    atmo = default_atmosphere(altitude_spacing=1000)
    geo = default_geometry()

    engine = sk.EngineCO(atmosphere=atmo, geometry=geo, wavelengths=ozone_wf_wavelengths,
                         options={'numssmoments': 2,
                                  'numdostreams': 2}
                         )

    atmo.wf_species = 'air'

    rad = numerical_wf(engine, atmo['air'].climatology, 'air', rel_eps=0.01)

    rad = rad.isel(perturbation=slice(0, NUM_WF_CALC))

    p_diff = (rad['wf_air'] - rad['wf_numerical_air']) / rad['wf_numerical_air'] * 100

    np.testing.assert_almost_equal(p_diff.values, 0.0, decimal=4)