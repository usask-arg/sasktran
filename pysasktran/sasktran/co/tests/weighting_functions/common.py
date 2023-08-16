import sasktran as sk
import numpy as np


def numerical_wf(engine: sk.EngineCO, climatology: sk.ClimatologyUserDefined, species_name: str, rel_eps=0.001):
    altitudes = climatology.altitudes

    start_rad = engine.calculate_radiance(output_format='xarray')

    numerical_wf = np.zeros((len(start_rad.wavelength), len(start_rad.los), len(altitudes)))

    for idx in range(len(altitudes)):
        dx = climatology[species_name][idx] * rel_eps
        climatology[species_name][idx] += dx

        pert_rad_above = engine.calculate_radiance(output_format='xarray')
        climatology[species_name][idx] -= 2*dx
        pert_rad_below = engine.calculate_radiance(output_format='xarray')
        climatology[species_name][idx] += dx

        numerical_wf[:, :, idx] = ((pert_rad_above['radiance'].values - pert_rad_below['radiance'].values) / (2*dx))

    start_rad['wf_numerical_' + species_name] = (['wavelength', 'los', 'perturbation'], numerical_wf)

    return start_rad


def validate_wf(analytic, numerical, decimal=6):

    max_by_alt = np.abs(analytic).max(dim='perturbation')

    max_by_alt.values[max_by_alt.values == 0] = 1e99

    rel_diff = (analytic - numerical) / max_by_alt

    np.testing.assert_array_almost_equal(rel_diff, 0, decimal=decimal)
