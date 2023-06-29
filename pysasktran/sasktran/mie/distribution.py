from scipy.stats import rv_continuous
import scipy.integrate as integrate
import sasktran as sk
import numpy as np
import xarray as xr


def sp_lognormal(median_radius: float, mode_width: float):
    """
    Creates a standard lognormal particle size distribution from median radius and mode width using the internal
    scipy lognorm distribution

    Parameters
    ----------
    median_radius : float
        Median radius of the distribution.  Can be in any units but note that the internal distribution units
        must be the same as wavelength if integrating over the distribution.
    mode_width : float
        Mode width of the distribution

    Returns
    -------
    rv_continuous
        An instance of a distribution
    """
    from scipy.stats import lognorm

    return lognorm(np.log(mode_width), scale=median_radius)


def _integrable_parameters(mie: sk.Mie, pdf):
    """

    """
    out = dict()
    s1 = mie.S1()
    s2 = mie.S2()

    # Phase function elements
    # Note that P22 = P11 and P44 = p33 so we have 6 indenpendent elements
    out['p11'] = np.abs(s1)**2 + np.abs(s2)**2
    out['p12'] = np.abs(s1)**2 - np.abs(s2)**2
    out['p33'] = np.real(s1 * np.conj(s2) + s2 * np.conj(s1))
    out['p34'] = np.real(-1j * (s1 * np.conj(s2) - s2 * np.conj(s1)))

    # Cross sections in units of input radius**2
    out['Cext'] = mie.Cext()
    out['Csca'] = mie.Csca()
    out['Cabs'] = mie.Cabs()

    # Phase moments
    out['Pmom'] = mie.Pmom()

    for key in out.keys():
        out[key] *= pdf

    return out


def _post_process(total: dict, wavelength):
    """
    Internal method to post process the Mie solution.  This involves transforming the quantities that we integrate over
    to the relevant bulk properties and applying the necessary normalizations.
    """
    k = 2 * np.pi / wavelength
    c = 4 * np.pi / (2 * k**2 * total['Csca'])

    total['p11'] *= c
    total['p12'] *= c
    total['p33'] *= c
    total['p34'] *= c

    total['Pmom'] *= (2 * np.array(list(range(0, total['Pmom'].shape[1]))) + 1)[np.newaxis, :]

    total['lm_p11'] = (total['Pmom'][0, :] + total['Pmom'][1, :]) * c
    total['lm_p12'] = (total['Pmom'][0, :] - total['Pmom'][1, :]) * c
    total['lm_p33'] = (total['Pmom'][2, :]) * c * 2
    total['lm_p34'] = -(total['Pmom'][3, :]) * c * 2

    del total['Pmom']


def integrated_mie(mie: sk.Mie, prob_dist: rv_continuous, refrac_index_fn, wavelengths,
                   num_quad=1024, maxintquantile=0.99999):
    """
    Integrates the Mie parameters over an arbitrary particle size distribution, returning scattering parameters such
    as cross sections, phase matrices, and legendre moments at a set of wavelengths.

    Note that the units of the input parameters are left arbitrary, but they must be consistent between each other,
    i.e., if wavelengths is specified in nm it is expected that the particle size distribution will be as a function
    of nm.

    Parameters
    ----------
    mie : sk.Mie
        Interal Mie algorithms to use
    prob_dist : rv_continuous
        Instance of a probability distribution to integrate over.  Distribution must be specified in the same units
        as wavelengths
    refrac_index_fn : fn
        Function taking in one argument (wavelength) and returning a complex number for refractive index at that
        wavelength.  Function argument should be in the same units as wavelengths.
    wavelengths : np.array
        Wavelengths to calculate for.  Units should be consistent with prob_dist
    num_quad : int, Optional
        Number of quadrature points to use when integrating over particle size.  Default is 1024
    maxintquantile : float
        Used to determine the maximum radius to use in integration.  A value of 0.99 means that at least 99% of the
        probability distribution area multipied by radius**2 will be included within the integration bounds.  Default
        0.99999
    Returns
    -------
    xr.Dataset
        Dataset containing the keys 'p11', 'p12', 'p33', 'p34', 'lm_p11', 'xs_total', 'xs_absorption', and 'xs_scattering'
        All scattering cross sections will be in units of wavelength**2
    """
    norm = integrate.quad(lambda x: prob_dist.pdf(x) * x**2, 0, 1e25, points=(prob_dist.mean()))[0]

    def pdf(x):
        return prob_dist.pdf(x) * x ** 2 / norm

    # We have to determine the maximum R to integrate to, this is apparently very challenging, what we do is
    # calculate the CDF of repeatadly large values until it is greater than our threshold
    max_r = prob_dist.mean()
    while (integrate.quad(pdf, 0, max_r * 2, points=(prob_dist.mean()))[0] - integrate.quad(pdf, 0, max_r, points=(prob_dist.mean()))[0]) > (1 - maxintquantile):
        max_r *= 2

    from scipy.special import roots_legendre

    x, w = roots_legendre(num_quad)
    x = 0.5 * (x + 1) * max_r
    w *= max_r / 2

    angles = mie.scattering_angles()
    max_moment = mie.max_legendre_moment()

    all_output = xr.Dataset({'p11': (['wavelength', 'angle'], np.zeros((len(wavelengths), len(angles)))),
                             'p12': (['wavelength', 'angle'], np.zeros((len(wavelengths), len(angles)))),
                             'p33': (['wavelength', 'angle'], np.zeros((len(wavelengths), len(angles)))),
                             'p34': (['wavelength', 'angle'], np.zeros((len(wavelengths), len(angles)))),
                             'lm_p11': (['wavelength', 'legendre'], np.zeros((len(wavelengths), max_moment))),
                             'lm_p12': (['wavelength', 'legendre'], np.zeros((len(wavelengths), max_moment))),
                             'lm_p33': (['wavelength', 'legendre'], np.zeros((len(wavelengths), max_moment))),
                             'lm_p34': (['wavelength', 'legendre'], np.zeros((len(wavelengths), max_moment))),
                             'xs_total': (['wavelength'], np.zeros(len(wavelengths))),
                             'xs_scattering': (['wavelength'], np.zeros(len(wavelengths))),
                             'xs_absorption': (['wavelength'], np.zeros(len(wavelengths)))
                             },
                            coords={'wavelength': wavelengths,
                                    'angle': angles})

    for idx, wavelength in enumerate(wavelengths):
        n = np.cdouble(refrac_index_fn(wavelength))

        total = dict()
        cext = []
        p11 = []
        mie.calculate(wavelength, x[0], n)
        start_params = _integrable_parameters(mie, prob_dist.pdf(x[0]))

        # Initialize the output
        for key, item in start_params.items():
            total[key] = np.zeros_like(item)

        # Integrate over every trapezoid
        for node, weight in zip(x, w):
            mie.calculate(wavelength, node, n)

            next_params = _integrable_parameters(mie, prob_dist.pdf(node))
            cext.append(next_params['Cext'])
            p11.append(next_params['p11'][0])

            for key in next_params:
                total[key] += weight * next_params[key]

        _post_process(total, wavelength)
        all_output['p11'].values[idx, :] = total['p11']
        all_output['p12'].values[idx, :] = total['p12']
        all_output['p33'].values[idx, :] = total['p33']
        all_output['p34'].values[idx, :] = total['p34']
        all_output['lm_p11'].values[idx, :] = total['lm_p11']
        all_output['lm_p12'].values[idx, :] = total['lm_p12']
        all_output['lm_p33'].values[idx, :] = total['lm_p33']
        all_output['lm_p34'].values[idx, :] = total['lm_p34']
        all_output['xs_total'].values[idx] = total['Cext']
        all_output['xs_scattering'].values[idx] = total['Csca']
        all_output['xs_absorption'].values[idx] = total['Cabs']
    return all_output
