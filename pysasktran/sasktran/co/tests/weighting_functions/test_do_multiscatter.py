import sasktran as sk
from sasktran.co.tests.common import default_atmosphere, default_geometry, standard_test_wavelengths
import numpy as np
import pytest
from .common import numerical_wf, validate_wf


NUM_WF_VERIFY = 60
NUMERICAL_EPS = 0.001


@pytest.mark.parametrize('delta_scale', [False, True], ids=['Delta Scale False', 'Delta Scale True'])
def test_absorption_wf(delta_scale):
    """
    Verifies that the weighting function for ozone absorption is equal to a numerical calculation up to machine precision
    In single scatter
    """
    ozone_wf_wavelengths = [310, 330, 600]

    atmo = default_atmosphere(altitude_spacing=1000)

    geo = default_geometry()

    engine = sk.EngineCO(atmosphere=atmo, geometry=geo, wavelengths=ozone_wf_wavelengths,
                         options={'msmode': 1,
                                  'numssmoments': 2,
                                  'numdostreams': 2,
                                  'applydeltascaling': delta_scale,
                                  'altitudegrid': atmo['o3'].climatology.altitudes
                                  })

    atmo.wf_species = ['o3']

    rad = numerical_wf(engine, atmo['o3'].climatology, 'o3', rel_eps=NUMERICAL_EPS*0.1)

    rad = rad.isel(perturbation=slice(0, NUM_WF_VERIFY))

    validate_wf(rad['wf_o3'], rad['wf_numerical_o3'], decimal=5)


@pytest.mark.parametrize('delta_scale', [True, False], ids=['Delta Scale True', 'Delta Scale False'])
def test_scattering_wf_one_scatterer(delta_scale):
    """
    Verifies that the weighting function for scattering is equal to a numerical calculation to 7 decimal places
    Only one scatterer is present so we are not checking dlegendre
    """
    ozone_wf_wavelengths = [310, 350, 600]

    atmo = default_atmosphere(altitude_spacing=1000)
    geo = default_geometry()

    engine = sk.EngineCO(atmosphere=atmo, geometry=geo, wavelengths=ozone_wf_wavelengths,
                         options={'msmode': 1,
                                  'numssmoments': 2,
                                  'numdostreams': 2,
                                  'applydeltascaling': delta_scale
                                  }
                         )

    atmo.wf_species = 'air'

    rad = numerical_wf(engine, atmo['air'].climatology, 'air', rel_eps=NUMERICAL_EPS)

    rad = rad.isel(perturbation=slice(0, NUM_WF_VERIFY))

    validate_wf(rad['wf_air'], rad['wf_numerical_air'])


@pytest.mark.parametrize('delta_scale', [True, False], ids=['Delta Scale True', 'Delta Scale False'])
def test_scattering_wf_multiple_scatterer(delta_scale):
    """
    Checks the weighting function for scattering when multiple scatterers are present.
    This checks dlegendre and the multiple weighting function capability of the model.
    """
    ozone_wf_wavelengths = [310, 350, 600]

    atmo = default_atmosphere(altitude_spacing=1000)

    atmo['air2'] = sk.Species(sk.SimpleRayleigh(), sk.MSIS90())

    geo = default_geometry()

    if delta_scale:
        ns = 2
    else:
        ns = 4

    engine = sk.EngineCO(atmosphere=atmo, geometry=geo, wavelengths=ozone_wf_wavelengths,
                         options={'msmode': 1,
                                  'numssmoments': ns,
                                  'numdostreams': ns,
                                  'applydeltascaling': delta_scale
                                  }
                         )

    atmo.wf_species = ['air2', 'air']

    rad = numerical_wf(engine, atmo['air'].climatology, 'air', rel_eps=NUMERICAL_EPS)

    rad = rad.isel(perturbation=slice(0, NUM_WF_VERIFY))

    validate_wf(rad['wf_air'], rad['wf_numerical_air'])


def test_wf_albedo():
    """

    """
    wavelengths = np.arange(280, 350, 1)

    atmo = default_atmosphere(altitude_spacing=1000)

    atmo['air2'] = sk.Species(sk.SimpleRayleigh(), sk.MSIS90())

    geo = default_geometry()

    engine = sk.EngineCO(atmosphere=atmo, geometry=geo, wavelengths=wavelengths,
                         options={'msmode': 1,
                                  'numssmoments': 4,
                                  'numdostreams': 4,
                                  'applydeltascaling': True
                                  }
                         )

    atmo.wf_species = ['brdf', 'air2']
    atmo.brdf = 0.3

    rad = engine.calculate_radiance('xarray')

    dbrdf = 0.001

    atmo.brdf = 0.3 + dbrdf

    rad_above = engine.calculate_radiance('xarray')

    atmo.brdf = 0.3 - dbrdf

    rad_below = engine.calculate_radiance('xarray')

    wf_brdf = (rad_above['radiance'] - rad_below['radiance']) / (2*dbrdf)

    pass