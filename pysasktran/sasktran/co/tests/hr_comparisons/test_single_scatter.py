import sasktran as sk
from sasktran.co.tests.common import default_atmosphere, default_geometry, standard_test_wavelengths
import numpy as np
import pytest

NUM_WF_VERIFY = 60
NUMERICAL_EPS = 0.0001


def test_singlescatter_scalar():
    """
    """
    wavelengths = standard_test_wavelengths()

    atmo = default_atmosphere(altitude_spacing=1000)

    geo = default_geometry()

    grid_spacing = 500
    engine_co = sk.EngineCO(atmosphere=atmo, geometry=geo, wavelengths=wavelengths,
                            options={'numssmoments': 16,
                                     'numdostreams': 16,
                                     'msmode': 0,
                                     'applydeltascaling': False,
                                     'altitudegrid': np.arange(0, 100001, grid_spacing)
                                     })

    engine_hr = sk.EngineHR(atmosphere=atmo, geometry=geo, wavelengths=wavelengths)
    engine_hr.grid_spacing = grid_spacing
    engine_hr.num_orders_of_scatter = 1

    rad_co = engine_co.calculate_radiance('xarray')
    rad_hr = engine_hr.calculate_radiance('xarray')

    p_diff = (rad_co - rad_hr) / rad_hr
