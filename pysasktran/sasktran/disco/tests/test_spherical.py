import sasktran as sk

from sasktran.disco.tests.common import default_atmosphere, default_geometry, standard_test_wavelengths
import sasktran.disco.interface as do
import numpy as np


def test_nadir_viewing_same_radiance():
    atmo = default_atmosphere()
    geo = sk.NadirGeometry()

    geo.from_zeniths_and_azimuth_difference(60, 0, 0)

    wavelengths = standard_test_wavelengths()

    engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
    engine.layer_construction = 'match_altitude_grid'

    rad_plane_parallel = engine.calculate_radiance()

    engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
    engine.viewing_mode = 'spherical'
    engine.layer_construction = 'match_altitude_grid'

    rad_spherical = engine.calculate_radiance()

    p_diff = (rad_spherical['radiance'] - rad_plane_parallel['radiance']) / rad_spherical['radiance'] * 100

    np.testing.assert_almost_equal(p_diff.values.flatten(), 0.0, decimal=2)
