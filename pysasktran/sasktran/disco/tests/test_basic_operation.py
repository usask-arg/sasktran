from sasktran.disco.tests.common import default_atmosphere, default_geometry
import sasktran.disco.engine as do
import sasktran as sk
import xarray as xr
import numpy as np


def test_simple_calculation():
    atmo = default_atmosphere()
    geo = default_geometry()

    wavelengths = np.arange(280, 800, 1)

    engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)

    rad = engine.calculate_radiance(output_format='xarray')


def test_multiple_lines_of_sight():
    atmo = default_atmosphere()
    geo = default_geometry()

    wavelengths = np.arange(280, 800, 10)

    engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
    engine.options['averagereferencepoint'] = 0
    rad = engine.calculate_radiance(output_format='xarray')

    los_rads = []
    for los in geo.lines_of_sight:
        new_geo = sk.Geometry()
        new_geo.lines_of_sight = [los]
        new_geo.sun = geo.sun

        engine = do.EngineDO(geometry=new_geo, atmosphere=atmo, wavelengths=wavelengths)

        los_rads.append(engine.calculate_radiance(output_format='xarray').isel(los=0))

    los_rads = xr.concat(los_rads, dim='los')
    diff = rad - los_rads

    np.testing.assert_array_equal(diff['radiance'].values, 0)


def test_num_streams():
    atmo = default_atmosphere()
    geo = default_geometry(mjd=54732, zenith=0, azimuth=0)

    wavelengths = [360]
    all_rad = []
    for num_stream in np.arange(2, 41, 2):
        engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
        engine.num_streams = num_stream

        all_rad.append(engine.calculate_radiance(output_format='xarray').isel(los=0, wavelength=0))

    all_rad = xr.concat(all_rad, dim='streams')

    true_values = np.array([0.10616244, 0.11598831, 0.11630997, 0.11622187, 0.11619075,
                            0.11618959, 0.11619159, 0.11619204, 0.11619181, 0.1161916,
                            0.11619153, 0.11619153, 0.11619153, 0.11619153, 0.11619153,
                            0.11619152, 0.11619152, 0.11619152, 0.11619152, 0.11619152])

    np.testing.assert_array_almost_equal(all_rad['radiance'].values, true_values)


def test_num_layers():
    atmo = default_atmosphere()
    geo = default_geometry(mjd=54732, zenith=0, azimuth=0)

    wavelengths = [360]
    all_rad = []
    for num_layer in [1, 2, 10, 50, 100]:
        engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
        engine.num_layers = num_layer

        all_rad.append(engine.calculate_radiance(output_format='xarray').isel(los=0, wavelength=0))

    all_rad = xr.concat(all_rad, dim='streams')

    true_values = np.array([0.11622926, 0.11620927, 0.11619428, 0.11619204, 0.11619162])

    np.testing.assert_array_almost_equal(all_rad['radiance'].values, true_values)


def test_inside_atmo():
    atmo = default_atmosphere()

    all_rad = []
    for obs_alt in [10, 20, 40, 50, 100, 120]:
        geo = sk.NadirGeometry()
        geo.from_zeniths_and_azimuths(60, 0, 54372, 0, 0, reference_point=(0, 0, 0, 54372),
                                      observer_altitudes=obs_alt*1000)

        engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=[350])

        all_rad.append(engine.calculate_radiance(output_format='xarray').isel(los=0, wavelength=0))

    all_rad = xr.concat(all_rad, dim='los')

    true_values = np.array([0.05544467, 0.06379937, 0.06471301, 0.06482461, 0.0658386, 0.0658386])

    np.testing.assert_array_almost_equal(all_rad['radiance'].values, true_values)


def test_layer_construction():
    atmo = default_atmosphere()
    geo = default_geometry(mjd=54732, zenith=0.001, azimuth=0)

    # Check uniform pressure layers
    wavelengths = [360]
    all_rad = []
    for num_layer in [1, 2, 10, 50, 100]:
        engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
        engine.num_layers = num_layer
        engine.layer_construction = 'uniform_pressure'

        all_rad.append(engine.calculate_radiance(output_format='xarray').isel(los=0, wavelength=0))

    all_rad = xr.concat(all_rad, dim='streams')

    true_values = np.array([0.11622926, 0.11620927, 0.11619428, 0.11619204, 0.11619162])

    np.testing.assert_array_almost_equal(all_rad['radiance'].values, true_values)

    # Check uniform height layers
    all_rad = []
    for num_layer in [1, 2, 10, 50, 100]:
        engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
        engine.num_layers = num_layer
        engine.layer_construction = 'uniform_height'

        all_rad.append(engine.calculate_radiance(output_format='xarray').isel(los=0, wavelength=0))

    all_rad = xr.concat(all_rad, dim='streams')

    true_values = np.array([0.11622914, 0.1162065, 0.1161917, 0.11619113, 0.11619111])

    np.testing.assert_array_almost_equal(all_rad['radiance'].values, true_values)

    # Check matching altitude grid layers
    all_rad = []
    for alt_grid_spacing in [20000, 10000, 5000, 2500, 100]:
        engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
        engine.alt_grid = np.arange(0, 100001, alt_grid_spacing)
        engine.layer_construction = 'match_altitude_grid'

        all_rad.append(engine.calculate_radiance(output_format='xarray').isel(los=0, wavelength=0))

    all_rad = xr.concat(all_rad, dim='streams')

    true_values = np.array([0.12477913, 0.11846253, 0.11674142, 0.11632324, 0.1161911])

    np.testing.assert_array_almost_equal(all_rad['radiance'].values, true_values)

    # Check manual layer placements
    all_rad = []
    for layer_grid_spacing in [20000, 10000, 5000, 2500, 100]:
        engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=wavelengths)
        engine.layer_construction = np.arange(0, 100001, layer_grid_spacing)

        all_rad.append(engine.calculate_radiance(output_format='xarray').isel(los=0, wavelength=0))

    all_rad = xr.concat(all_rad, dim='streams')

    true_values = np.array([0.11619436, 0.11619192, 0.11619151, 0.11619137, 0.11619132])

    np.testing.assert_array_almost_equal(all_rad['radiance'].values, true_values)
