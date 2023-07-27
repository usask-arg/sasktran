import sasktran as sk

from sasktran.disco.tests.common import default_atmosphere, default_geometry, standard_test_wavelengths
import sasktran.disco.interface as do
import numpy as np
import xarray as xr


def test_coulsen_tables():
    # This recreates the Coulsen table tests that are done in the C++ code to make sure that the SasktranIF
    # layer is working correctly
    atmo = sk.Atmosphere()

    rayleigh = sk.SimpleRayleigh()

    od = 0.5
    w = 350

    xs = float(rayleigh.calculate_cross_sections(sk.MSIS90(), 0, 0, 10000, 54372, [w]).total)

    nd = od / (xs) / (100000 * 100)

    clim = sk.ClimatologyUserDefined([0, 100000], {'air': [nd, nd]})

    atmo['air'] = sk.Species(rayleigh, clim)
    atmo.brdf = 0

    # Set the zen=1 case to 0.9999999 instead because of azimuth ambiguity in the stokes vector
    # definition
    zen = np.rad2deg(np.arccos([0.02, 0.4, 0.9999999, 0.02, 0.4, 0.9999999]))
    azi = np.array([180, 180, 180, 120, 120, 120])

    sza = np.rad2deg(np.arccos(0.2))

    geo = sk.NadirGeometry()
    geo.from_zeniths_and_azimuth_difference(sza, zen, azi, observer_alt=1e8, solar_azimuth=180)

    correct_rad = np.array([0.44129802, -0.01753141, 0,
                            0.16889020, 0.01119511, 0,
                            0.05300496, 0.03755859, 0,
                            0.30091208, -0.15965601, 0.07365528,
                            0.12752450, -0.06066038, 0.05293867,
                            0.05300496, -0.01877930, 0.03252669])

    engine = do.EngineDO(geo, atmo, [w])

    engine.num_layers = 1
    engine.num_streams = 40
    engine.num_stokes = 3
    engine.options['usepseudospherical'] = 0
    engine.options['forcenumazimuthterms'] = 10

    rad = engine.calculate_radiance()
    rad['radiance'] *= np.pi

    np.testing.assert_almost_equal(correct_rad, rad['radiance'].values.flatten(), decimal=4)

    # Repeat with more layers, should be identical results
    engine = do.EngineDO(geo, atmo, [w])

    engine.num_layers = 20
    engine.num_streams = 40
    engine.num_stokes = 3
    engine.options['usepseudospherical'] = 0
    engine.options['forcenumazimuthterms'] = 10

    rad = engine.calculate_radiance()
    rad['radiance'] *= np.pi

    np.testing.assert_almost_equal(correct_rad, rad['radiance'].values.flatten(), decimal=4)


def test_inside_atmo():
    atmo = default_atmosphere()

    all_rad = []
    for obs_alt in [10, 20, 40, 50, 100, 120]:
        geo = sk.NadirGeometry()
        geo.from_zeniths_and_azimuths(60, 0, 54372, 0.001, 0, reference_point=(0, 0, 0, 54372),
                                      observer_altitudes=obs_alt*1000)

        engine = do.EngineDO(geometry=geo, atmosphere=atmo, wavelengths=[350])
        engine.num_stokes = 3

        all_rad.append(engine.calculate_radiance(output_format='xarray').isel(los=0, wavelength=0))

    all_rad = xr.concat(all_rad, dim='los')

    true_values = np.array([5.42336818e-02, 1.08188359e-02, -3.77652886e-06, 6.22770341e-02,
                            1.53645841e-02, -5.36331177e-06, 6.31686981e-02, 1.60992561e-02,
                            -5.61976340e-06, 6.32779067e-02, 1.61629650e-02, -5.64200223e-06,
                            6.42701851e-02, 1.65945662e-02, -5.79266110e-06, 6.42701851e-02,
                            1.65945662e-02, -5.79266106e-06])

    np.testing.assert_array_almost_equal(all_rad['radiance'].values.flatten(), true_values)
