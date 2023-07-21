import sasktran as sk
import numpy as np


def default_atmosphere(altitude_spacing=500, ground_altitude=0, toa_altitude=100000,
                       lat=0, lon=0, mjd=54372) -> sk.Atmosphere:
    altitudes = np.arange(ground_altitude, toa_altitude + 0.1, altitude_spacing)

    atmo = sk.Atmosphere()

    o3_clim = sk.ClimatologyUserDefined(altitudes, {'o3': sk.Labow().get_parameter('SKCLIMATOLOGY_O3_CM3', lat, lon, altitudes, mjd)})
    air_clim = sk.ClimatologyUserDefined(altitudes, {'air': sk.MSIS90().get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', lat, lon, altitudes, mjd)})

    background_clim = sk.ClimatologyUserDefined(altitudes,
                                                {'SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3': sk.MSIS90().get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', lat, lon, altitudes, mjd),
                                                 'SKCLIMATOLOGY_PRESSURE_PA': sk.MSIS90().get_parameter('SKCLIMATOLOGY_PRESSURE_PA', lat, lon, altitudes, mjd),
                                                 'SKCLIMATOLOGY_TEMPERATURE_K': sk.MSIS90().get_parameter('SKCLIMATOLOGY_TEMPERATURE_K', lat, lon, altitudes, mjd)})

    atmo['o3'] = sk.Species(sk.O3DBM(), o3_clim)
    atmo['air'] = sk.Species(sk.Rayleigh(), air_clim)
    atmo.atmospheric_state = background_clim

    return atmo


def default_geometry(solar_zenith=20, solar_azimuth=0, mjd=None, zenith=None, azimuth=None,
                     lat=0, lon=0, obs_alt=600e3) -> sk.Geometry:
    geo = sk.NadirGeometry()

    if mjd is None:
        zenith = [0.0001, 20, 40, 40, 80, 80]
        azimuth = [0, 0, 0, 45, 0, 45]
        mjd = np.ones_like(zenith) * 54372

    geo.from_zeniths_and_azimuths(solar_zenith, solar_azimuth, mjd, zenith, azimuth, reference_point=(lat, lon, 0, mjd),
                                  observer_altitudes=obs_alt)

    return geo


def standard_test_wavelengths():
    return np.array([310, 330, 350, 525, 600, 675, 800])
