import appdirs
from pathlib import Path
import shlex
import xarray as xr
import numpy as np
import sasktran as sk


def rtm_comparison_file():
    data_dir = Path(appdirs.user_data_dir('sasktran'))

    if not data_dir.exists():
        data_dir.mkdir(parents=True)

    file = data_dir.joinpath('zawada_AMT_rtm_comparison_data_v1.nc')

    if not file.exists():
        from zenodo_get import zenodo_get
        zenodo_get(shlex.split('--record 4292303 -o "{}"'.format(data_dir.as_posix())))

    return file


def rtm_atmosphere(anc_data, albedo_scen, atmo_scen):
    atmo = sk.Atmosphere()
    air_den = sk.ClimatologyUserDefined(anc_data.altitude.values,
                                        {'SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3': anc_data.air_numden.values})

    ozone_den = sk.ClimatologyUserDefined(anc_data.altitude.values,
                                          {'SKCLIMATOLOGY_O3NUMBERDENSITY_CM3': anc_data.ozone_numden.values})

    atmo['air'] = sk.Species(sk.SimpleRayleigh(), air_den)

    if atmo_scen > 0:

        opt_prop = sk.OpticalProperty('USERDEFINED_TABLES')
        opt_prop.skif_object().SetProperty('WavelengthTruncation', 1)

        opt_prop.skif_object().AddUserDefined(100, anc_data.wavelength.values, anc_data.ozone_absorption_cross_section.values)
        opt_prop.skif_object().AddUserDefined(400, anc_data.wavelength.values, anc_data.ozone_absorption_cross_section.values)

        atmo['ozone'] = sk.Species(opt_prop, ozone_den)

    if atmo_scen == 2:
        altitudes = anc_data.altitude.values
        values = anc_data.aerosol_numden.values
        species = sk.SpeciesAerosol(altitudes,
                                    {'SKCLIMATOLOGY_AEROSOL_CM3': values},
                                    {'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': np.ones_like(altitudes) * 0.08,
                                     'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': np.ones_like(altitudes) * 1.6}, 'H2SO4')

        f = '/Users/dannyz/Downloads/mie_sulfate_strat.nc'
        ds = xr.open_dataset(f)
        opt_prop = sk.UserDefinedScatterConstantHeight(ds.wavelength.values,
                                                       ds.xs_scattering.values,
                                                       ds.xs_absorption.values,
                                                       lm_a1=ds.lm_a1.values.T,
                                                       lm_a2=ds.lm_a2.values.T,
                                                       lm_a3=ds.lm_a3.values.T,
                                                       lm_a4=ds.lm_a4.values.T,
                                                       lm_b1=ds.lm_b1.values.T,
                                                       lm_b2=ds.lm_b2.values.T)

        species._optical_property = opt_prop

        atmo['aerosol'] = species

    if albedo_scen == 0:
        atmo.brdf = 0
    elif albedo_scen == 1:
        atmo.brdf = 0.3
    elif albedo_scen == 2:
        atmo.brdf = 1

    return atmo


def test_case_geometry(sza, saa):
    # First define the tangent points
    # Latitude 35.411411411411414 is very close to 6371 km Earth radius
    # Use longitude 0 always
    geometry = sk.Geometry()
    alts_m = np.arange(500, 80000, 1000)

    tangent_geo = sk.Geodetic()

    for a in alts_m:
        tangent_geo.from_lat_lon_alt(35.41133597835978, 0.0, float(a))

        spher_up = tangent_geo.local_up
        spher_south = tangent_geo.local_south
        spher_west = tangent_geo.local_west

        look = np.cos(np.deg2rad(saa)) * spher_south + np.sin(np.deg2rad(saa)) * spher_west

        observer = tangent_geo.location - (2000 * 1000) * look

        geometry.lines_of_sight.append(sk.LineOfSight(54372, observer, look))

    sun = spher_south * np.sin(np.deg2rad(sza)) + np.cos(np.deg2rad(sza)) * spher_up

    geometry.sun = sun
    # geometry.reference_point = [35.411411411411414, 0.0, 0.0, 54372]

    return geometry


def load_scenario(geometry_index: int, atmosphere_index: int, albedo_index: int, test_case: int, numlegendre: int,
                  altitude_spacing: float=500):
    rtm_file = rtm_comparison_file()

    geo = xr.open_dataset(rtm_file, group='geometry_data')
    anc = xr.open_dataset(rtm_file, group='ancillary_data')
    model = xr.open_dataset(rtm_file, group='model_data')

    geo = geo.isel(solar_condition=geometry_index)

    alts = np.arange(0, 100001, altitude_spacing)
    anc = anc.interp(altitude=alts)

    model = model.isel(solar=geometry_index, composition=atmosphere_index, albedo=albedo_index, test_case=test_case)

    result = dict()

    result['atmo'] = rtm_atmosphere(anc, albedo_index, atmosphere_index)
    result['geo'] = test_case_geometry(float(geo.tangent_sza), float(geo.tangent_saa))
    result['wavelengths'] = anc.wavelength.values
    result['model_data'] = model
    result['altitudes'] = alts

    return result
