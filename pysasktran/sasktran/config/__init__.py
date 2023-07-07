import sys
import sysconfig
import os
import yaml
import appdirs
import urllib.request
import logging
import tarfile
from pathlib import Path
from packaging import version

# Location of the sasktranif registry yaml file
REGISTRY_FILE = os.path.join(sysconfig.get_path('data'), 'share', 'usask-arg', 'registry', 'sasktranif', 'globalkey.yaml')

# Mapping from option in the user config file to sasktranif registry setting
OPTION_MAPPINGS = {'hitran_directory': ['software', 'usask-arg', 'skopticalproperties',
                                        'hitran', 'basedirectory'],
                   'hitran_cache': ['software', 'usask-arg', 'skopticalproperties', 'hitran',
                                    'spectrallinecachedir'],
                   'baum_directory': ['software', 'usask-arg', 'skopticalproperties', 'baum_icecrystals', 'storage',
                                      '2014database'],
                   'ecmwf_directory': ['software', 'usask-arg', 'climatology', 'ecmwf', 'storage', 'erainterimdir'],
                   'pratmo_directory': ['software', 'usask-arg', 'climatology', 'pratmo', 'storage'],
                   'grib_samples_directory': ['software', 'usask-arg', 'grib', 'samples_folder'],
                   'grib_definitions_directory': ['software', 'usask-arg', 'grib', 'definitions_folder'],
                   'mie_aerosol_cache_directory': ['software', 'usask-arg', 'skopticalproperties', 'mieaerosol', 'aerosol_cache_directory']}

# Configuration management
APPDIRS = appdirs.AppDirs(appname='sasktran', appauthor='usask-arg')


def load_sasktran_registry():
    """
    Opens the sasktranif registry as a dictionary
    """
    with open(REGISTRY_FILE) as f:
        if version.parse(yaml.__version__) > version.parse('5'):
            config = yaml.load(f, Loader=yaml.FullLoader)
        else:
            config = yaml.load(f)
    return config


def user_config_file_location():
    """
    File location of the user config file
    """
    user_config_dir = APPDIRS.user_config_dir

    user_config_file = os.path.join(user_config_dir, 'config.yml')

    return user_config_file


def load_user_config():
    """
    Opens the user config file as a dictionary
    """
    user_config_file = user_config_file_location()

    try:

        with open(user_config_file) as f:
            if version.parse(yaml.__version__) > version.parse('5'):
                config = yaml.load(f, Loader=yaml.FullLoader)
            else:
                config = yaml.load(f)
        return config
    except FileNotFoundError:
        return dict()


def save_user_config(user_config: dict):
    """
    Saves the user config from a dictionry
    """
    user_config_file = user_config_file_location()

    Path(user_config_file).parent.mkdir(exist_ok=True, parents=True)

    with open(user_config_file, 'w') as f:
        yaml.dump(user_config, f, default_flow_style=False)


def merged_registry(user_config, sasktran_registry):
    """
    Merges the user config file into the sasktran registry dictionary.  Returns a dictionary of merged values.
    """
    for key, user_val in user_config.items():
        try:
            mapping = OPTION_MAPPINGS[key]
        except KeyError:
            # Either typo or an option that does not go into the registry
            continue

        # Traverse the registry except for the last entry
        base = sasktran_registry
        for m in mapping[:-1]:
            base = base[m]
        # Assign the last entry to the user value
        base[mapping[-1]] = user_val

    return sasktran_registry


def save_registry(registry):
    """
    Rewrites the sasktranif registry according to an input dictionary
    """
    old_registry = load_sasktran_registry()

    if old_registry != registry:
        with open(REGISTRY_FILE, 'w') as f:
            yaml.dump(registry, f, default_flow_style=False)


def update_registry_from_config():
    """
    Rewrites the sasktranif registry using the user config file to merge in new values.
    """
    try:
        user_config = load_user_config()
    except FileNotFoundError:
        # No config file, we dont have to do anything
        return
    sasktran_reg = load_sasktran_registry()
    merged = merged_registry(user_config, sasktran_reg)
    save_registry(merged)


def download_pratmo():
    """
    Downloads the PRATMO database file from the ARG web server and updates the registry to point to this new
    location
    """
    user_config = load_user_config()

    if 'pratmo_directory' in user_config:
        data_directory = Path(user_config['pratmo_directory'])
    else:
        data_directory = Path(APPDIRS.user_data_dir)

    if not data_directory.exists():
        data_directory.mkdir(parents=True)
    output_file = data_directory.joinpath('skClimatology_PratmoNO2.bin')

    if output_file.exists():
        logging.info('PRATMO file already exists, skipping download')

    else:
        try:
            urllib.request.urlretrieve('https://arg.usask.ca/sasktranfiles/skClimatology_PratmoNO2.bin',
                                       filename=output_file.as_posix())
        except Exception as e:
            logging.exception(e)

    user_config['pratmo_directory'] = data_directory.as_posix()
    save_user_config(user_config)

    update_registry_from_config()


def download_grib_def():
    """
    Downloads the Grib definitions files from the ARG web server and updates the registry path to point to the
    new values
    """
    user_config = load_user_config()

    if 'grib_samples_directory' in user_config and 'grib_definitions_directory' in user_config:
        grib_sample_folder = Path(user_config['grib_samples_directory'])
        grib_definitions_folder = Path(user_config['grib_definitions_directory'])
        grib_base_folder = grib_sample_folder.parent
    else:
        data_directory = Path(APPDIRS.user_data_dir)
        if not data_directory.exists():
            data_directory.mkdir(parents=True)

        grib_base_folder = data_directory.joinpath('grib')
        grib_sample_folder = grib_base_folder.joinpath('samples')
        grib_definitions_folder = grib_base_folder.joinpath('definitions')

    grib_sample_folder.mkdir(parents=True, exist_ok=True)
    grib_definitions_folder.mkdir(parents=True, exist_ok=True)

    tarfilename = grib_base_folder.joinpath('grib_definitions.tar.gz')

    if tarfilename.exists():
        tarfilename.unlink()

    try:
        urllib.request.urlretrieve('https://arg.usask.ca/sasktranfiles/grib_definitions.tar.gz',
                                   filename=tarfilename.as_posix())
    except Exception as e:
        logging.exception(e)

    tar = tarfile.open(tarfilename.as_posix())
    tar.extractall(grib_base_folder.as_posix())
    tar.close()

    if tarfilename.exists():
        tarfilename.unlink()

    user_config['grib_samples_directory'] = grib_sample_folder.as_posix()
    user_config['grib_definitions_directory'] = grib_definitions_folder.as_posix()
    save_user_config(user_config)

    update_registry_from_config()


def download_hitran():
    """
    Downloads the HITRAN files and updates the registry
    """
    from sasktran_core.hitran_manager import Hitran_HAPI_Manager

    user_config = load_user_config()

    data_directory = Path(user_config.get('hitran_directory', APPDIRS.user_data_dir + '/hitran/'))
    cache_directory = Path(user_config.get('hitran_cache', APPDIRS.user_data_dir + '/hitran/cache'))

    cache_directory.mkdir(exist_ok=True, parents=True)

    hitran = Hitran_HAPI_Manager(data_directory.as_posix())

    hitran.make_molparam_file()
    hitran.load_hitran_lines()

    user_config['hitran_directory'] = data_directory.as_posix()
    user_config['hitran_cache'] = cache_directory.as_posix()
    save_user_config(user_config)
    update_registry_from_config()
