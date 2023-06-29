import numpy as np
import logging
from io import BytesIO
from pathlib import Path
import appdirs
import re
import zipfile
import urllib3


def _package_data_folder() -> Path:
    return Path(__file__).resolve().parent.joinpath('data')


def _load_sorce(prefix: str, data_folder: str):
    files = list(Path(data_folder).glob(prefix + '*'))

    if len(files) == 0:
        return None, None, None
    else:
        fname = files[0]

        m = re.search('(\d+)_(\d+)$', fname.stem)

        return files[0], m.group(1), m.group(2)


def _download_sorce(url):
    http = urllib3.PoolManager()
    response = http.request('GET', url)

    zip = zipfile.ZipFile(BytesIO(response.data))

    return zip


def hitran_o2():
    file = _package_data_folder().joinpath('HITRANo2.txt')

    return file.as_posix()


def incoming_flux(date):
    data_dir = appdirs.user_data_dir('sasktran')
    if not Path(data_dir).exists():
        Path(data_dir).mkdir(parents=True)

    # Try to load in the FUV file
    fuv, start_date, end_date = _load_sorce('sorce_fuv', data_dir)

    if fuv is None or int(end_date) < date:
        for f in Path(data_dir).glob('sorce_fuv*'):
            f.unlink()
        # Download the fuv data from SORCE
        zip = _download_sorce(r'lasp.colorado.edu/data/sorce/ssi_data/solstice/fuv/txt/sorce_fuv_latest.txt.zip')
        zip.extractall(data_dir)
        fuv, start_date, end_date = _load_sorce('sorce_fuv', data_dir)

    if fuv is None:
        raise ValueError('Error downloading the FUV file from SORCE')

    # Try to load in the MUV file
    muv, start_date, end_date = _load_sorce('sorce_muv', data_dir)

    if muv is None or int(end_date) < date:
        for f in Path(data_dir).glob('sorce_muv*'):
            f.unlink()
        # Download the fuv data from SORCE
        zip = _download_sorce(r'http://lasp.colorado.edu/data/sorce/ssi_data/solstice/muv/txt/sorce_muv_latest.txt.zip')
        zip.extractall(data_dir)
        muv, start_date, end_date = _load_sorce('sorce_muv', data_dir)

    if muv is None:
        raise ValueError('Error downloading the MUV file from SORCE')

    # Try to load in the SIM file
    sim, start_date, end_date = _load_sorce('sorce_sim', data_dir)

    if sim is None or int(end_date) < date:
        for f in Path(data_dir).glob('sorce_sim*'):
            f.unlink()
        # Download the fuv data from SORCE
        zip = _download_sorce(r'http://lasp.colorado.edu/data/sorce/ssi_data/sim/txt/sorce_sim_latest.txt.zip')
        zip.extractall(data_dir)
        # TODO: truncate SIM data

        sim, start_date, end_date = _load_sorce('sorce_sim', data_dir)

    if sim is None:
        raise ValueError('Error downloading the SIM file from SORCE')

    try:
        import pandas as pd
        fuv_data = pd.read_csv(fuv.as_posix(), skiprows=85, header=None, delimiter='\s+').values
        muv_data = pd.read_csv(muv.as_posix(), skiprows=85, header=None, delimiter='\s+').values
        sim_data = pd.read_csv(sim.as_posix(), skiprows=85, header=None, delimiter='\s+').values
    except ImportError:
        logging.warning('Falling back to numpy to load in SORCE data.  Loading SORCE data can be sped up significantly by installing pandas')
        fuv_data = np.genfromtxt(fuv.as_posix(), skip_header=85)
        muv_data = np.genfromtxt(muv.as_posix(), skip_header=85)
        sim_data = np.genfromtxt(sim.as_posix(), skip_header=85)

    ix = np.where(fuv_data[:, 0] == date)
    if len(ix[0]) == 0:
        raise ValueError('Requested date not found in SORCE FUV data')
    fuv_data = fuv_data[ix[0], :]

    ix = np.where(muv_data[:, 0] == date)
    if len(ix[0]) == 0:
        raise ValueError('Requested date not found in SORCE MUV data')
    muv_data = muv_data[ix[0], :]

    ix = np.where(sim_data[:, 0] == date)
    if len(ix[0]) == 0:
        raise ValueError('Requested date not found in SORCE SIM data')
    sim_data = sim_data[ix[0], :]

    return fuv_data, muv_data, sim_data


def msis():
    file = _package_data_folder().joinpath('msis2.txt')

    return file.as_posix()


def xs_har():
    file = _package_data_folder().joinpath('datHar.txt')

    return file.as_posix()


def xs_lya():
    file = _package_data_folder().joinpath('datLya.txt')

    return file.as_posix()


def xs_src():
    file = _package_data_folder().joinpath('datSRC.txt')

    return file.as_posix()


if __name__ == "__main__":
    incoming_flux(5)