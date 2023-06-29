import pandas as pd
from pathlib import Path
from scipy import interpolate


def _from_osiris_file(path: Path):
    data = pd.read_csv(path.as_posix(), header=None)

    return interpolate.interp1d(data.values[:, 0], data.values[:, 1] - 1j * data.values[:, 2])


def refractive_index_fn_h2so4(source='osiris'):
    if source.lower() == 'osiris':
        data_file = Path(__file__).parent.joinpath('data/refrac_h2so4_osiris.txt')
        return _from_osiris_file(data_file)


def refractive_index_fn_ice(source='osiris'):
    if source.lower() == 'osiris':
        data_file = Path(__file__).parent.joinpath('data/refrac_ice_osiris.txt')
        return _from_osiris_file(data_file)


def refractive_index_fn_water(source='osiris'):
    if source.lower() == 'osiris':
        data_file = Path(__file__).parent.joinpath('data/refrac_water_osiris.txt')
        return _from_osiris_file(data_file)


def refractive_index_fn_dust(source='osiris'):
    if source.lower() == 'osiris':
        data_file = Path(__file__).parent.joinpath('data/refrac_dust_osiris.txt')
        return _from_osiris_file(data_file)
