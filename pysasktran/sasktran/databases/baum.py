import sasktran as sk
import xarray as xr
from scipy.special import legendre
from pathlib import Path
import numpy as np
from tqdm import tqdm


def generate_baum_legendre(numcoeff: int):
    config = sk.config.load_sasktran_registry()

    baum_dir = Path(config['software']['usask-arg']['skopticalproperties']['baum_icecrystals']['storage']['2014database'])

    baum_file = baum_dir.joinpath('GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc')

    if not baum_file.is_file():
        raise ValueError('Cannot find Baum database file')

    baum_data = xr.open_dataset(baum_file.as_posix())

    num_size = baum_data['nDeff'].size
    num_wavel = baum_data['nWaveLen'].size

    all_coeff = np.zeros((numcoeff, num_size, num_wavel))

    for i in range(0, num_size):
        print(i)
        for j in tqdm(range(0, num_wavel)):
            sliced = baum_data.isel(nDeff=i, nWaveLen=j)

            phase_vals = sliced['p11_phase_function'].values
            phase_angles = sliced['phase_angles'].values

            cos_angle = np.cos(np.deg2rad(phase_angles))
            # Reverse so cos_angle is ascending
            cos_angle = cos_angle[::-1]
            phase_vals = phase_vals[::-1]

            for k in range(0, numcoeff):
                lp = legendre(k)
                lp_vals = lp(cos_angle)

                all_coeff[k, i, j] = np.trapz(lp_vals * phase_vals, cos_angle) * (2 * k + 1) / 2.0

    # Normalize the coefficients so m=0 is 1
    all_coeff /= all_coeff[0, :, :][np.newaxis, :, :]

    baum_data['phase_moments'] = (['nPhaseMoment'], list(range(0, numcoeff)))
    baum_data['p11_phase_moments'] = (['nPhaseMoment', 'nDeff', 'nWaveLen'], all_coeff)

    out_file_name = baum_dir.joinpath(baum_file.stem + '_LegendreAdded.nc')
    baum_data.to_netcdf(out_file_name.as_posix())


if __name__ == "__main__":
    generate_baum_legendre(64)