import sasktran as sk
from . import load_scenario
import numpy as np


def test_all_multiscatter():
    for atm_scenario in range(1, 2):
        for geo_scenario in range(8, 9):
            scen = load_scenario(geo_scenario, atm_scenario, 0, 1, 16, altitude_spacing=1000)

            engine = sk.EngineCO(atmosphere=scen['atmo'], geometry=scen['geo'], wavelengths=scen['wavelengths'],
                                 options={'msmode': 1,
                                          'altitudegrid': scen['altitudes'],
                                          'applydeltascaling': False
                                          }
                                 )

            engine.nstokes = 3

            rad = engine.calculate_radiance('xarray')

            scia_vals = scen['model_data']['radiance'].isel(model=2).values[:, :, :3]
            p_diff = (rad['radiance'].values - scia_vals) / scia_vals * 100

            max_error = np.max(np.abs(p_diff.flatten()))
