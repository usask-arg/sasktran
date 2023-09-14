import sasktran as sk
from . import load_scenario
import numpy as np
import matplotlib.pyplot as plt


def test_all_multiscatter():
    for atm_scenario in range(2, 3):
        for geo_scenario in range(8, 9):
            scen = load_scenario(geo_scenario, atm_scenario, 1, 1, 16, altitude_spacing=1000)

            engine = sk.EngineCO(atmosphere=scen['atmo'], geometry=scen['geo'], wavelengths=scen['wavelengths'],
                                 options={'msmode': 2,
                                          'altitudegrid': scen['altitudes'],
                                          'numhriterations': 1,
                                          'numhrincoming': 350,
                                          'numhroutgoing': 350,
                                          'initializehrwithdo': 1,
                                          'applydeltascaling': True,
                                          'numthreads': 8
                                          }
                                 )

            engine.nstokes = 3

            rad = engine.calculate_radiance('xarray')

            scia_vals = scen['model_data']['radiance'].isel(model=2).values[:, :, :3]
            hr_vals = scen['model_data']['radiance'].isel(model=0).values[:, :, :3]

            mmm_vals = scen['model_data']['mmm'][:, :, :3]

            mmm_vals = scen['model_data']['radiance'].isel(model=1).values[:, :, :3]


            p_diff_I = ((rad['radiance'].values - mmm_vals) / mmm_vals * 100)[:, :, 0]
            p_diff_scia = ((scia_vals - mmm_vals) / mmm_vals * 100)[:, :, 0]
            p_diff_hr = ((hr_vals - mmm_vals) / mmm_vals * 100)[:, :, 0]


            rad_lp = np.sqrt((rad['radiance'].isel(stokes=1)**2 + rad['radiance'].isel(stokes=2)**2))
            scia_lp = np.sqrt(scen['model_data']['radiance'].isel(model=2).values[:, :, 1]**2 + scen['model_data']['radiance'].isel(model=2).values[:, :, 2]**2)

            p_diff_lp = ((rad_lp- scia_lp) / scia_lp * 100)

            max_error = np.max(np.abs(p_diff_I.flatten()))
