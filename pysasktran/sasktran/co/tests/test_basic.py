import sasktran as sk
import numpy as np
import numpy as np
import time
import logging


class Timer(object):
    """
    Simple wrapper to time things easier.
    """

    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            logging.info(self.name)
        print('Elapsed: ' + str(time.time() - self.tstart) + 's')




tanalts_km = np.arange(10, 50, 0.2)

# First recreate our geometry and atmosphere classes
geometry = sk.VerticalImage()
geometry.from_sza_saa(sza=20, saa=20, lat=0, lon=0, tanalts_km=tanalts_km, mjd=54372, locallook=0,
                      satalt_km=600, refalt_km=20)

atmosphere = sk.Atmosphere()
atmosphere.brdf = 0

atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
#atmosphere['aerosol'] = sk.SpeciesAerosolGloSSAC()

atmosphere.wf_species = ['air', 'ozone']

# And now make the engine
engine = sk.EngineCO(geometry=geometry, atmosphere=atmosphere, options={'msmode': 0})

wavel = [340, 600]
# Choose some wavelengths to do the calculation at
engine.wavelengths = wavel

# And do the calculation
with Timer():
    radiance = engine.calculate_radiance('xarray')

alts = np.arange(0, 100001, 1000)
# And now make the engine
engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)
engine.wavelengths = wavel

engine.options['wfheights'] = alts
engine.options['wfwidths'] = np.gradient(alts)
engine.grid_spacing = 500

engine.num_orders_of_scatter = 1
#engine.prefill_solartransmission_table = True
#engine.use_solartable_for_singlescatter = True

# And do the calculation
with Timer():
    radianceHR = engine.calculate_radiance('xarray')

(radiance['radiance'] / radianceHR['radiance']).isel(wavelength=0).plot(y='los')
(radiance['radiance'] / radianceHR['radiance']).isel(wavelength=1).plot(y='los')

#import matplotlib.pyplot as plt
#plt.plot(radiance['wf_ozone'].isel(wavelength=1, los=30).values)
#plt.plot(radianceHR['wf_ozone'].isel(wavelength=1, los=30).values)