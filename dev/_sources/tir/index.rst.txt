.. _tir_main:

SASKTRAN Thermal InfraRed Engine
================================

The TIR engine solves the radiative transfer equation in the thermal infra-region on a true spherical Earth.

.. important::

  Unlike other SASKTRAN engines, the TIR engine does not produce radiances with sun-normalized units. Solar
  radiation is negligible in the thermal infrared regime and using sun-normalized units would require
  the thermal emissions to be divided by a solar spectrum, which would add uncertainty to the calculation
  result.

  The TIR engine returns radiances with units of :math:`\mathrm{photon}/(\mathrm{s \: cm^2 \: nm \: sr})`.

  Wavelengths must be passed to the engine in increasing order.

The following example shows how to calculate radiance using the TIR engine. As with other SASKTRAN engines, there
are three inputs required by the TIR engine object:

* An atmospheric object which defines the species contained in the model atmosphere. Each species is defined by
  a climatology and optical property.
* A geometry object which defines the lines of sight along for which the radiance will be computed.
* An array of wavelengths to calculate the radiance at.

::

    from sasktran.tir.engine import EngineTIR
    from sasktran.tir.opticalproperty import HITRANChemicalTIR
    import sasktran as sk
    import numpy as np

    # select wavelengths
    wavelen = np.arange(13077.02, 13097.58, 0.005)

    # specify a limb line of sight
    geometry = sk.VerticalImage()
    geometry.from_sza_saa(sza=0, saa=0, lat=45, lon=0, tanalts_km=[20], mjd=56300, locallook=0.0)

    # create an atmosphere containing ozone
    atmosphere = sk.Atmosphere()
    atmosphere.atmospheric_state = sk.MSIS90()
    atmosphere['O3'] = sk.Species(HITRANChemicalTIR('O3'), sk.Labow())

    # create the engine
    engine = EngineTIR(geometry=geometry, atmosphere=atmosphere, wavelengths=wavelen)

    # do the calculation
    radiance = engine.calculate_radiance()

See :ref:`tir_examples` for more examples of usage.


Contents:

.. toctree::
   :maxdepth: 2

   sasktran_tir
   weighting_functions
   hitran
   climatology
   properties
