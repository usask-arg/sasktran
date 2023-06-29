.. _aband:

A Band Emission Model
*********************
SASKTRAN includes a model to include effects from photochemical emission in the oxygen A band.
This is included as :py:class:`ABandEmission <sasktran.ABandEmission>`, and can be created through ::

    import sasktran as sk
    import numpy as np

    emissions_alts_km = np.arange(0, 200)

    aband_emission = sk.ABandEmission(emissions_alts_km)

Running the Model
=================
To run the model we need three things

  1.  A representation of the atmospheric state, containing pressure, temperature, ozone, and O2.
  2.  A location to run the model at.
  3.  The solar geometry

The representation of the atmospheric state is an instance of the :py:class:`Atmosphere <sasktran.Atmosphere>` object, for example, ::

    atmosphere = sk.Atmosphere()

    # Set the background pressure/temperature
    atmosphere.atmospheric_state = sk.MSIS90()

    # Set the constituent species
    atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
    atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
    atmosphere['O2'] = sk.Species(sk.HITRANChemical('O2'), sk.MSIS90(include_o2=True), species='SKCLIMATOLOGY_O2_CM3')

Note that some of the climatologies used here are focused on the stratosphere, and may not be suitable for the mesosphere.
For more accurate calculations using the A band emission model it is recommended to create :py:class:`ClimatologyUserDefined <sasktran.ClimatologyUserDefined>` objects with realistic mesospheric profiles.

Having created the atmosphere we can then run the model, usually running the model takes around 40 seconds depending how many altitudes are requested. ::

    ref = [40, 0, 0, 54372]    # Lat/lon/altitude/mjd
    sun = np.array([0, 0, 1])

    emission_source, emission_wavel = aband_emission.run_model(atmosphere, ref, sun)

Including the Model in a Radiative Transfer Calculation
=======================================================
There are two ways we can use the emission model and include them in a radiative transfer calculation.
The first is to take the values we have already calculated and create a :py:class:`EmissionTable <sasktran.EmissionTable>` object. ::

    # Convert emission altitudes to m, and source function from /cm to /m
    emission_table = sk.EmissionTable(emissions_alts_km * 1000, emission_wavel, emission_source * 100)

The :py:class:`EmissionTable <sasktran.EmissionTable>` can then be added to the :py:class:`Atmosphere <sasktran.Atmosphere>` object that goes into the radiative transfer model.
We also have to tell the radiative transfer model to specifically include emissions in the calculation. ::

    atmosphere.emissions['aband'] = emission_table

    geo = sk.VerticalImage()
    geo.from_sza_saa(70, 60, 0, 0, np.array([30, 40, 50, 60, 70]), 54372, 0)

    engine = sk.EngineHR(geometry=geo, atmosphere=atmosphere, wavelengths=np.arange(760, 780, 1e-2))
    engine.include_emissions = True

    radiance = engine.calculate_radiance()

This method of including emissions in the model is a bit clunky since we have to first run the model at a specific location and solar geometry.

The second method uses the :py:class:`ABandEmission <sasktran.ABandEmission>` object directly instead of creating an intermediate :py:class:`EmissionTable <sasktran.EmissionTable>` ::

    atmosphere.emissions['aband'] = sk.ABandEmission(emissions_alts_km)

    geo = sk.VerticalImage()
    geo.from_sza_saa(70, 60, 0, 0, np.array([30, 40, 50, 60, 70]), 54372, 0)

    engine = sk.EngineHR(geometry=geo, atmosphere=atmosphere, wavelengths=np.arange(760, 780, 1e-2))
    engine.include_emissions = True

    radiance = engine.calculate_radiance()

The main advantage to this method is that the location and solar geometry used to run the A band emission model are taken directly from SASKTRAN to match the radiative calculation.
The downside is that the A band emission model is executed every time the radiative transfer calculation is performed.