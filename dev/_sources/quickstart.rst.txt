.. _quickstart:

Quickstart Guide
****************

Performing a calculation with ``sasktran`` can be broken up into three steps:

* Defining the geometry of the problem.  Typically this boils down to calculating an observer location,
  a look vector, and a timestamp for the measurements.  ``sasktran`` contains several tools to aid with this computation.

* Creating a representation of the atmospheric state for the calculation.

* Choosing a radiative transfer model and setting it up for the calculation.

Defining the Geometry
=====================
The raw input to ``sasktran`` for each measurement is the observer location, a unit look vector, and a timestamp
represented with a Modified Julian Date.  If you already have these things, excellent, we can directly create a
:class:`sasktran.Geometry` object, ::

    import sasktran as sk

    geometry = sk.Geometry()

    los_1 = sk.LineOfSight(observer=[3.676013154788849600e+005, 1.009976313640051500e+006, -6.871601202127538600e+006],
                           look_vector=[2.884568631765662100e-001, 7.925287180643269000e-001,  5.372996083468238900e-001],
                           mjd=54832.5)
    los_2 = sk.LineOfSight(observer=[3.692808540679614500e+005, 1.014590807988641800e+006, -6.870844156040793300e+006],
                           look_vector=[2.884568631765662100e-001, 7.925287180643269000e-001,  5.372996083468238900e-001],
                           mjd=54832.5)

    geometry.lines_of_sight = [los_1, los_2]

If you want to define the geometry based on solar angles, there are several convenience methods that can help::

    from sasktran.geometry import VerticalImage

    geometry = VerticalImage()
    geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=[10, 20, 30, 40], mjd=54372, locallook=0,
                          satalt_km=600, refalt_km=20)

Nadir viewing geometries can also be configured in a similar way::

    from sasktran.geometry import NadirGeometry

    geometry = NadirGeometry()
    geometry.from_zeniths_and_azimuths(solar_zenith=60, solar_azimuth=30, observer_mjds=[54372], observer_zeniths=[90], observer_azimuths=[60])


Representing the Atmosphere
===========================
Each atmospheric constituent is represented by two things:

* A climatology which specifies the amount and distribution of the constituent, the class :class:`sasktran.Climatology`
  helps here.  All climatologies are specified with respect to particle number density in units of /cm3
* An optical property of type :class:`sasktran.OpticalProperty` which defines the species cross sections.

The :class:`sasktran.Atmosphere` combines these things together for multiple species to make up the full atmospheric
state.  Many precomputed climatologies and optical properties are avaialable, or you can define your own.  But for now
let's create an atmosphere consisting of Rayleigh scattering and ozone absorption with precomputed climatologies::

    import sasktran as sk

    atmosphere = sk.Atmosphere()

    atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
    atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

Setting up the Engine
=====================
The final step is creating the :class:`sasktran.Engine` object which actually performs the radiative transfer
calculation.  `sasktran` contains multiple radiative transfer models, but for most standard applications it is
recommended to use the ``hr`` engine which is suitable for all viewing geometries in the UV-VIS-NIR spectral regime::

    import sasktran as sk
    from sasktran.geometry import VerticalImage

    # First recreate our geometry and atmosphere classes
    geometry = VerticalImage()
    geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=[10, 20, 30, 40], mjd=54372, locallook=0,
                          satalt_km=600, refalt_km=20)

    atmosphere = sk.Atmosphere()

    atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
    atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())

    # And now make the engine
    engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere)

    # Choose some wavelengths to do the calculation at
    engine.wavelengths = [340, 600]

    # And do the calculation
    radiance = engine.calculate_radiance()