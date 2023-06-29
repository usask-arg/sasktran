.. _ISKEngine:

*********
ISKEngine
*********
   Create a radiative transfer engine. Several engines are available,
   
   * :ref:`enginehr`, *High Resolution* our latest successive orders, spherical radiative transfer model.  We recommend this model for most limb retrieval applications.
   * **MC**, *Monte Carlo*, a Monte-Carlo spherical radiative transfer model. Slow but accurate. Great for testing the accuracy of the successive orders models.
   * **SO**, *Successive Orders*, our first successive orders, spherical radiative transfer model. Still works and we still maintain the code but we have stopped developing new features.
   * **OCC**, *Occultation*, our first cut at an occultation engine. It is a faithfull reproduction of the optical depth model used in ACE-FTS data processing.  Optimized for use with HITRAN in the IR.
   
   Our goal has been to build a plug and play interface where the user can easily swap engines, atmospheres, cross-sections in their radiative transfer code without having to rewrite large amounts of code.
   In addition we have tried to keep the code fast so it can be used in operational scenarios. We think we have made great strides towards this goal but hope you can help us continue to improve the software.
   
   The basic radiative transfer calculation is a relatively simple process,
   
   #. Choose and create the engine of choice with :py:class:`ISKEngine`
   #. Create a description of the relevant atmospheric optical constituents with :meth:`~ISKEngine.AddSpecies` and :meth:`~ISKEngine.AddEmission`
   #. Describe the background atmosphere (usually **P** and **T**) with :meth:`~ISKEngine.SetAtmosphericState` 
   #. Describe the ground albedo with :meth:`~ISKEngine.SetAlbedo`
   #. Specify if polarized calculations are required with :meth:`~ISKEngine.SetPolarizationMode`
   #. Optional. Set and model specific options with :meth:`~ISKEngine.SetProperty`
   #. Define lines of sight by specifying the location of the observer and the look direction with :meth:`~ISKEngine.AddLineOfSight`
   #. Calculate the observed radiance with :meth:`~ISKEngine.CalculateRadiance` or :meth:`~ISKEngine.CalculateStokesVector`

   Different engines are used for different radiative transfer problems but they all use the same interface::
   
      import sasktranif.sasktranif as skif
      
      engine  = skif.ISKEngine('HR')


.. py:class:: ISKEngine(name)

   Constructor for the class.    

   :param str name:
      The name of the engine to be created. The name must correspond to 
      an installed engine,
      
      +---------+---------------------------------------------------------------------------------+
      | Value   | Description                                                                     |
      +=========+=================================================================================+
      |  \'HR\' | The high resolution, successive orders model. A fast, highly configurable model.|
      +---------+---------------------------------------------------------------------------------+
      |  \'MC\' | The Monte Carlo spherical radiative transfer engine                             |
      +---------+---------------------------------------------------------------------------------+
      |  \'SO\' | Old successive orders model. Still available but no longer upgraded             |
      +---------+---------------------------------------------------------------------------------+
      |  \'OCC\'| The occultation engine. A first cut at an occultation engine.                   |
      +---------+---------------------------------------------------------------------------------+


AddEmission
^^^^^^^^^^^
.. py:method:: ISKEngine.AddEmission(species, emission) -> ok

      Adds an emission object to the engine’s internal table of atmospheric optical
      properties. The emission objects are used to model light internally generated
      by the atmosphere, for example, photochemical, auroral and thermal emissions::

         ok = engine.AddEmission( species, emission)

      :param CLIMATOLOGY_HANDLE speciesid:
         The unique identification code of the emission object being added to the engine.
         The engine will use this identification to either add this emission object to its internal table or replace an existing entry.

      :param ISKEmission emission:
         The emission object to be added to the engines internal table. Note that the engine
         will divide the emission objects isotropic radiance by the irradiance of the top-of-the atmosphere sun

      :param boolean ok:
         The returned value. Evaluates to true if the call succeeds, false if it does not.

      :return:
         Returns true if successful

AddLineOfSight
^^^^^^^^^^^^^^
.. py:method:: ISKEngine.AddLineOfSight(mjd, observer, lookvector) -> ok,losindex

      Adds an observer’s location, time, and look vector to the list of rays stored in the engine.
      Radiance will be calculated along this ray, as well as all of the others previously defined,
      upon the next call to :meth:`~ISKEngine.CalculateRadiance` and :meth:`~ISKEngine.CalculateStokesVector`::

          ok,losindex = engine.AddLineOfSight( mjd, observer,lookvector);

      :param double mjd:
         The time of the line-of-sight measurement. Most of the engines implicitly assume
         that all of the rays associated with a given calculation are approximately at the same time.
         The time of the measurement typically affects the position of the sun chosen to represent
         all of the rays and the time used for the number densities, pressures and temperatures
         extracted from climatologies. As a rule of thumb we note that the Sun, which has an
         angular diameter of 0.5°, moves 0.5° across the sky in 2 minutes, thus
         all measurements should occur within a 2 minute window for most of the engines.

      :param nxVector observer:
         The location of the observer expressed as a 3 element array of X, Y, Z geographic distances from the centre of the
         Earth location in meters. The X and Y coordinates are in the equatorial plane and X is aligned 
         with the Greenwich Meridian at 0° longitude. Y points to 90° longitude and Z is parallel to the
         spin axis of the Earth and points northward.

      :param nxVector lookvector:
         The unit look vector of the ray away from observer expressed as a 3 element array of X, Y, Z.
         The geographic X, Y, Z directions are the same as for observer. Note that the engine makes
         the light flow towards the observer even though the look-vector is away from the observer.

      :param boolean ok:
         The first element of the returned list. Evaluates to true if the call succeeds, false if it does not.

      :param integer losindex:
         The second element of the returned list. Returns the numeric index of the line of sight as stored within the engine. This can be used to index
         the 2-D radiance array returned by :meth:`~ISKEngine.CalculateRadiance` and :meth:`~ISKEngine.CalculateStokesVector`         

      :return: A two element list [*ok,losindex*]

AddSpecies
^^^^^^^^^^
.. py:method:: ISKEngine.AddSpecies(species,climatology,opticalproperty) -> ok

      Adds a species to the engine’s internal table of atmospheric optical properties. This table is used
      to fully describe the optical properties of the atmosphere.  The engine uses the table as it executes radiative
      transfer calculations.
      
      The method will add a new species to the internal table if an entry does not exist for that species but 
      replaces existing entries that do exist. The new entries are used in the next call to CalculateRadiance. The method
      allows the caller to send an empty instance of ISKOpticalProperty in which case an existing entry must exist and 
      only the climatology is replaced.  This is convenient for updating species profiles from within retrieval code.
      
      The calculation of atmospheric absorption, scattering or extinction by a specific atmospheric species requires knowledge about the
      cross-sections of the individual atoms/molecules/particles and knowledge about the
      number density of the atoms/molecules/particles. The opticalproperty parameter provides
      knowledge about the individual cross-sections while the climatology parameter provides
      knowledge about the number density::

         ok = engine.AddSpecies( species, climatology, opticalproperty)        # Add new or replace existing entry. Updates both climatology and optical properties.
         ok = engine.AddSpecies( species2, climatology, ISKOpticalProperty() ) # replace only the climatology of existing entry for species2

      :param string species:
         The identification code of the species being defined. The engine will add this species
         to its internal optical properties table or it will overwrite an existing entry with the
         same identification code. The identification code must be supported by the climatology
         parameter and must return the species number density

      :param ISKClimatology climatology:
         The climatology of the species being considered. This climatology should return the number
         density (cm-3) of the atoms/molecules/particles at a specific point in time and space when 
         passed the species parameter.  The cross-section returned by the opticalproperty parameter
         will be multiplied by this number density.

      :param ISKOpticalProperty opticalproperty:
         The object used to calculate the cross-sections and optical properties of individual
         atoms/molecules/particles of the species being considered. Many atoms/molecules/particles
         have cross-sections that depend upon atmospheric state, e.g. pressure and temperature and
         this is controlled by a call to SetAtmosphericState.

      :param boolean ok:
         The returned value. Evaluates to true if the call succeeds, false if it does not.

      :return: returns true if successful

GetWeightingFunctions
^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKEngine.GetWeightingFunctions()-> (ok,wf)

      Returns the 3-D weighting functions for the current wavelength, lines of sight and volume. The :ref:`enginehr`
      engine is the only engine that currently supports the feature::

         ok,wf = engine.GetWeightingFunctions()

      :param array wf:
         The ``wf`` object s returned as a three dimensional ``numpy.ndarray`` with dimensions corresponding to
         ``[wavelength, line of sight, volume]``,

         .. math::
            \texttt{wf[i, j, k]} = \frac{\partial I(\lambda_i, \text{LOS}_j)}{\partial x_k},

         and has units of ``[radiance/cm^{-3}]``.  In our example, the quantity :math:`x_k` is the ozone number density
         over a finite volume.  Since we set ``calcwf = 2`` the finite volumes are uniform spherical shells spaced ``1 km`` apart evenly
         from ``0.5 km`` to ``99.5 km``.  Therefore, the quantity ``wf[i, j, 10]`` can be thought of as the derivative of
         the radiance (for wavelength ``i`` and line of sight ``j``) with respect to changing ozone number density in the     
         ``10.5 km`` shell.

      :return: status ok and parameter wf are returned as a tuple.

SetAlbedo
^^^^^^^^^
.. py:method:: ISKEngine.SetAlbedo( albedo) -> ok

      Sets the albedo of the ground or minimum height considered in the model. The albedo
      is implemented as the ratio of upward flux to downward flux and is normally 
      assumed to be Lambertian. The current interface uses the same albedo for all wavelengths::

         ok = engine.SetAlbedo ( albedo );

      :param double albedo:
         The ground albedo used in subsequent radiative transfer calculations. The value
         will be typically a number between 0 and 1. Most engines can handle albedos greater
         than 1 but may not properly handle negative albedos as many engines truncate negative
         intermediate radiances to zero.  

      :param bool ok:
         Returns true if successful

      :return: True if successful

SetAtmosphericState
^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKEngine.SetAtmosphericState(climatology) -> ok

      Sets the climatology used to calculate background atmospheric state at all locations
      in the atmosphere. The atmospheric state is used by the optical properties
      in the engine (see `~ISKEngine.AddSpecies`) to calculate as-needed atmospheric state parameters such
      as pressure and temperature. The atmospheric state climatology will typically support
      pressure and temperature as a minimum although it may have to support other parameters 
      depending upon the needs of the individual optical properties. The user is responsible
      for ensuring the atmospheric state used by the engine is appropriate for the 
      optical properties used::

         ok = engine.SetAtmosphericState( climatology )


      :param ISKClimatology climatology:
         The climatology that will be used for background atmospheric state in subsequent calls
         to `~ISKEngine.CalculateRadiance` and `~ISKEngine.CalculateStokesVector`

      :param bool ok:
         Returns true if successful

      :return: True if successful

SetPolarizationMode
^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKEngine.SetPolarizationMode( mode ) -> ok

      Sets the polarization mode used by the engines.The HR and MC engines can perform polarized or scalar 
      calculations. The SO and OCC engines only support scalar calculations.

      Note that if polarized calculations are requested they will be performed even if the 
      user requests scalar radiances by calling :meth:`~ISKEngine.CalculateRadiance`. Similarly, if scalar calculations
      are requested they will be performed even if the user requests full Stokes vectors by calling 
      :meth:`~ISKEngine.CalculateStokesVectors`::

         ok = engine.SetPolarizationMode(mode) 

      :param integer mode:
         Set the polarization mode used in subsequent calls to :meth:`~ISKEngine.CalculateRadiance` or :meth:`~ISKEngine.CalculateStokesVector`.
         There are several values for the polarization mode used inside various engines which are outlined
         in the table below, 

         +-------+--------------+-----------+------------------------------------------------+
         | Value | Polarization | Engine    |      Description                               |
         |       |    Type      | Support   |                                                |
         +=======+==============+===========+================================================+
         |  0    | none         | All       | Scalar calculations. Default value             |
         +-------+--------------+-----------+------------------------------------------------+
         |  1    | pseudopol1   | MC, HR    | Polarization to 1st order scatter              |
         +-------+--------------+-----------+------------------------------------------------+
         |  2    | pseudopol2   | MC, HR    | Polarization to 2nd order scatter              |
         +-------+--------------+-----------+------------------------------------------------+
         |  3    | pseudopol3   | MC, HR    | Polarization to 3rd order scatter              |
         +-------+--------------+-----------+------------------------------------------------+
         |  4    | pseudopol3c  | MC, HR    | Polarization to 3rd order scatter. C method.   |
         +-------+--------------+-----------+------------------------------------------------+
         |  5    | pseudopol4c  | MC, HR    | Polarization to 4th order scatter. C method    |
         +-------+--------------+-----------+------------------------------------------------+
         |  6    | pseudopol4d  | MC, HR    | Polarization to 4th order scatter. D method    |
         +-------+--------------+-----------+------------------------------------------------+
         |  99   | fullpol      | MC        | Full polarization. NOT YET IMPLEMENTED         |
         +-------+--------------+-----------+------------------------------------------------+

SetWavelengths
^^^^^^^^^^^^^^
.. py:method:: ISKEngine.SetWavelengths( wavelengths)

      Sets the wavelengths for subsequent radiative transfer calculations. The next call to
      :meth:`~ISKEngine.CalculateRadiance or :meth:`~ISKEngine.CalculateStokesVector` will calculate 
      radiance along each line of sight for each of the wavelengths. The wavelengths replace any
      wavelengths defined in earlier calls::

         ok = engine.SetWavelengths( wavelengths )

      :param array wavelenarray:
         An array of wavelengths expressed in nanometers. Radiance will be calculated at all of these wavelengths for each line of sight.

      :param bool ok:
         Returns true if successful

      :return: True if successful

CalculateRadiance
^^^^^^^^^^^^^^^^^
.. py:method:: ISKEngine.CalculateRadiance()->[ok, radiance]

      Calculates the scalar radiance for each ray defined by each previous call to :meth:`~ISKEngine.AddLineOfSight` and for each wavelength
      defined by the last call to :meth:`~ISKEngineSetWavelengths`. Note that the engine may use a polarized
      radiative transfer model even though a scalar radiance depending upon the value of :meth:`~ISKEngine.SetPolarizationMode`.
      The radiance returned assumes that the top-of-the-atmosphere solar irradiance is 1.0 at all wavelengths. 
      The engine returns the radiance to the user as a 2-D matrix of doubles:: 

         ok,radiance = engine.CalculateRadiance()

      :param boolean ok:
         The first element of the returned list. Evaluates to true if the call succeeds, false if it does not.

      :param matrix radiance:
         The second element of the returned list. Returns the scalar radiance as 2-D matrix of doubles (wavelengths, lines-of-sight).  Note that the occultation
         engine returns optical depth rather than radiance.

      :return: A two element list [*ok,radiance*]

CalculateStokesVector
^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKEngine.CalculateStokesVector()->[ok, stokes]

      Calculates the vector radiance for each ray defined by each previous call to AddLineOfSight and for each wavelength
      defined by the last call to SetWavelengths. Note that the engine may (or may not) use a scalar 
      model internally even though a vector radiance is returned. The calculation may take substantial time.
      The radiance is calculated assuming that the top-of-the-atmosphere solar irradiance is 1.0 at all wavelengths. 
      The engine returns the radiance to the user as a 2-D matrix of ISKStokesVector objects:: 

         ok,stokes = engine.CalculateStokesVector()

      :param boolean ok:
         The first element of the returned list. Evaluates to true if the call succeeds, false if it does not.

      :param matrix stokes:
         The second element of the returned list. Returns the scalar radiance as 2-D matrix of ISKStokesVector (wavelengths, lines-of-sight). 

      :return: A two element list [*ok,stokes*]


IsValidObject
^^^^^^^^^^^^^
.. py:method:: ISKEngine.IsValidObject()

      Checks to see if the underlying C++ engine is properly created. The function is primarily intended for internal usage::

         ok = engine.IsValidObject()

      :param boolean ok:
         Returns true if the call succeeds, false if it does not.

      :return: returns true if successful

SetProperty
^^^^^^^^^^^
.. py:method:: ISKEngine.SetProperty(args)

      Sets a property

GetProperty
^^^^^^^^^^^
.. py:method:: ISKEngine.GetProperty(args)

      Fetches a property



