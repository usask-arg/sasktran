.. _properties:

******************
Full Property List
******************

When using the ``sasktran`` package properties are set through::

    engine.options['PropertyName'] = value

or with::

    ok = isk_engine.SetProperty('PropertyName', value)

when using the ``sasktranif`` package.
The property name is not case sensitive
Depending on the property, value may be a scalar, list, or array.

A complete list of properties is provided for reference.


Alphabetical List of Properties 
--------------------------------

===================================== =======================
Property                              Comments
===================================== =======================
:ref:`tir_addwfspecies`               -
:ref:`tir_calctemperaturewf`          -
:ref:`tir_geoidmodel`                 -
:ref:`tir_groundemissivity`           -
:ref:`tir_look`                       -
:ref:`tir_minextinctionratioofcell`   -
:ref:`tir_maxopticaldepthofcell`      -
:ref:`tir_numlinesofsight`            -
:ref:`tir_numthreads`                 -
:ref:`tir_observer`                   -
:ref:`tir_opticalpropertiesheightres` -
:ref:`tir_opticaltabledimensions`     -
:ref:`tir_raytracingshells`           -
:ref:`tir_referencepoint`             -
:ref:`tir_sourcetermorder`            -
:ref:`tir_surfaceheight`              -
:ref:`tir_tangentalts`                -
:ref:`tir_toaheight`                  -
:ref:`tir_useadaptiveintegration`     -
:ref:`tir_usecache`                   -
:ref:`tir_uselinearextinction`        -
:ref:`tir_userefraction`              -
:ref:`tir_usevmrwfunit`               -
:ref:`tir_wavel`                      -
:ref:`tir_wavelengths`                -
:ref:`tir_wfheights`                  -
:ref:`tir_wfspecies`                  -
:ref:`tir_wfwidths`                   -
===================================== =======================


Weighting Functions
-------------------

General weighting function description goes here.

..  _tir_addwfspecies:

AddWFSpecies
^^^^^^^^^^^^

.. option:: AddWFSpecies (string s) [Default: None]

    Tells the engine to compute weighting functions for a species. The input is a string with the ``CLIMATOLOGY_HANDLE``
    of the desired species. E.g. to compute CO2 weighting functions, ``s="SKCLIMATOLOGY_CO2_CM3"``.

.. _tir_wfheights:

WFHeights
^^^^^^^^^

.. option:: WFHeights (array v) [Default: linspace(500, 99500, 100)]

    Heights, in meters, to compute analytic weighting functions at.

.. _tir_wfspecies:

WFSpecies
^^^^^^^^^

.. option:: WFSpecies (string s) [Default: None]

    Tells the engine which species to compute analytic weighting functions for. The input string is a list of comma- or
    space-separated climatology handles. E.g. to compute weighting functions for CO2 and O3,
    ``s="SKCLIMATOLOGY_CO2_CM3, SKCLIMATOLOGY_O3_CM3"``. Setting this option overwrites any previous weighting function
    species setting.

.. _tir_wfwidths:

WFWidths
^^^^^^^^

.. option:: WFWidths (array v) [Default: ones(100) * 1000]

    Weighting function widths. The weighting functions determine the effect of perturbing a species population at
    each height specified by :option:`WFHeights` on the radiance. The property :option:`WFWidths` species the vertical distance
    from each value in :option:`WFHeights` to the height where the perturbation is zero. Therefore, :option:`WFHeights` and :option:`WFWidths`
    must have identical lengths.

    i.e. for a weighting function calculation at ``wfheights[i]``, the perturbation is maximum at ``wfheights[i]`` and
    decreases linearly (with height) to 0 at ``(wfheights[i] - wfwidths[i])`` and ``(wfheights[i] + wfwidths[i])``.

.. _tir_usevmrwfunit:

UseVMRWFUnit
^^^^^^^^^^^^

.. option:: UseVMRWFUnit (int n) [Default: n=0]

    Tells the TIR engine to return weighting functions with VMR units. This means the weighting functions will have
    the form :math:`\mathrm{d}I / \mathrm{d(VMR)}`. By default this option is not set and the weighting functions are
    computed as the change in radiance due to a change in the species number density (with units :math:`\mathrm{cm^{-3}}`).

.. _tir_calctemperaturewf:

CalcTemperatureWF
^^^^^^^^^^^^^^^^^

.. option:: CalcTemperatureWF (int n) [Default: n=0]

    A boolean property, that if set to true (n=1), instructs the TIR engine to calculate weighting functions for
    temperature. If weighting functions of gas number density or VMR were also computed, the temperature weighting
    functions are appended to the gas species weighting functions.


Integration
-----------

The integration options affect the performance and accuracy of the line of sight integration.

.. _tir_minextinctionratioofcell:

MinExtinctionRatioOfCell
^^^^^^^^^^^^^^^^^^^^^^^^

.. option:: MinExtinctionRatioOfCell (double d) [Default: d=0.9]

    Set the minimum ratio between the extinction at the endpoints of single cell that is allowed before the cell
    is split, when adaptive integration is enabled. If the option :option:`UseAdaptiveIntegration` is set to False, this
    option has no effect. The ratio is calculated as the smaller of the two extinction values divided by the
    larger.

.. _tir_maxopticaldepthofcell:

MaxOpticalDepthOfCell
^^^^^^^^^^^^^^^^^^^^^

.. option:: MaxOpticalDepthOfCell (double d) [Default: d=0.1]

    Set the maximum optical depth of a single cell is allowed to have before it is split into two cells, when
    adaptive integration is enabled. If the option :option:`UseAdaptiveIntegration` is set to False, this option has no
    effect.

.. _tir_sourcetermorder:

SourceTermOrder
^^^^^^^^^^^^^^^

.. option:: SourceTermOrder (int n) [Default: n=0]

    Determines how the source term is represented within each integration cell. Allowed values are 0 or 2.
    The default 0th order option makes the source term constant across the cell. If a value of 2 is used,
    the source function is quadratically splined along the line of sight.

    NOTE: The analytic weighting functions calculated by TIR assume a constant source term in each cell,
    regardless of this setting. Using a 2nd order source term will only effect the radiance returned by
    the model.

.. _tir_useadaptiveintegration:

UseAdaptiveIntegration
^^^^^^^^^^^^^^^^^^^^^^

.. option:: UseAdaptiveIntegration (int n) [Default: n=1]

    A boolean property that, if True, enables dynamic cell splitting during the radiative transfer integration.
    During the ray tracing process, the line of sight is split into cells based on the intersection of the line of
    sight with geometric shapes used to split the atmosphere into homogeneous layers. When the optical depth of a
    cell is computed, the cell is split if one or both of the following conditions are satisfied:

      1. The optical depth of the cell exceeds the value specified by the property :option:`MaxOpticalDepthOfCell`
      2. The ratio of the smaller of the two values of extinction at the endpoints of the cell to the larger
         value is less than the property :option:`MinExtinctionRatioOfCell`

.. _tir_uselinearextinction:

UseLinearExtinction
^^^^^^^^^^^^^^^^^^^

.. option:: UseLinearExtinction (int n) [Default: n=1]

    If this property is set to True, the engine calculates the optical depth of a cell by allowing the extinction
    to vary linearly with height, between the values at the endpoints of the cell. If this property is set to False,
    the extinction of a cell is constant, calculated as the average of the values at the endpoints.


Geometry and Atmosphere
-----------------------

..  _tir_geoidmodel:

GeoidModel
^^^^^^^^^^

.. option:: GeoidModel (string s) [Default: s="WGS84"]

    Manually sets the reference ellipsoid used in the model coordinate system. Allowed values are
    ``"GEOID_SPHERE"``, ``"IAU1976"``, ``"GRS80"``, ``"MERIT83"``, and ``"WGS84"``. Mainly useful for comparing with
    other radiative transfer models.

.. _tir_groundemissivity:

GroundEmissivity
^^^^^^^^^^^^^^^^

.. option:: GroundEmissivity (double d) [Default: d=1.0]

    Sets the scalar emissivity of the surface of the Earth. If the line of sight ends at the ground, the emission
    from the Earth's surface must be taken into account. In the TIR engine this emission is modelled as a
    black body with constant emissivity. For limb-viewing geometries this setting will have no effect as there are
    no scattering effects in the TIR engine.

.. _tir_opticalpropertiesheightres:

OpticalPropertiesHeightRes
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. option:: OpticalPropertiesHeightRes (double d) [Default: d=1000.0]

    Sets the height resolution (in meters) of the optical properties table. The optical properties table stores the
    values of species cross sections, absorption, and thermal emissions at discrete heights. To obtain the value
    of one of these properties at any height, the model performs an interpolation. The table values are distributed
    between :option:`SurfaceHeight` and :option:`TOAHeight` with a spacing of :option:`OpticalPropertiesHeightRes`.
    Decreasing this spacing will improve model accuracy as there will be less interpolation error, but will result
    in increased computation times.

.. _tir_opticaltabledimensions:

OpticalTableDimensions
^^^^^^^^^^^^^^^^^^^^^^

.. option:: OpticalTableDimensions (int n) [Default: n=1]

    Sets the number of dimensions for the atmosphere in the radiative transfer calculation.

    ========== =======================================================
    Input      Setting
    ========== =======================================================
    1          One-dimensional in altitude
    2          Two-dimensional in altitude and along the line of sight
    ========== =======================================================

.. _tir_raytracingshells:

RayTracingShells
^^^^^^^^^^^^^^^^

.. option:: RayTracingShells (double d) [Default: d=500.0]

   Sets the spacing (in meters) between the spherical shells which divide the atmosphere into homogeneous layers.
   During the ray tracing process the line of sight is divided into segments by determining its intersection with
   these spherical shells. Atmospheric properties are assumed to be constant across each segment. Decreasing the
   spacing will improve accuracy at the cost of increased computation time.

.. _tir_surfaceheight:

SurfaceHeight
^^^^^^^^^^^^^

.. option:: SurfaceHeight (double d) [Default: d=0.0]

    Sets the surface height in meters.

.. _tir_toaheight:

TOAHeight
^^^^^^^^^

.. option:: TOAHeight (double d) [Default: d=100000.0]

    Sets the altitude of the top of the atmosphere. Above this altitude it is assumed there is no atmosphere.
    Due to the assumption of local thermodynamic equilibrium within the Thermal InfraRed engine, increasing
    the value above 100000 m is not recommended because this simplification can not be used at higher altitudes.

.. _tir_userefraction:

UseRefraction
^^^^^^^^^^^^^

.. option:: UseRefraction (int n) [Default: n=0]

    If set to true, then refractive ray tracing is performed for the observer line of sight rays.


Performance
-----------

.. _tir_numthreads:

NumThreads
^^^^^^^^^^

.. option:: NumThreads (int n) [Default: n=0]

    Controls the number of threads to use when multithreading the radiative transfer calculation.  The default value
    of 0 indicates that the number of threads used will be equal to the number of available logical cores on
    the target machine.

    Setting this value to a lower number than the number of available cores can be useful when running a radiative
    transfer calculation in the background and the computer is too slow to multitask.

.. _tir_usecache:

UseCache
^^^^^^^^

.. option:: UseCache (int n) [Default: n=0]

    Designed to enable faster iterative radiance calculations, this setting instructs the TIR engine to compute
    absorption using cached cross sections. This can only be done if:

      1. :meth:`CalculateRadiance()` has been called at least once
      2. The only change to the engine object since the last call to :meth:`CalculateRadiance()` is the number density
         climatology of one or more species added to the engine. The optical properties used for these
         species must the same; if a new optical property is specified, the resulting radiance will be as though
         the old optical property were used.
      3. The species whose climatologies have been changed must also be set as the :option:`WFSpecies` of the engine
         , i.e. only the species for which weighting functions are computed may have their climatologies
         changed


Other
-----

.. _tir_wavelengths:

Wavelengths
^^^^^^^^^^^

.. option:: Wavelengths (array v)

    This option exists to allow the TIR engine to use the :py:class:`sasktran.engine` interface. Unlock other engines,
    TIR requires wavelengths to be set prior to model configuration where the optical properties table is allocated.
    The user does not need to use this option and should instead set wavelengths with the :meth:`EngineTIR.wavelength`
    property when using the ``sasktran`` package or the :meth:`ISKEngine.SetWavelengths` method when using ``sasktranif``.


Diagnostics
-----------

The following options can be used to return information about the engine using the :meth:`ISKEngine.GetProperty` method.
Some parameters such as :ref:`tir_look` and :ref:`tir_observer` require an index to be added as part of
the argument string passed to :meth:`GetProperty`. For example::

    ok, look = engine.GetProperty('Look[0]')  # Get look vector at index 0

.. _tir_look:

Look
^^^^

.. option:: Look (int n) [Default: n=-1]

    Returns the look unit vector of the line of sight at with index ``n``.

.. _tir_numlinesofsight:

NumLinesOfSight
^^^^^^^^^^^^^^^

.. option:: NumLinesOfSight

    Returns the number of lines of sight added to the engine, as a float value.

.. _tir_observer:

Observer
^^^^^^^^

.. option:: Observer (int n) [Default: n=-1]

    Returns the observer position for the line of sight at with index ``n``.

.. _tir_referencepoint:

ReferencePoint
^^^^^^^^^^^^^^

.. option:: ReferencePoint

    Return the :py:class:`GEODETIC_INSTANT` ([latitude, longitude, height, mjd]) which the engine determined from the
    input line(s) of sight and used as the time and location to draw from climatologies to populate the atmospheric
    tables. Can only be called after :meth:`CalculateRadiance`.

.. _tir_tangentalts:

TangentAlts
^^^^^^^^^^^

.. option:: TangentAlts

    Returns a list of tangent altitudes (in meters) for each line of sight added to the engine. Can only
    be called after :meth:`CalculateRadiance`. For lines of sight not looking through the atmospheric limb the
    values returned may be negative or nonsensical.

.. _tir_wavel:

Wavel
^^^^^

.. option:: Wavel

    Returns the wavelengths (in nm) of the most recent radiance calculation.
