.. _enginehr:


HR
===
The ``HR`` engine solves the radiative transfer equation in a region of interest using the successive orders technique on a true spherical Earth.
The theoretical basis for the technique can be found in [Bourassa2008]_ and upgrades specific to the ``HR`` engine are described in
[Zawada2015]_. The engine is well-suited for calculating limb radiances for **EUV** to **NIR** (nominally 200 nm to 5 micron) wavelengths measured from high altitudes balloons, aircraft and spacecraft
as well as ground-based slant angle observations. It is particularly suitable for conditions where the curvature of the Earth is a significant factor and solar illumination is an important term. The model
features the following highlights:

  * Support to calculate the weighting functions of target species.
  * Adapative integration within optically thick regions.
  * Polarized or scalar radiance calculations.
  * Support for gradients across a region of interest using high resolution diffusive fields.
  * Support for terminator conditions.
  * Support for thermal emissions.
  * Support for elevated surface topology.

Details:

..  toctree::
    :maxdepth: 2

    hr_details

Example:
--------

Here is a simple example that shows basic usage of the radiative transfer engine::

      engine             = skif.ISKEngine('HR')                                                 # Create the engine
      albedo             = 0.3                                                                  # Set the albedo to 0
      usermjd            = 54832.000000000000;                                                  # Set the time of the measurements
      satpos             = [366082.10014671006,  1005802.3038196121,  -6874430.9971755166]      # Set the satellite position
      lookv              = [0.28845686317656621, 0.79252871806432690, 0.53729960834682389]      # Set the look unit-vector
      wavelen            = [294.0,       295.0,       296.0]                                    # Set the wavelengths in nanometers
      rayleigh           = skif.ISKOpticalProperty('RAYLEIGH')                                  # Get Rayleigh scattering optical properties of dry air
      atmosphere         = skif.ISKClimatology    ('MSIS90')                                    # Get a climatology of air number density
      o3_opticalprops    = skif.ISKOpticalProperty('O3_DBM')                                    # Get ozone cross-sections
      o3numberdensity    = skif.ISKClimatology    ('O3LABOW')                                   # Get a climatology of ozone number density

      engine.AddSpecies( 'SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', atmosphere, rayleigh );          # Add dry air to the radiative transfer engine
      engine.AddSpecies( 'SKCLIMATOLOGY_O3_CM3', o3numberdensity, o3_opticalprops);             # add ozone to the radiative transfer engine
      engine.SetAlbedo ( albedo );                                                              # Set the Lambertian albedo of the surface
      engine.SetWavelengths( wavelen );                                                         # Add the three wavelengths

      engine.AddLineOfSight( usermjd,  satpos, lookv);                                          # Add one line of sight, add more if necessary
      ok, radiance = engine.CalculateRadiance();                                                # Calculate the radiance for all lines of sight

It is possible to customize and adjust the performance of the radiative transfer model but
you must specify all of the settings,other than wavelength and species, before calling :meth:`~ISKEngine.CalculateRadiance` or :meth:`~ISKEngine.CalculateStokesVector` as all
settings for each instance, other than wavelength and species, are locked in at this time.
You must create a new instance of the engine if you wish to change any parameters after this time. The
ability to change change wavelengths and species at any time allows the engine to be efficiently adjusted by
atmopsheric retrieval codes.

.. note::

  The ``HR`` model should be used with caution for:

  * Extremely optical thick conditions like clouds. The code does not handle the large optical depths very well.
  * Nadir geometries. They are supported but do not adapt well to the strengths of the ``HR`` model.
  * Scattering and surface BRDF properties with sharp spikes at angles other than forward-scatter are not well supported.


..  py:module:: HR

Properties
----------
All properties are set through::

    ok = engine.SetProperty('PropertyName', value)

The property name is not case sensitive, and value may be a scalar, list, or array.

A complete list of properties is provided for reference.

========================================== =======================
Property                                   Comments
========================================== =======================
aerosolwfsizepercentchange                 undocumented
aerowfmode                                 undocumented
albedoclim                                 undocumented
basis                                      undocumented
:py:func:`~HR.CalcWF`                          -
diagnosticscatterorders                    undocumented
diagnosticdiffuseprofiles                  undocumented
diffuseheightres                           undocumented
:py:func:`~HR.DiffuseIncomingResolution`       -
:py:func:`~HR.DiffuseMaxAngleOffPlane`         -
:py:func:`~HR.DiffusePlacementType`            -
:py:func:`~HR.earthradius`                     -
forcelineardiffuseplacement                undocumented
:py:func:`~HR.ForceOpticalCacheUpdates`        -
:py:func:`~HR.ForceV21DiffuseIncoming`         -
:py:func:`~HR.GroundShiftAlt`                  -
:py:func:`~HR.HorizonSize`                     -
:py:func:`~HR.IntegrationTechnique`            -
look                                       undocumented
:py:func:`~HR.ManualDiffuseHeights`            -
manualdiffuselocations                     undocumented
:py:func:`~HR.ManualOpticalHeights`            -
:py:func:`~HR.manualraytracingshells`          -
:py:func:`~HR.ManualSolarRayTracingShells`     -
:py:func:`~HR.MaxOpticalDepthOfCell`           -
:py:func:`~HR.MinExtinctionRatioOfCell`        -
numlinesofsight                            undocumented
:py:func:`~HR.NumThreads`                      -
:py:func:`~HR.NumDiffuseProfilesInPlane`       -
:py:func:`~HR.NumOrdersOfScatter`              -
:py:func:`~HR.NumDiffuseOutgoing`              -
:py:func:`~HR.NumDiffusePlanes`                -
maxdiffuseheight                           undocumented
observer                                   undocumented
:py:func:`~HR.OpticalAngleGrid`                -
opticalpropertiesheightres                 undocumented
:py:func:`~HR.OpticalNormalAndReference`       -
:py:func:`~HR.OpticalTableType`                -
polarizationtype                           undocumented
precachewavel                              undocumented
:py:func:`~HR.RayTracingCustomShells`          -
:py:func:`~HR.RayTracingShells`                -
referencepoint                             undocumented
:py:func:`~HR.surfaceheight`                   -
:py:func:`~HR.toaheight`                       -
scattermatrixstoragemethod                 undocumented
:py:func:`~HR.setsun`                          -
setreferencepoint                          undocumented
:py:func:`~HR.SolarRayTracingShells`           -
solarspectrum                              undocumented
solartablenumcossza                        undocumented
spectralalbedo                             undocumented
stokesvec                                  undocumented
sun                                        undocumented
tangentalts                                undocumented
:py:func:`~HR.ThreeDOpticalTableParam`         -
useemissions                               undocumented
usesolartransmission                       undocumented
userefraction                              undocumented
:py:func:`~HR.UseShellRayTracer`               -
wavel                                      undocumented
wf                                         undocumented
:py:func:`~HR.WFHeights`                       -
:py:func:`~HR.WFPrecision`                     -
:py:func:`~HR.WFSpecies`                       -
:py:func:`~HR.WFWidths`                        -
========================================== =======================


Geophysical Properties
----------------------

A list of geophysical properties

..  _hr_setsun:

SetSun
^^^^^^
.. py:function:: SetSun (array v) [Default: Determined from MJD]

    Specifies the sun direction in the model.  ``v`` is a unit vector from the Earth to the Sun.  By default the sun
    location is determined from the average MJD of the input lines of sight.

..  _hr_surfaceheight:

SurfaceHeight
^^^^^^^^^^^^^

.. py:function:: SurfaceHeight (double d) [Default: d=0]

    The height of the Earth's surface can be changed to account for topography. This value is the altitude of the surface
    above mean sea level expressed in meters. The model uses this height for the whole of the Earth's surface across the
    ray tracing region used in the model. It is useful for consistently high elevation areas such as as Greenland, Antarctica
    and the Himalayan plateau.

    We currently have not implemented any internal topological maps such as Etopos in to Sasktran and the user must choose
    the surface height. As a guide we note that rays used for limb geometry in the Earth's stratospheric region span
    several degrees of angle around the Earth and this scale size is appropriate for selecting a surface height from
    a topography model.

    This method must be called during the initialization phase of the HR RT model. It cannot be set
    after an instance has called :meth:`~ISKEngine.CalculateRadiance`.

..  _hr_toaheight:

TOAHeight
^^^^^^^^^

.. py:function:: TOAHeight (double d) [Default: d=100000.0]

    The height of the top of the atmosphere in meters. It is assumed there is no atmosphere above this altitude. All climatologies
    used in the radiative transfer model should support valid (not NaN) values up to this height. We have experienced stange
    warning messages during the calculations when this requirement was not met.

    This method must be called during initialization of the HR RT model. It cannot be set after an instance has
    called :meth:`~ISKEngine.CalculateRadiance`.

..  _hr_earthradius:

EarthRadius
^^^^^^^^^^^

.. py:function:: EarthRadius (double d) [Default: Osculating Sphere Calculation]

    Manually sets the radius of the Earth in meters. Normally not required but useful if you want to try and use the same setting
    as another radiative transfer model. This value cannot be set after an instance has called :meth:`~ISKEngine.CalculateRadiance`.


Weighting Functions
-------------------

The HR engine supports analytical computation of weighting functions for absorbing and scattering species. Prior to calculating weighting functions it is
necessary to tell the engine that is should calculate weighting functions, :func:`CalcWF`, and also specify which species
need weighting functions, :func:`WFSpecies`. More details can be found in this description of the HR :func:`hrweightingfunctions`.

The default configuration does not calculate weighting functions.

CalcWF
^^^^^^

.. py:function:: CalcWF (int n)

    [Default: n=0]

    Sets the weighting function calculation mode. For more details see :func:`hrweightingfunctions`.

    ================ ===========
      n              Setting
    ================ ===========
      0              No weighting function calculation
      1              Line of sight weighting function calculation (1D)
      2              User specified (1D)
      3              Two Dimensional weighting function calculation
    ================ ===========


WFPrecision
^^^^^^^^^^^

.. py:function:: WFPrecision (int n)

    [Default: n=0]

    Weighting functions are calculated considering line of sight effects, direct solar effects, and higher order terms.
    Setting this option to 1 disables the direct solar and higher order terms, greatly improving weighting function
    calculation speed at the cost of accuracy.

    ================ ===========
      n              Setting
    ================ ===========
      0              Weighting functions are calculated with full precision/
      1              Weighting functions are calculated neglecting direct solar and higher order effects.
    ================ ===========


WFHeights
^^^^^^^^^

.. py:function:: WFHeights (array v)

    [Default: linspace(500, 99500, 100)]

    See :func:`hrweightingfunctions`.

WFWidths
^^^^^^^^

.. py:function:: WFWidths (array v)

    [Default: ones(100) * 1000]

    See :func:`hrweightingfunctions`.

WFSpecies
^^^^^^^^^

.. py:function:: WFSpecies (ISKOpticalProperty o)

    [Default: None]

    Sets the species to calculate weighting functions for. Species must have scattering or absorbing cross
    sections.  See :func:`hrweightingfunctions` for detailed information on weighting functions.


Diffuse Field
-------------
The primary option that controls the speed and accuracy of the ``HR`` radiative transfer calculation is the
number of *diffuse profiles*.  Each diffuse profile represents a single solar zenith angle where the multiple
scattering source function is calculated.  By default, only one diffuse profile is used in the calculation and
is placed at the solar zenith angle of the average tangent point.  Whether or not this is appropriate is heavily
dependant on a variety of factors. More details at HR :func: `hrdiffuse`,

DiffusePlacementType
^^^^^^^^^^^^^^^^^^^^

.. py:function:: DiffusePlacementType (int n)

    [Default: See Text]

    Changes the placement of diffuse profiles inside the engine.  The default value changes depending on how the
    model is initialized.  It is recommended not to change this setting unless there is good reason to.

    ================ ===========
      n              Setting
    ================ ===========
      0              Diffuse profiles are placed linearly in angle along the line of sight plane
      1              Diffuse profiles are placed linearly in solar zenith angle, forcing one profile at the tangent point
      2              Diffuse profiles are placed linearly in solar zenith angle
    ================ ===========

ForceV21DiffuseIncoming
^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: ForceV21DiffuseIncoming (int n)

    [Default: n=0]

    Changes the spacing of incoming rays on diffuse spheres.  SASKTRANV21 by default had a set of hard coded
    incoming rays, but this could be changed to have higher resolution spacing immediately before the horizon.
    ``HR`` by default has higher resolution incoming rays both above and below the horizon.

    ================ ===========
      n              Setting
    ================ ===========
      0              Incoming rays have high resolution both above and below the horizon
      1              Incoming rays are hard coded to the SASKTRANV21 values
      2              Incoming rays have high resolution only below the horizon
    ================ ===========

NumDiffuseProfilesInPlane
^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: NumDiffuseProfilesInPlane (int n)

    [Default: n=1]

    Sets the number of diffuse profiles in either the solar zenith angle or line of sight plane, depending on the
    setting of :func:`DiffusePlacementType`.  See :func:`hrdiffuse` for more information.

NumDiffuseOutgoing
^^^^^^^^^^^^^^^^^^

.. py:function:: NumDiffuseOutgoing (int n)

    [Default: 169]

    Sets the number of outgoing rays on diffuse spheres.  ``n`` must be a perfect square and is limited to up to
    approximately 1000.

DiffuseIncomingResolution
^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: DiffuseIncomingResolution (array v)

    [Default: See Text]

    Sets the incoming ray resolution for diffuse spheres

    ================ ===========
      Array Index    Setting
    ================ ===========
      0              Number of zenith discretizations before (above) the horizon
      1              Number of zenith discretizations in the horizon region
      2              Number of zenith discretizations below the horizon (towards the ground)
      3              Number of azimuth discretizations in all regions
    ================ ===========

    Default value is ``v=[6, 8, 10, 12]``.

ManualDiffuseHeights
^^^^^^^^^^^^^^^^^^^^

.. py:function:: ManualDiffuseHeights (array v)

    [Default: See Text]

    Sets the altitude for diffuse points.  This setting applies to every profile.  It is recommended to place
    a diffuse point just above the ground.  The default setting for this is::

        v = np.concatenate(([10^-6], np.linspace(500, 99500, 100))]


NumDiffusePlanes
^^^^^^^^^^^^^^^^

.. py:function:: NumDiffusePlanes (int n)

    [Default: n=1] [EXPERIMENTAL]

    Diffuse profiles can also be placed outside the line of sight plane. This option sets the number of planes
    to create, they are spacing based off of :func:`DiffuseMaxAngleOffPlane`.

DiffuseMaxAngleOffPlane
^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: DiffuseMaxAngleOffPlane (double d)

    [Default: d=5.0] [EXPERIMENTAL]

    Sets the maximum angle in degrees that diffuse planes are placed off of the line of sight plane.


HorizonSize
^^^^^^^^^^^

.. py:function:: HorizonSize (double n)

    [Default: n=20]

    Changes the horizon size in degrees for the incoming resolution on diffuse points.  ``HR`` places more incoming rays
    above and below the horizon.


Adaptive Integration
--------------------
Properties that set the adaptive integration parameters.

IntegrationTechnique
^^^^^^^^^^^^^^^^^^^^

.. py:function:: IntegrationTechnique (int n)

    [Default: n=1]

  ``HR`` has added a new line integration technique where the source function is quadratically splined over each
  ray.  Cells with areas of large extinction are also recursively split through an adaptive integration procedure.
  This option allows the integration improvements to be turned off, reverting to the old technique of assuming the source
  function is constant within each cell.  Turning this improvement off offers a small performance increase, however
  unless necessary it is not recommended to do so.

================ ===========
  n              Setting
================ ===========
  0              Uses the legacy SASKTRANV21 integration technique where the source function is assumed to be constant inside each cell.  Adaptive integration is disabled.
  1              The source function is quadratically splined along the line of sight.  Adaptive integration is enabled.
================ ===========

MaxOpticalDepthOfCell
^^^^^^^^^^^^^^^^^^^^^

.. py:function:: MaxOpticalDepthOfCell (double d)

    [Default: n=1000000]

    ``HR`` supports adaptive splitting of cells in the model when certain conditions are met.
    For a cell to be split, the optical depth of the cell must be less than :func:`MaxOpticalDepthOfCell` and
    the ratio of scattering extinction at the start of the cell to the end of the cell must be less than :func:`MinExtinctionRatioOfCell`.
    Note that the high default value means the by default, adaptive integration is essentially disabled.

MinExtinctionRatioOfCell
^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: MinExtinctionRatioOfCell (double d)

    [Default: d=0.9]

    ``HR`` supports adaptive splitting of cells in the model when certain conditions are met.
    For a cell to be split, the optical depth of the cell must be less than :func:`MaxOpticalDepthOfCell` and
    the ratio of scattering extinction at the start of the cell to the end of the cell must be less than :func:`MinExtinctionRatioOfCell`.



Ray Tracing Properties
-----------------------

Properties that affect how the engine traces rays through the atmosphere.

UseShellRayTracer
^^^^^^^^^^^^^^^^^

.. py:function:: UseShellRayTracer (int n)

    [Default: n=0]

    The ``HR`` engine uses a new ray tracer which calculates intersections with a set of arbitrary geometry
    primitives (for one dimensional atmospheres this is a set of spherical shells).  The old ray tracer from
    the ``SO`` engine can be used by calling this option with ``n=1``, however there is no reason to.

    ================ ===========
      n              Setting
    ================ ===========
      0              Uses the standard generic ray tracer
      1              Uses the old shell ray tracer
    ================ ===========

SolarRayTracingShells
^^^^^^^^^^^^^^^^^^^^^

.. py:function:: SolarRayTracingShells (double d)

    [Default: d=1000]

    Sets the vertical shell spacing used for tracing solar rays in [m].  This setting should be changed if higher
    vertical resolution calculations are desired.

ManualSolarRayTracingShells
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: ManualSolarRayTracingShells (array v)

    [Default: See Text]

    Sets custom shell heights for tracing solar rays in [m]. The default is ``v = linspace(0, 100000, 101)``. 
    The lowest shell should be equal to SurfaceHeight (default 0). Only used if SolarRayTracingShells
    is not set.

RayTracingShells
^^^^^^^^^^^^^^^^

.. py:function:: RayTracingShells (double d)

    [Default: d=1000]

    Sets the vertical shell spacing used for tracing rays other than the solar rays in [m].  This setting should
    be changed if higher vertical resolution calculations are desired.

ManualRayTracingShells
^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: RayTracingCustomShells (array v)

    [Default: See Text]

    Sets custom shell heights for tracing non-solar rays in [m]. The default is ``v = linspace(0, 100000, 101)``.
    Only used if RayTracingShells is not set.


Optical Table Properties
------------------------

Properties that affect the internal optical property table.

OpticalTableType
^^^^^^^^^^^^^^^^
.. py:function:: OpticalTableType (int n)

    [Default: n=0]

    Sets the type of atmosphere used internally in the model.  This does not change how the user specifies
    the atmosphere through climatologies, this setting only affects how the model handles the atmosphere it is given.

    ================ ===========
      n              Setting
    ================ ===========
      0              One dimensional atmosphere in height
      1              Three dimensional atmosphere using a Delaunay triangulation.
      2              Two dimensional atmosphere in height and angle along the line of sight plane.
      3              Two dimensional atmosphere in solar zenith angle and height [EXPERIMENTAL]
    ================ ===========

ManualOpticalHeights
^^^^^^^^^^^^^^^^^^^^

.. py:function:: ManualOpticalHeights (array v)

    [Default: See Text]

    Sets the heights in the model where cross sections and number densities are specified.  Typically this
    is twice the resolution of the ray tracing grid.  The default value is ``v = linspace(0, 100000, 201)``.


OpticalNormalAndReference
^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: OpticalNormalAndReference (array v)

    [Default: See Text]

    Only has an effect when :func:`OpticalTableType` is 2.  ``v`` is a 6 element array
    where v[0:3] is the normal and v[3:6] is the reference.  The normal vector is a unit vector normal to the optical
    table plane.  The reference is the x axis of the plane, and the normal cross the reference is the y axis.


OpticalAngleGrid
^^^^^^^^^^^^^^^^

.. py:function:: OpticalAngleGrid (array v)

    [Default: See Text]

    Sets the angular grid used when :func:`OpticalTableType` is set to 2 in degrees.
    The default value is ``v = linspace(-10, 10, 30)``

ThreeDOpticalTableParam
^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: ThreeDOpticalTableParam (array v)

    [Default: See Text]

    When :func:`OpticalTableType` is set to 1, a Delaunay triangulation is constructed around the reference point
    creating the three dimensional optical properties table.  The construction is done through placing ``v[1]`` concentric
    cones around the tangent point spaced ``v[0]`` degrees apart.  Then, ``v[2]`` points are placed uniformly on each
    cone.  The default value is ``v=[1, 10, 10]``.

ForceOpticalCacheUpdates
^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: ForceOpticalCacheUpdates (int n)

    [Default: n=0] [EXPERIMENTAL]

    Debugging option.  If set to 1 the engine forces all climatologies and optical properties to update their
    caches every time the model is evaluated. This is useful for debugging purposes but will slow the code down substantially.


Miscellaneous Properties
------------------------

Miscellaneous properties

NumOrdersOfScatter
^^^^^^^^^^^^^^^^^^

.. py:function:: NumOrdersOfScatter (int n)

    [Default: n=50]

    Sets the number of scatter orders to ``n``.  For more details see :func:`hrscatterorder`.

NumThreads
^^^^^^^^^^

.. py:function:: NumThreads (int n)

    [Default: n=0]

    The number of threads the model will run on.  A setting of ``n=0`` indicates to run on as many threads as the machine
    has logical cores.


Deprecated Properties
---------------------

A list of properties that have been deprecated. The properties are still supported but users are urged to use alternatives

GroundShiftAlt
^^^^^^^^^^^^^^

.. py:function:: GroundShiftAlt (double d)

    [Default: d=0] [DEPRECATED]

    This name is deprecated. Please use SurfaceHeight instead.


References
----------

.. [Bourassa2008] | A.E. Bourassa, D.A. Degenstein, E.J. Llewellyn, SASKTRAN: A spherical geometry radiative transfer code for efficient estimation of limb scattered sunlight, Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 109, Issue 1, January 2008, Pages 52-73, ISSN 0022-4073, http://dx.doi.org/10.1016/j.jqsrt.2007.07.007.

.. [Zawada2015]   | Zawada, D. J., Dueck, S. R., Rieger, L. A., Bourassa, A. E., Lloyd, N. D., and Degenstein, D. A.: High resolution and Monte Carlo additions to the SASKTRAN radiative transfer model, Atmos. Meas. Tech. Discuss., 8, 3357-3397, doi:10.5194/amtd-8-3357-2015, 2015.


