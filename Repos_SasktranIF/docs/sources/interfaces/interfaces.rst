
.. _interfaces:

*********************
SasktranIF Interfaces
*********************

Radiative transfer calculations within the SasktranIF framework consist of the following core conceptual components,

==================================  =====================================================
Class                               Description
==================================  =====================================================
:class:`ISKOpticalProperty`         Optical properties of individual atoms/molecules or particles present in the atmosphere.
:class:`ISKClimatology`             Number density and other properties of different species at different points in the atmosphere.
:class:`ISKEngine`                  Mathematical techniques to trace rays through the described atmosphere and calculate observed radiances.
:class:`ISKBrdf`                    Implements surface reflectance functions
:class:`ISKSolarSpectrum`           Conversion of calculations to an observed radiance needs an accurate solar spectrum at all wavelengths/wavenumbers.
:class:`ISKEmission`                Background radiating brightness of different locations in the atmosphere due to thermal or photo-chemical emissions.
:class:`ISKGeodetic`                Specify observers, lines of sight, tangent points and look directions we need a model of the oblate Earth.
==================================  =====================================================

..  toctree::
    :maxdepth: 2

    iskopticalproperties
    iskclimatology
    iskengine
    isksolarspectrum
    iskstokesvector
    iskbrdf
    iskemission
    iskgeodetic



