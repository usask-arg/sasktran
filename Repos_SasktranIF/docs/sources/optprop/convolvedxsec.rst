.. _optical_convolvedxsec:

CONVOLVED_CROSSSECTION
=======================

A special optical property object that allows a user to convolve another optical property object 
to a lower resolution. The object will normally convolve and cache optical properties as a function
of position as many cross-section objects have different cross-sections for different temperatures and pressures.

Users should try to set the base cross-section object early in the life-cycle of the CONVOLVED_CROSSSECTION object.
We recommend that new CONVOLVED_CROSSSSECTON objects are created for new base cross-section objects rather than
try to reuse the same object.

Note that many of our :class:`ISKOpticalProperty` objects know the resolution and sampling properties of the
instrument that measured them and these are taken into account during the convolution process.

Properties
^^^^^^^^^^
.. option::SetFWHM ( double fwhm)
    Set the resolution of the convolved cross-section. The convolution is too a Gaussian whose FWHM is specified in nanometers.

.. option::SetBaseCrossSection( ISKOpticalProperty object)
    Set the base cross-section object. This object is an instance of :class:`ISKOpticalProperty`. 
