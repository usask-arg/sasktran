.. _sasktran_core_solar:

..  _solar_sao2010:

SAO2010
=======
Implements the SAO2010 solar spectrum reported by Chance & Kurucz in 2010.
The spectrum is a high quality. high resolution solar spectrum going from
200 nm to 1 micron at 0.01 nm spacing and a resolution of 0.04 nm FWHM.


Example
-------
::

   import sasktranif.sasktranif as skif
   sun = skif.ISKSolarSpectrum('SAO2010')

References
-----------
**K. Chance** and R.L. Kurucz, "An improved high-resolution solar reference spectrum for earth's atmosphere measurements in the ultraviolet, visible, and near infrared ",
*Journal of Quantitative Spectroscopy and Radiative Transfer*, 111, 9, pp. 1289-1295, 2010, doi: 10.1016/j.jqsrt.2010.01.036,
url = `http://www.sciencedirect.com/science/article/pii/S0022407310000610 <http://www.sciencedirect.com/science/article/pii/S0022407310000610>`_

..  _solar_fontela_uvis_3micron:

FONTELA_UVIS_3MICRON
====================
Implements the solar spectrum presented by Fontela et al 2001. We use the 200 nm to 3 micron table
at 0.02 nm spacing with 0.1 nm resolution., file -> NUVisIr3irrad1001Lowion00.conv1a.fits.

Example
-------
::

   import sasktranif.sasktranif as skif
   sun = skif.ISKSolarSpectrum('FONTELA_UVIS_3MICRON')

..  _solar_fontela_uvis_100micron:

References
----------
**Fontenla, J. M.**, Harder, J., Livingston, W., Snow, M. and Woods, T. "High-resolution solar spectral irradiance from extreme ultraviolet to far infrared", *Journal of Geophysical Research: Atmospheres*, 116, D20, pp 2156-2202, 2011
doi: 10.1029/2011JD016032, url: `http://dx.doi.org/10.1029/2011JD016032 <http://dx.doi.org/10.1029/2011JD016032>`_

FONTELA_UVIS_100MICRON
======================
Implements the solar spectrum presented by Fontela et al 2001. We use the 200 nm to 100 micron table
at 0.2 nm spacing with 1 nm resolution., file -> NUVisIr100irrad1001Lowion00.conv1nm.fits.

Example
-------
::

   import sasktranif.sasktranif as skif
   sun = skif.ISKSolarSpectrum('FONTELA_UVIS_100MICRON')

References
----------
**Fontenla, J. M.**, Harder, J., Livingston, W., Snow, M. and Woods, T. "High-resolution solar spectral irradiance from extreme ultraviolet to far infrared", *Journal of Geophysical Research: Atmospheres*, 116, D20, pp 2156-2202, 2011
doi: 10.1029/2011JD016032, url: `http://dx.doi.org/10.1029/2011JD016032 <http://dx.doi.org/10.1029/2011JD016032>`_