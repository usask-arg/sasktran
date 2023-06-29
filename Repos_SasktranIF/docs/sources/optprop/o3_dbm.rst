.. _optical_o3dbm:


O3_DBM
======
Tabulated high resolution cross-sections of O3 measured by Daumont, Brion and Malicet in the early
1990's (Daumont et al. 1992). The wavelength range slightly varies with temperature but covers the
entire UV to NIR region, from 194.50 nm to 830.00 nm at 0.01 to 0.02 nm resolution. The cross-section 
data were collected at 0.01-0.02 nm resolution and each wavelength/cross-section table varies in size
from 22,052 to 63,501 entries. The data consists of 5 tables of wavelength versus cross-section 
for 5 temperatures.

Temperature Range
^^^^^^^^^^^^^^^^^
Measurements are provided at 5 temperatures covering typical stratospheric and tropospheric conditions:

* 218 K
* 228 K
* 243 K
* 273 K
* 295 K

Wavelength Range
^^^^^^^^^^^^^^^^
The wavelength range of each temperature table is slightly different and is given below. Note that most of
the temperature variation occurs in the Huggins band between 315 and 360 nm,

* 218K -> 194.50nm to 650.01nm
* 228K -> 194.50nm to 520.01nm
* 243K -> 194.50nm to 519.01nm
* 273K -> 299.50nm to 520.01nm
* 295K -> 195.00nm to 830.00nm

We looked into temperature interpolation and while DBM suggest that a quadratic interpolation scheme (page 269 of paper 3.)
they do not indicate an explicit technique. We tested several quadratic fitting routines and found that a truncated linear 
fit in temperature was visually more appealing than any of the quadratic fits and had none of the undesirable artifacts 
(excessive curvature etc.) that naturally arise with quadratic curve fitting. Consequently this object uses a truncated linear fit
in temperature.

Properties
^^^^^^^^^^
.. option:: SetTemperature (double n)
   
   Sets the temperature in Kelvins that will be used in the next calculation of cross-sections.

Data Source
^^^^^^^^^^^
These data are an exact replication of the data files,

* O3_CRS_BDM_218K.dat
* O3_CRS_BDM_228K.dat
* O3_CRS_BDM_243K.dat
* O3_CRS_BDM_273K.dat
* O3_CRS_BDM_295K.dat

on the IGACO site, http://igaco-o3.fmi.fi/ACSO/files/cross_sections. The files were copied on July 16-July 25 2012.

References
^^^^^^^^^^

1. Daumont Brion Malicet. O3 cross-sections as published in, D. Daumont et al. Ozone UV spectroscopy I: Absorption cross-sections at room temperature Journal of Atmospheric Chemistry Vol. 15, 1992.
2. Brion J. et al. High-resolution laboratory absorption cross section of O3. Temperature effect. Chemical Physics Letters Vol. 213, No. 5,6, 1993.
3. Malicet J. et al. Ozone UV spectroscopy II: Absorption cross-sections and temperature dependence.  Journal of Atmospheric Chemistry Vol. 21, 1995.
4. Brion J. et al. Absorption spectra measurements for the ozone molecule in the 350-830 nm region.  Journal of Atmospheric Chemistry Vol. 31, 1998.
