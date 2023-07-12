.. _optical_rayleigh:

RAYLEIGH
========
An optical property object that provides Rayleigh molecular scattering in dry-air. The code closely follows the
algorithm published by Bates 1984 and exactly replicates his cross-section calculations to the 4 significant digits 
in his Table 1. The cross-section is weighted to account for the different gas ratios in standard atmospheric
composition. No attempt is made to track changes in CO2 composition.  Note that water vapour effects are implicitly
ignored as it only considers dry air. The only difference with Bates is that this object takes into account the
tiny fraction of gas that is not N2, O2, Argon or CO2.  Bates ignores this component while this object assumes 
it the residual gas with properties similar to Argon::

	optprop = ISKOpticalProperty('RAYLEIGH');
	[ok, abs,ext,sca] = optprop.CalculateCrossSections( 1.0E7/600.0, [0.0, 0.0, 25000,52393.378298]);

References
^^^^^^^^^^^

1. Rayleigh Scattering by Air: D.R. Bates. Planet. Space. Sci., 32, 6, 785-790, 1984.	
2. The Refractive indices and Verdet constants of the inert gases: A. Dalgarno and A.E. Kingston., Proc.  Roy. Soc. A., 259, 424-429, 1960.
3. Interferometric determination of the refractive index of carbon dioxide in the ultra violet region. A. Bideua-Mehu, Y. Guern, R. Adjean and A. Johannin-Gilles. Optics Communications, 9, 432-434, 1973.
