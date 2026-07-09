.. _other_engines:

Choosing an Engine
******************
SASKTRAN is a framework containing modules to handle all aspects of the radiative transfer problem, including setting up the geometry, defining the atmosphere, and solving the radiative transfer equation.
The radiative transfer equation can be solved in multiple ways, and SASKTRAN contains multiple engines, or solvers, to perform this task.
Which engine to use depends on the specific task at hand.

SASKTRAN-HR
===========
SASKTRAN-HR is the primary engine used and is described in Zawada et. al. (2015).

Suitable For
^^^^^^^^^^^^
Problems requiring fast and accurate radiative transfer in the UV-VIS-NIR spectral region where scattering is dominant.

Advantages
^^^^^^^^^^
* Fully spherical geometry
* Fast
* Approximate, analytic Jacobian calculation
* Support for two- and three-dimensionally varying atmospheres

Disadvantages
^^^^^^^^^^^^^
* Accuracy is controlled by a wide variety of settings which may require specific tuning for your problem.
* Support for atmospheric emissions is available, but experimental and limited.


SASKTRAN-OCC
============
SASKTRAN-OCC is an engine designed for fast occultation applications.


Suitable For
^^^^^^^^^^^^
Occultation viewing geometry in all spectral regions.

Advantages
^^^^^^^^^^
* Extremely fast
* Refraction is treated within a fully spherical atmosphere

Disadvantages
^^^^^^^^^^^^^
* Scattering and emission terms are neglected

SASKTRAN-MC
===========
SASKTRAN-MC performs Monte-Carlo integration to solve the RTE and is described in Zawada et. al. (2015).

Suitable For
^^^^^^^^^^^^
Problems requiring accurate radiative transfer in the UV-VIS-NIR spectral region where scattering is dominant and runtime is of secondary concern.

Advantages
^^^^^^^^^^
* Fully spherical geometry
* Accurate, requiring much less configuration than SASKTRAN-HR

Disadvantages
^^^^^^^^^^^^^
* Can be slow depending on the application 
* Results contain random noise

SASKTRAN-DO
===========
SASKTRAN-DO is a discrete ordinates solver based upon the VLIDORT technique.

Suitable For
^^^^^^^^^^^^
Nadir viewing applications in the UV-VIS-NIR spectral region.

Advantages
^^^^^^^^^^
* Extremely fast 
* Perfectly linearized for most inputs

Disadvantages
^^^^^^^^^^^^^
* Not suitable for limb viewing applications
* Errors that increase as the viewing zenith angle increases.

SASKTRAN-TIR
============
SASKTRAN-TIR is a thermal radiative transfer model designed for applications in the IR and beyond spectral regions.

Suitable For
^^^^^^^^^^^^
Applications in the IR and beyond spectral regime where scattering can be neglected.

Advantages
^^^^^^^^^^
* Fast
* Spherical geometry
* Linearized for most atmospheric parameters 

Disadvantages
^^^^^^^^^^^^^
* Scattering is not included 

References
==========
Zawada, D. J., Dueck, S. R., Rieger, L. A., Bourassa, A. E., Lloyd, N. D., and Degenstein, D. A.: High-resolution and Monte Carlo additions to the SASKTRAN radiative transfer model, Atmos. Meas. Tech., 8, 2609-2623, https://doi.org/10.5194/amt-8-2609-2015, 2015.

Bourassa, A. E., Degenstein, D. A., and Llewellyn, E. J.: SASKTRAN: A Spherical Geometry Radiative Transfer Code for Efficient Estimation of Limb Scattered Sunlight, J Quant Spectrosc Radiat Trans, Volume 109, Issue 1, 52-73, https://doi.org/10.1016/j.jqsrt.2007.07.007, 2008.