.. SASKTRAN documentation master file, created by
   sphinx-quickstart on Tue Feb  7 14:08:28 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The SASKTRAN Radiative Transfer Framework
=========================================
The SASKTRAN radiative transfer framework is a radiative transfer tool developed at the University of Saskatchewan.
Originally designed for use with the OSIRIS instrument (https://research-groups.usask.ca/osiris/) it has since
evolved to be applicable to a large variety of applications.  SASKTRAN is a full framework and not just a radiative
transfer model, as such it contains databases or interfaces to standard climatologies and species optical properties.
It also contains multiple methods to solve the radiative transfer equation that we call `engines`, each one of which has a
specific domain that it is applicable to.  The primary two engines are ``SASKTRAN-HR`` which is a fully spherical
scattering model and contained within this package, and ``SASKTRAN-DO`` which is a discrete ordinates model suitable for
nadir viewing applications.  More information on ``SASKTRAN-DO`` can be found at https://arg.usask.ca/docs/sasktran_do/.

Development
-----------
The SASKTRAN code is open-source and is in active development at https://github.com/usask-arg/sasktran

License
-------
SASKTRAN is made available under the MIT license restricted for academic and educational use, for commercial use please contact the package authors.
Commerical level support may also be available for specific applications.
We request that users of the model contact the authors before publishing results using SASKTRAN, and that
the following publications are acknowledged:

Zawada, D. J., Dueck, S. R., Rieger, L. A., Bourassa, A. E., Lloyd, N. D., and Degenstein, D. A.: High-resolution and Monte Carlo additions to the SASKTRAN radiative transfer model, Atmos. Meas. Tech., 8, 2609-2623, https://doi.org/10.5194/amt-8-2609-2015, 2015.

Bourassa, A. E., Degenstein, D. A., and Llewellyn, E. J.: SASKTRAN: A Spherical Geometry Radiative Transfer Code for Efficient Estimation of Limb Scattered Sunlight, J Quant Spectrosc Radiat Trans, Volume 109, Issue 1, 52-73, https://doi.org/10.1016/j.jqsrt.2007.07.007, 2008.


Contact
-------
Inquiries and bug reports can be posted at https://github.com/usask-arg/sasktran.


.. toctree::
   :maxdepth: 2

   installing
   quickstart
   constructing_atmosphere
   constructing_geometry
   other_engines
   user_aerosol
   examples
   config
   api_reference
   advanced
   experimental
   changelog

Specific Engine Information
---------------------------
.. toctree::
   :maxdepth: 2
   :glob:

   sasktran_hr
   tir/index
   do/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

