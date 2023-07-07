# sasktran
The SASKTRAN radiative transfer framework is a radiative transfer tool developed at the University of Saskatchewan.
Originally designed for use with the OSIRIS instrument (https://research-groups.usask.ca/osiris/) it has since
evolved to be applicable to a large variety of applications.  SASKTRAN is a full framework and not just a radiative
transfer model, as such it contains databases or interfaces to standard climatologies and species optical properties.
It also contains multiple methods to solve the radiative transfer equation that we call `engines`, each one of which has a
specific domain that it is applicable to.  The primary two engines are ``SASKTRAN-HR`` which is a fully spherical
scattering model and contained within this package, and ``SASKTRAN-DO`` which is a discrete ordinates model suitable for
nadir viewing applications.  More information on ``SASKTRAN-DO`` can be found at https://arg.usask.ca/docs/sasktran_do/.

## Installation
SASKTRAN is designed to work with the Anaconda scientific Python distribution (https://www.anaconda.com/products/distribution), but it may work with your system Python distribution as long as the necessary dependencies are also available.
Currently we support Python versions 3.8, 3.9, 3.10, 3.11 on Windows and Linux x64 environments.
Mac wheels for x86_64 are also built for the same Python versions, but support is considered experimental, and problems are known to exist.
Mac ARM builds are planned, but not currently available.

Installing SASKTRAN can be done with

```
pip install sasktran
```

This will install this package as well as the C++ code.

## Usage
Detailed documentation can be found at https://usask-arg.github.io/sasktran/.
The primary way to use SASKTRAN is through the Python interface.


## Development
The SASKTRAN code is open-source and is in active in this github repository (https://github.com/usask-arg/sasktran)

## License
SASKTRAN is made available under the MIT license subject to the Commons Clause condition (see license.md).
Effectively this is a  MIT license restricted for academic and educational use, for commercial use please contact the package authors.
Commerical level support may also be available for specific applications.

## Acknowledgement
We request that users of the model contact the authors before publishing results using SASKTRAN, and that
the following publications are acknowledged:

Zawada, D. J., Dueck, S. R., Rieger, L. A., Bourassa, A. E., Lloyd, N. D., and Degenstein, D. A.: High-resolution and Monte Carlo additions to the SASKTRAN radiative transfer model, Atmos. Meas. Tech., 8, 2609-2623, https://doi.org/10.5194/amt-8-2609-2015, 2015.

Bourassa, A. E., Degenstein, D. A., and Llewellyn, E. J.: SASKTRAN: A Spherical Geometry Radiative Transfer Code for Efficient Estimation of Limb Scattered Sunlight, J Quant Spectrosc Radiat Trans, Volume 109, Issue 1, 52-73, https://doi.org/10.1016/j.jqsrt.2007.07.007, 2008.

## Compiling the code
To compile the code you need the following dependencies that are findable within cmake:

```
boost
netcdf
catch2
yaml-cpp
blas/lapack (recommended to use OpenBLAS)
```

Then the code can be compiled and installed with

```
cmake -S . -B build  -DPythonEnvironments="current" -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --target install
```

Which will build the `c++` code and create a Python wheel for the current Python interpreter on the path.