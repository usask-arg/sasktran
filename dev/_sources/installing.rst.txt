
.. _installing:

Installing SASKTRAN
*******************
SASKTRAN is designed to work with the Anaconda scientific Python distribution (https://www.anaconda.com/products/distribution), but it may work with your system Python distribution as long as the necessary dependencies are also available.
Currently we support Python versions 3.8, 3.9, 3.10, 3.11 on Windows, Linux, and Mac (x86_64) environments, as well as experimental support for Mac ARM (m-series processors).

SASKTRAN can be installed through::

   pip install sasktran

This will install this package as well as the C++ code.

.. important::

  If you are on an ARM Mac you must install SASKTRAN through::

    pip install sasktran -f https://arg.usask.ca/wheels/

  instead.  Note that these wheels are currently experimental.

Mac Issues
==========
Some incompatibilities are known to exist with other Python packages that use a different OMP runtime library, notably MKL based numpy builds.
To create a new conda environment that is known to be compatible with SASKTRAN you can use::

    conda create -n ENVNAME -c conda-forge python=3.11 numpy nomkl xarray

and then::

    conda activate ENVNAME
    pip install sasktran

where you can change ENVNAME to a name of your choosing, as well as use any supported version of Python.
