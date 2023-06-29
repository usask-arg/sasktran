..  _sasktran_core:

sasktran_core
--------------
The `sasktran_core` extension is provided with the `sasktran` package. It provides a variety of sasktran objects
that allow most users to perform precision atmospheric radiative transfer calculations out-of-the-box.

Configuration
^^^^^^^^^^^^^
The `sasktran_core` extension has two configuration aspects. First, it must specify folders for caching of information
generated as objects execute their code. Default values based upon the users profile are chosen for the caches during installation but the user may
wish to review and change these locations. Users should ensure there are several gigabytes of free space for caching activity.

Second, Some sasktran objects within `sasktran_core` use databases which, due to size considerations, are not distributed with the
package. These databases must be downloaded by the user if they wish to use these objects.  Currently `sasktran_core` objects may use
the following database

    * Hitran spectral line database.
    * Pratmo climatology
    * GRIB definitions and samples (only required on Windows systems)

The `sasktran_core` package provides a method to download these databases which,

    #. fetches the Hitran 160 column database for every molecule from the `Hitran <https://hitran.org>`_ website using their HAPI interface.
    #. fetches the pratmo climatology binary file from the `ARG <https://arg.usask.ca>`_ website. This is currently a low level binary file
       and will only work on Windows and Intel based Linux. Eventually we will convert over netcdf or similar.
    #. fetches the GRIB definitions used for proper use of the GRIB file formats.  This is only required for Windows operating systems.

The method can be invoked in python with sasktran_core.update_settings.Modify()::

    import sasktran
    from sasktran_core.update_settings import modify

    modify()

The code will place up a simple menu which the user can use to modify settings.  Note that the settings may need to be re-entered if
the sasktran_core package is re-installed or upgraded (sorry).



