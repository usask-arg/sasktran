.. _configuration:

Configuration
*************
Some settings, particularly the location of large databases, are stored in a global configuration file.  The
location of this file on your system can be found by running::

    import sasktran as sk

    print(sk.config.user_config_file_location())

The file is a standard `yaml` format file that can be edited with any text editor.  The format of the file is
`key: value`, you will probably have to create the file the first time.


Databases
=========
Several databases can be downloaded using functions within the `sasktran` package.  Note that after downloading
any of these databases it might be necessary to restart the interpreter for the settings to take effect.

PRATMO
^^^^^^
The :class:`sasktran.Pratmo` climatology object requires a ~50 mb file to be downloaded to use which can be done
by running::

    import sasktran as sk

    sk.config.download_pratmo()


HITRAN
^^^^^^
In order to use HITRAN you must download a large (few gb) database with::

    import sasktran as sk

    sk.config.download_hitran()


This might take upwards of 20 minutes.


Key Reference
=============

hitran_directory
^^^^^^^^^^^^^^^^
Location of the HITRAN 2016 database.  Necessary to use the :class:`sasktran.HITRANChemical` object. Folder
should contain `By-Molecules` and `Global_Data` directories which contains the database.  The standard value for this key on the
University of Saskatchewan network is `\\\\utls\\SasktranFiles\\hitranhapi`

hitran_cache
^^^^^^^^^^^^
Location of a HITRAN cache folder.  Should be set by default, but can be changed to a different location
if desired.

baum_directory
^^^^^^^^^^^^^^
Location of the Baum 2014 ice crystal optical property database.  Necessary to use the :class:`sasktran.BaumIceCrystal`
object.  Folder should contain the file `GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc`. The
standard value for this key on the University of Saskatchewan network is `\\\\utls\\SasktranFiles\\BaumIceCrystals`

pratmo_directory
^^^^^^^^^^^^^^^^
Location of a folder containing a climatology based upon the Pratmo photochemical box model.
Necessary to use the :class:`sasktran.Pratmo` climatology object.  This file can be directly downloaded with::

    import sasktran as sk

    sk.config.download_pratmo()

which will automatically download the file and place it in an appropriate user location.  If you run this command
it is not necessary to set this key.

glossac_file
^^^^^^^^^^^^
Full path to a NetCDF file containing the GloSSAC climatology.  Necessary to use the
:class:`sasktran.GloSSAC` climatology object.
