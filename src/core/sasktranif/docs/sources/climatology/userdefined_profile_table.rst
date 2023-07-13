.. _clim_userdefined_profile_table:

USERDEFINED_PROFILE_TABLE
=========================
A climatology class that allows users to define their own height profiles of scalar
parameters.  The class is very similar in functionality to climatology :ref:`clim_userdefined_profile`
but wraps an older C++ implementation that was used for the OSIRIS MART retrievals. This class
provides a useful property for reading height profiles from two column text files.

This climatology only provides linear interpolation in altitude. We recommend that users 
only use this class for backward compatibility. New uses should employ the :ref:`clim_userdefined_profile`
class.


Python extension
----------------
The USERDEFINED_PROFILE_TABLE climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The USERDEFINED_PROFILE_TABLE climatology needs no external preparation.

Properties
----------

.. py:module:: USERDEFINED_PROFILE_TABLE

Climatology_Handle
^^^^^^^^^^^^^^^^^^
.. py:function:: Climatology_Handle ( handle: str )

    Sets the Climatology Handle that will be used for the profile loaded in the in the next call to property FromTextFile.
    The current implemenattion requires that the user use a climatology_handle that is already defined within the Sasktran framework.
    Future versions may ease this restriction.

FromTextFile
^^^^^^^^^^^^
.. py:function:: FromTextFile( filename: str )

    Loads a height profile from the given file. It is assumed that the file is organized as a two column array of numbers. The first
    column is height above sea level in meters. The second column is the height profile of the given species. The species is defined by the 
    last call to Property 'Climatology_Handle'.

Heights
^^^^^^^
.. py:function::   Heights

    Sets the height grid that will be used in subsequent calls to :meth:`~ISKClimatology.SetPropertyUserDefined`.
    The heights must be in ascending order and specify the height of the grid point above sea level in meters.
    The default is an empty array.
    
 
