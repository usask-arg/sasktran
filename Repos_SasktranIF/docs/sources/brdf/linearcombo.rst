.. _brdf_linearcombo:

Linear Combination
==================
This BRDF allows an arbitrary number of BRDF kernels to be used in 	a linear combination.

Example
-------
::

   import sasktranif.sasktranif as skif
   import math

   brdf = skif.ISKBrdf('LINEAR_COMBINATION')
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)

Properties
-----------
..  py:module:: BRDF_LINEAR_COMBINATION

AddKernel
^^^^^^^^^
..  py:function:: AddKernel( object kernel)

    Appends the given BRDF object to the current list of kernels. A linear weight can be applied to the kernel by calling
    :py:func:`~BRDF_LINEAR_COMBINATION.KernelWeights` after all the kernels have been added.

KernelWeights
^^^^^^^^^^^^^
..  py:function:: KernelWeights(array weights)

    Sets the linear weight for each element in the current list of kernels.

RemoveKernel
^^^^^^^^^^^^
..  py:function:: RemoveKernel( int index )

    Removes the given kernel from the current list of kernels.

    :param int index: The zero based index of the requested kernel.
