.. _brdf_li_sparse_kernel:

Li Sparse Kernel
=================

Implements the Li Sparse Kernel BRDF

Example
-------
::

   import sasktranif.sasktranif as skif
   import math

   brdf = skif.ISKBrdf('LI_SPARSE_KERNEL')
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)

