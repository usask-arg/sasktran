.. _brdf:

****
Brdf
****

The BRDF objects provide the surface reflectance of the ground. The interface currently supports scattering from one incoming
direction to another. It is a function of the incoming zenith angle, outgoing zenith angle and the azimuthal difference between the
the two directions.

==============================================      ======================  =====================================================
BRDF                                                Extension               Description
==============================================      ======================  =====================================================
:ref:`brdf_lambertian`                              :ref:`sasktran_core`    Lambertian BRDF
:ref:`brdf_cox_munk`                                :ref:`sasktran_core`    Cox-Munk ocean surface BRDF.
:ref:`brdf_kokhanovsky`                             :ref:`sasktran_core`    Snow BRDF
:ref:`brdf_hapke`                                   :ref:`sasktran_core`    Hapke BRDF
:ref:`brdf_rahman`                                  :ref:`sasktran_core`    Rahman BRDF
:ref:`brdf_modis`                                   :ref:`sasktran_core`    Modis BRDF
:ref:`brdf_roujean`                                 :ref:`sasktran_core`    Roujean BRDF
:ref:`brdf_roujean_kernel`                          :ref:`sasktran_core`    Roujean Kernel BRDF
:ref:`brdf_li_sparse_kernel`                        :ref:`sasktran_core`    Li Sparse Kernel BRDF
:ref:`brdf_li_dense_kernel`                         :ref:`sasktran_core`    Li Dense Kernel BRDF
:ref:`brdf_li_sparse_recip_kernel`                  :ref:`sasktran_core`    Li Sparse Reciprocal Kernel
:ref:`brdf_ross_thin_kernel`                        :ref:`sasktran_core`    Ross Thin Kernel
:ref:`brdf_ross_thick_kernel`                       :ref:`sasktran_core`    Ross Thick kernel
:ref:`brdf_linearcombo`                             :ref:`sasktran_core`    Linear combination of two BRDF's
==============================================      ======================  =====================================================

..  toctree::
    :maxdepth: 2

    brdf/brdf_lambertian
    brdf/cox_munk
    brdf/brdf_snow_kokhanovsky
    brdf/hapke
    brdf/rahman
    brdf/modis
    brdf/brdf_roujean
    brdf/roujean_kernel
    brdf/li_sparse_kernel
    brdf/li_dense_kernel
    brdf/li_sparse_recip_kernel
    brdf/ross_thin_kernel
    brdf/ross_thick_kernel
    brdf/linearcombo
