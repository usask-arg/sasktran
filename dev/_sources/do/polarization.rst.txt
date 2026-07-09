.. _polarization:

Polarization
************
SK-DO contains the ability to simulate polarized radiances as well as scalar radiances.  Scalar mode is the default,
however polarization can be enabled by setting

.. code-block:: python

    engine.num_stokes = 3


Currently SK-DO assumes that there is no circular polarization (i.e. V=0), which is generally an excellent approximation
within the Earth's atmosphere.