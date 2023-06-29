.. _ISKStokesVector:

****************
ISKStokesVector
****************
The ISKStokesVector is used to store polarized radiance. This includes the 4
elements of radiance,

1. **I**
2. **Q**
3. **U**
4. **V**

and also includes 3 vectors that define the basis/reference frame for the polarization,

1. The propagation vector parallel to the direction of the ray.
2. The theta direction perpendicular to the ray direction.
3. The phi direction perpendicular to both thea and the ray propagation.

.. py:class:: ISKStokesVector()

   Create a new empty instance os ISKStokesVector(). Yu should never need to create your own instances of
   ISKSTokesVector, if you do you will find there are no methods available to assign them meaningful values.
   Instances of ISKStokesVector are normally returned by calls to :meth:`ISKEngine::CalculateStokesVector`

I
^^^^
.. py:method:: ISKStokesVector.I()

      Returns the I component of the polarized radiance::
      
         I = stokes.I()
      
      :return: the I component of polarized radiance. Always a scalar value

Q
^^^^
.. py:method:: ISKStokesVector.Q()

      Returns the Q component of the polarized radiance::
      
         Q = stokes.Q()
      
      :return: the Q component of polarized radiance. Always a scalar value

U
^^^^^
.. py:method:: ISKStokesVector.U()

      Returns the U component of the polarized radiance::
      
         U = stokes.U()
      
      :return: the U component of polarized radiance. Always a scalar value

V
^^^^^
.. py:method:: ISKStokesVector.V()

      Returns the V component of the polarized radiance::
      
         V = stokes.V()
      
      :return: returns the V component of polarized radiance. Always a scalar value

propagation_direction
^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKStokesVector.propagation_direction()

      Returns the 3 element vector parallel to the ray propagation::
      
         prop = stokes.propagation_direction()
      
      :return: the 3 element vector (X,Y,Z) where the 3 vectors are expressed in
               a geographic, geocentric system: origin at centre of Earth, X in equatorial plane
               points to Greenwich meridian, Z parallel to spin axis and Y in equatorial plane points to
               90 degrees East.


theta_direction
^^^^^^^^^^^^^^^
.. py:method:: ISKStokesVector.theta_direction()

      Returns the 3 element vector perpendicular to the ray propagation::
      
         theta = stokes.theta_direction()
      
      :return: the 3 element vector (X,Y,Z) where the 3 vectors are expressed in
               a geographic, geocentric system: origin at centre of Earth, X in equatorial plane
               points to Greenwich meridian, Z parallel to spin axis and Y in equatorial plane points to
               90 degrees East.

phi_direction
^^^^^^^^^^^^^
.. py:method:: ISKStokesVector.phi_direction()

      Returns the 3 element vector perpendicular to the ray propagation::
      
         phi = stokes.phi_direction()
      
      :return: the 3 element vector (X,Y,Z) where the 3 vectors are expressed in
               a geographic, geocentric system: origin at centre of Earth, X in equatorial plane
               points to Greenwich meridian, Z parallel to spin axis and Y in equatorial plane points to
               90 degrees East.

to_new_basis
^^^^^^^^^^^^
.. py:method:: ISKStokesVector.to_new_basis( prop, theta, phi) -> ok
   
      Rotates the current polarization from the current reference frame to a new reference frame
      given by the 3 vectors **prop**, **theta** and ** phi**. The polarized intensities I,Q,U and V are
      transformed to trhis new coordinate system::
      
         ok = stokes.to_new_basis( prop, theta, phi)
   
      :param array prop:
         A 3 element array of doubles. Defines the new propagation direction.
         
      :param array theta:
         A 3 element array of doubles. Defines the new theta direction.

      :param array phi:
         A 3 element array of doubles. Defines the new phi direction.
         
      :param boolean ok:
         returns true if successful

      :return: returns true if successful
