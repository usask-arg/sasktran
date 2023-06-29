.. _advanced:

Advanced Topics
***************

Creating SasktranIF Like Objects
================================

Motivation
----------
Often it is desired to create a new optical property or climatology for a project, but modifying the SasktranIF
internals is not an option.  In the climatology case solution for this kind of problem is to use the UserDefined classes available
to make a climatology for example.  However this has the downside that the created object no longer acts like a
`Climatology`, it is often fixed to an instant in time and location.  With the `sasktran` package we are able to create
pseudo-sasktranif objects where the logic is implemented in Python, but are still able to be used in the same way as
regular `sasktranif` style objects.  This allows a decoupling of the configuration and execution stages.

How It Works
------------
All climatologies, optical properties, etc, that are input to a SASKTRAN engine must still be raw `sasktranif` objects.
However since in the `sasktran` package we operate one level above the raw `sasktranif` interface, only converting
to `ISK` style objects when necessary, we can defer this conversion until the last possible moment, typically right
before the radiative transfer calculation.  At the time of this conversion we ideally have a variety of information
available to us, for example, the wavelengths of the calculation, and the internal model reference point.  Therefore
any class that has information that can be converted into a `sasktranif` base object, such as a `UserDefinedProfile` can
be integrated very closely with the `sasktran` package.

Example - Scaling Climatology
-----------------------------
As a first example, we will implement a scaling climatology, something which takes in an existing climatology and scales
it by a constant value. ::

    from sasktran.climatology import ClimatologyBase
    import numpy as np

    class ScalingClimatology(ClimatologyBase):
        def __init__(self, other_climatology, scale_factor):
            self._other_climatology = other_climatology
            self._scale_factor = scale_factor

        def get_parameter(self, species: str, latitude: float, longitude: float, altitudes: np.ndarray, mjd: float):
            return self._other_climatology.get_parameter(species, latitude, longitude, altitudes, mjd) * self._scale_factor

        def skif_object(self, **kwargs):
            pass

All we do is call the other climatology, and scale it by a number.  Let's try it out::

    from sasktran.climatology import Labow

    other_clim = Labow()
    scaled_clim = ScalingClimatology(other_clim, scale_factor=2.0)

    other_clim.get_parameter('SKCLIMATOLOGY_O3_CM3', 0, 0, [10000, 30000], 54372)
    # array([  4.78059418e+11,   3.66810991e+12])
    scaled_clim.get_parameter('SKCLIMATOLOGY_O3_CM3', 0, 0, [10000, 30000], 54372)
    # array([  9.56118837e+11,   7.33621983e+12])


Now we have a climatology that we can get parameters out of it, but can we can't actually use this in a radiative
transfer calculation.  The climatology is passed through to the SASKTRAN engine by using the `skif_object` method,
which we have not yet implemented.  In this function we need to convert our Climatology to a form that SASKTRAN
can use.  This can be done in a variety of ways, but lets take advantage of an existing class, the ClimatologyUserDefined
class. We also need to know what the reference point for the model is, this is where the \*\*kwargs argument saves the day.
All standard engines when calling skif_object will pass a reference to themselves inside the \*\*kwargs argument. ::

    from sasktran.climatology import ClimatologyBase, ClimatologyUserDefined
    import numpy as np

    class ScalingClimatology(ClimatologyBase):
        def __init__(self, other_climatology, scale_factor):
            self._other_climatology = other_climatology
            self._scale_factor = scale_factor

        def get_parameter(self, species: str, latitude: float, longitude: float, altitudes: np.ndarray, mjd: float):
            return self._other_climatology.get_parameter(species, latitude, longitude, altitudes, mjd) * self._scale_factor

        def skif_object(self, **kwargs):
            altitudes = np.arange(0, 100000, 500)
            clim_values = dict()

            engine = kwargs['engine']
            reference_point = engine.model_parameters['referencepoint']
            latitude = reference_point[0]
            longitude = reference_point[1]
            mjd = reference_point[3]

            for species in self._other_climatology.supported_species():
                clim_values[species] = self.get_parameter(species, latitude, longitude, altitudes, mjd)

            userdef_clim = ClimatologyUserDefined(altitudes, clim_values)

            return userdef_clim.skif_object()

Here we have taken the reference point from the engine right before the radiative transfer calculation is about to happen,
and then we create our climatology dynamically based upon this.  For this example we have hard coded the altitudes in,
however we could also look at the engines option map and figure out what the optical property table altitudes are.
We now have an object where the logic is implemented inside python, but we can still use it the same as any
`sasktranif` climatology object.

Why Not Just Make a User Defined Climatology?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We could have just as easily make a ClimatologyUserDefined object with the scaled values above and used that instead,
what is better about what we have done?  The primary advantage is that we do not need to know the configuration ahead
of time.  If we wanted to make a UserDefined object from the beginning we would need to know what latitude/longitude/mjd
to use for our values, and if we wanted to change them later we would need to remake the object.  Another advantage
is that we can change things in our climatology relatively easy.  If we wanted to add a method to change the scaling factor
we could do that very easily.  If we wanted to change the scaling factor of the UserDefined climatology we would have
to remake the object, and reset its reference in other places such as the atmosphere class.  Writing these kinds of
objects allows easy integration with the core of the `sasktran` package.