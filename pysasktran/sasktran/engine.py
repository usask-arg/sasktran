import sasktranif.sasktranif as skif
import numpy as np
import sasktran as sk
import logging
from copy import copy
from sasktran.exceptions import wrap_skif_functionfail
from collections import namedtuple
from sasktran.util import DictWithCallback, LowercaseKeyDictWithCallback, to_iter
from sasktran.stokesvector import StokesVector


class Engine(object):
    """
    Base class for objects which can perform radative transfer calculations.
    """

    @wrap_skif_functionfail
    def __init__(self, name: str, geometry: sk.Geometry = None, atmosphere: sk.Atmosphere = None,
                 wavelengths: list = None, options: dict = None):
        """ Base class constructor. The base class cannot be instantiated.
        """
        if name.lower() == 'hr' and not isinstance(self, EngineHR):
            raise ValueError("Using sasktran.Engine('hr') is no longer supported, use sasktran.EngineHR() instead ")

        self._name = name
        self._iskengine = skif.ISKEngine(self._name)
        self.wavelengths = wavelengths
        if geometry is not None:
            self._geometry = geometry
        else:
            self._geometry = sk.Geometry()

        if atmosphere is not None:
            self._atmosphere = atmosphere
        else:
            self._atmosphere = sk.Atmosphere()
        self._needs_reconfigure = True
        self._options = LowercaseKeyDictWithCallback(options) if options is not None else LowercaseKeyDictWithCallback()
        self._options.add_modified_callback(self._set_needs_reconfigure)

    def __getstate__(self):
        state = {k: v for k, v in self.__dict__.items() if k != '_iskengine'}
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._iskengine = skif.ISKEngine(self._name)
        self._needs_reconfigure = True

    def skif_object(self, **kwargs):
        """ Get the underlying ISK object.

        Returns
        -------
        skif.ISKEngine
            Underlying ISKEngine object
        """
        return self._iskengine

    @property
    def options(self):
        """
        dict : Dictionary map of additional engine options.

        Examples
        --------
        >>> import sasktran as sk
        >>> engine = sk.EngineHR()
        >>> engine.options['numordersofscatter'] = 1
        """
        return self._options

    @property
    def atmosphere(self):
        """
        sasktran.Atmosphere : The engines atmosphere definition object.
        """
        return self._atmosphere

    @atmosphere.setter
    def atmosphere(self, value):
        self._atmosphere = value
        self._needs_reconfigure = True

    @property
    def geometry(self):
        """
        sasktran.Geometry : The engines geometry definition object.
        """
        return self._geometry

    @geometry.setter
    def geometry(self, value):
        self._geometry = value
        self._needs_reconfigure = True

    @property
    def wavelengths(self):
        """
        np.ndarray or float : Wavelengths in nm to perform the radiative transfer calculation at.

        Examples
        --------
        >>> import sasktran as sk
        >>> engine = sk.EngineHR()
        >>> engine.wavelengths = 600 # Single wavelength
        >>> engine.wavelengths = [600, 340] # Multiple wavelengths
        """
        return self._wavelengths

    @wavelengths.setter
    def wavelengths(self, value):
        self._wavelengths = np.atleast_1d(value)
        self._needs_reconfigure = True

    @wrap_skif_functionfail
    def calculate_radiance(self, output_format='numpy', full_stokes_vector=False, stokes_orientation='geographic'):
        """ Calculate radiances.

        Parameters
        ----------
        output_format : str, optional
            Determines the output format, one of ['numpy' or 'xarray'].  If output_format='numpy', the radiance
            and optionally the weighting functions are returned as numpy arrays.  If output_format='xarray' the
            radiance and weighting functions are returned inside of an xarray Dataset object.
        full_stokes_vector : bool, optional
            If True, the full stokes vector is output.  If False only the first element of the stokes vector
            is included. Note that setting this does not automatically enable vector calculations
            within the model, it only changes the output format, every model requires a different setting
            to be changed to change the calculation mode.  Default: False
        stokes_orientation : str, optional
            One of 'geographic.  The geographic basis is defined relative to the plane formed by the
            observer position and look vector.  Here, Q, would be the linear polarization component for the nominal
            'up' direction of an instrument.
            Default: 'geographic'.

        Returns
        -------
        np.ndarray **or** np.ndarray, np.ndarray **or** xarray.Dataset
            If output_format='numpy':
            The first np.ndarray is always the computed radiances. The shape of the computed
            radiances is (W, L) where W is the number of wavelengths, and L is the number of
            lines of sight.
            If ``atmosphere.wf_species != None`` then a second np.ndarray will be returned
            containing the computed weighting functions. The shape of this np.ndarray is
            engine specific but typically (W, L, A) where A is the number of altitudes at which
            at which the weighting functions were calculated.
            If full_stokes_vector is set to True, then what is returned is an array of `sasktran.StokesVector`
            objects.

            If output_format='xarray':
            A dataset containing the radiance and weighting functions with appropriate coordinates.
        """
        if self._needs_reconfigure or self._atmosphere.has_changed or self._geometry.has_changed:
            self._iskengine = skif.ISKEngine(self._name)
            self._add_atmosphere_wf()
            self._add_lines_of_sight()
            self._add_options()
            self._initialize_model()
            self._add_atmosphere()
            self._needs_reconfigure = False
            self._geometry._added_to_engine()
            self._atmosphere._added_to_engine()

        if self._wavelengths is not None:
            self._iskengine.SetWavelengths(np.atleast_1d(self._wavelengths))
        else:
            raise ValueError("Property wavelengths needs to be set before calling Engine.calculate_radiance()")

        if full_stokes_vector:
            raw_data = self._iskengine.CalculateStokesVector()[1]

            data = np.zeros(np.shape(raw_data), dtype=StokesVector)
            for idw in range(raw_data.shape[0]):
                for idlos in range(raw_data.shape[1]):
                    data[idw, idlos] = StokesVector.from_skif_object(raw_data[idw, idlos])

        else:
            data = self._iskengine.CalculateRadiance()[1]

        if output_format.lower() == 'numpy':
            return data
        elif output_format.lower() == 'xarray':
            try:
                import xarray as xr
            except ModuleNotFoundError:
                raise ValueError('Output Format "xarray" requires the package xarray to be installed')

            los_vectors = np.vstack([l.look_vector for l in self.geometry.lines_of_sight])
            obs_positions = np.vstack([l.observer for l in self.geometry.lines_of_sight])
            mjds = np.vstack([l.mjd for l in self.geometry.lines_of_sight]).flatten()

            if full_stokes_vector:
                prop_dir = np.vstack([s.propagation_direction for s in data[0, :]])
                theta = np.vstack([s.theta_direction for s in data[0, :]])
                phi = np.vstack([s.phi_direction for s in data[0, :]])

                data = xr.Dataset({'I': (['wavelength', 'los'], np.vectorize(lambda x: x.I)(data)),
                                   'Q': (['wavelength', 'los'], np.vectorize(lambda x: x.Q)(data)),
                                   'U': (['wavelength', 'los'], np.vectorize(lambda x: x.U)(data)),
                                   'V': (['wavelength', 'los'], np.vectorize(lambda x: x.V)(data)),
                                   'propagation_direction': (['los', 'xyz'], prop_dir),
                                   'theta_direction': (['los', 'xyz'], theta),
                                   'phi_direction': (['los', 'xyz'], phi),
                                   'mjd': (['los'], mjds),
                                   'los_vector': (['los', 'xyz'], los_vectors),
                                   'observer_position': (['los', 'xyz'], obs_positions)},
                                  coords={'wavelength': self.wavelengths,
                                          'xyz': ['x', 'y', 'z']})
            else:
                data = xr.Dataset({'radiance': (['wavelength', 'los'], data),
                                   'mjd': (['los'], mjds),
                                   'los_vector': (['los', 'xyz'], los_vectors),
                                   'observer_position': (['los', 'xyz'], obs_positions)},
                                  coords={'wavelength': self.wavelengths,
                                          'xyz': ['x', 'y', 'z']})

            return data
        else:
            raise ValueError('Supported values for output_format are "numpy" and "xarray"')

    @wrap_skif_functionfail
    def _initialize_model(self):
        self._iskengine.InitializeModel()

    def _set_needs_reconfigure(self):
        self._needs_reconfigure = True

    @wrap_skif_functionfail
    def _add_lines_of_sight(self):
        for los in to_iter(self._geometry.lines_of_sight):
            self._iskengine.AddLineOfSight(los.mjd, los.observer, los.look_vector)

        if self._geometry.sun is not None:
            self.options['setsun'] = self._geometry.sun

    @wrap_skif_functionfail
    def _add_atmosphere(self):
        for _, species in self._atmosphere.species.items():
            if species.climatology is not None and species.optical_property is not None:
                self._iskengine.AddSpecies(species.skif_species, species.climatology.skif_object(engine=self, opt_prop=species.optical_property),
                                           species.optical_property.skif_object(engine=self))
        if self._atmosphere.atmospheric_state is not None:
            self._iskengine.SetAtmosphericState(self._atmosphere.atmospheric_state.skif_object())
        if self._atmosphere.brdf is not None:
            self._iskengine.SetBRDF(self._atmosphere.brdf.skif_object(engine=self))

    @wrap_skif_functionfail
    def _add_atmosphere_wf(self):
        if self._atmosphere.wf_species is not None:
            if self._atmosphere.wf_species in self._atmosphere.species:
                op = self._atmosphere.species[self._atmosphere.wf_species].optical_property.skif_object(engine=self)
                self._iskengine.SetProperty('wfspecies', op)
                if 'calcwf' not in [key.lower() for key in self._options]:
                    self._iskengine.SetProperty('calcwf', 2)
            else:
                raise ValueError("WFspecies '%s' not found in species list" % self._atmosphere.wf_species)

    @wrap_skif_functionfail
    def _add_options(self):
        for key, val in self._options.items():
            self._iskengine.SetProperty(key, val)

    @wrap_skif_functionfail
    def _reset_engine(self):
        self._iskengine = skif.ISKEngine(self._name)

        self._needs_reconfigure = True


class EngineHR(Engine):
    """
    Examples
    --------
    >>> import sasktran as sk
    >>> # configure your geometry
    >>> geometry = sk.VerticalImage()
    >>> geometry.from_sza_saa(60, 60, 0, 0, [10, 20, 30, 40], 54372, 0)
    >>> # configure your atmosphere
    >>> atmosphere = sk.Atmosphere()
    >>> atmosphere['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
    >>> atmosphere.atmospheric_state = sk.MSIS90()
    >>> atmosphere['o3'] = sk.Species(sk.O3DBM(), sk.Labow())
    >>> atmosphere['no2'] = sk.Species(sk.NO2Vandaele1998(), sk.Pratmo())
    >>> atmosphere.brdf = sk.Lambertian(0.3)
    >>> # make your engine
    >>> engine = sk.EngineHR(geometry=geometry, atmosphere=atmosphere, wavelengths=[350, 400, 450])
    >>> rad = engine.calculate_radiance()
    >>> print(rad) #doctest: +ELLIPSIS
    [[ 0.11674...  0.10872...  0.04904...  0.01456...]
     [ 0.10837...  0.08537...  0.02897...  0.00804...]
     [ 0.10318...  0.06396...  0.01829...  0.00483...]]
    """
    def __init__(self, geometry: sk.Geometry = None, atmosphere: sk.Atmosphere = None,
                 wavelengths: list = None, options: dict = None):
        super().__init__('hr', geometry, atmosphere, wavelengths, options)

        self._model_parameters = dict()
        self._wf_shape = None

    def _initialize_model(self):
        super()._initialize_model()

        self._model_parameters['referencepoint'] = self._iskengine.GetProperty('ReferencePoint')[1]
        self._model_parameters['sun'] = self._iskengine.GetProperty('sun')[1]

    @wrap_skif_functionfail
    def _add_atmosphere_wf(self):
        if self.atmosphere.wf_species is not None:
            wf_species_handles = []
            for a_species in np.atleast_1d(self.atmosphere.wf_species):
                if not (a_species in self.atmosphere.species):
                    err_msg = "Weighting function species '{}' was not in the species list".format(a_species)
                    # Could be an aerosol special mode
                    special_modes = ['_lognormal_medianradius', '_lognormal_modewidth']
                    for special_mode in special_modes:
                        if a_species.endswith(special_mode):
                            primary_species = a_species.replace(special_mode, '')
                            wf_species_handles.append(self.atmosphere.species[primary_species].species + special_mode)
                            err_msg = None
                    if err_msg is not None:
                        raise ValueError(err_msg)
                else:
                    wf_species_handles.append(self.atmosphere.species[a_species].species)
            self.options['wfspecies'] = " ".join(wf_species_handles)
            if 'calcwf' not in self.options:
                self.options['calcwf'] = 2
            if hasattr(self.atmosphere.wf_species, '__iter__') and not isinstance(self.atmosphere.wf_species, str):
                self._wf_shape = (
                    len(self.wavelengths),
                    len(self.geometry.lines_of_sight),
                    len(self.atmosphere.wf_species),
                    -1
                )
            else:
                self._wf_shape = (
                    len(self.wavelengths),
                    len(self.geometry.lines_of_sight),
                    -1
                )

    def _add_atmosphere(self):
        super()._add_atmosphere()

        # Add in any atmospheric emissions
        for species, emission in self._atmosphere.emissions.items():
            if emission.skif_species_override():
                guid = emission.skif_species_override()
            else:
                guid = species
                skif.AddGeneratedGlobalClimatologyHandleIfNotExists(guid)
            self._iskengine.AddEmission(guid, emission.skif_object(atmosphere=self._atmosphere, engine=self))

    def _add_lines_of_sight(self):
        super()._add_lines_of_sight()

        if self._geometry.reference_point is not None:
            self._iskengine.SetProperty('setreferencepoint', self._geometry.reference_point)

    @property
    def model_parameters(self):
        return self._model_parameters

    @wrap_skif_functionfail
    def calculate_radiance(self, output_format='numpy', full_stokes_vector=False, stokes_orientation='geographic'):
        rad = super().calculate_radiance(output_format=output_format, full_stokes_vector=full_stokes_vector,
                                         stokes_orientation=stokes_orientation)
        if self._atmosphere.wf_species is not None:
            wf = self._iskengine.GetWeightingFunctions()[1]
            wf = wf.reshape(self._wf_shape)
            if output_format.lower() == 'numpy':
                output = namedtuple('EngineOutput', ['radiance', 'weighting_function'])
                return output(radiance=rad, weighting_function=wf)
            elif output_format.lower() == 'xarray':
                if len(self._wf_shape) == 4:
                    if full_stokes_vector:
                        wf = wf.reshape((wf.shape[0], wf.shape[1], wf.shape[2], int(wf.shape[3]/4), 4))
                        for species_idx, species in enumerate(self.atmosphere.wf_species):
                            rad['wf_' + species] = (['wavelength', 'los', 'perturbation', 'stokes'], wf[:, :, species_idx, :, :])
                    else:
                        for species_idx, species in enumerate(self.atmosphere.wf_species):
                            rad['wf_' + species] = (['wavelength', 'los', 'perturbation'], wf[:, :, species_idx, :])
                else:
                    if full_stokes_vector:
                        wf = wf.reshape((wf.shape[0], wf.shape[1], int(wf.shape[2]/4), 4))
                        rad['wf_' + self.atmosphere.wf_species] = (['wavelength', 'los', 'perturbation', 'stokes'], wf)
                    else:
                        rad['wf_' + self.atmosphere.wf_species] = (['wavelength', 'los', 'perturbation'], wf)
                if self.options.get('calcwf', 2) == 2:
                    rad.coords['altitude'] = ('altitude', self.options.get('wfheights', np.arange(500, 99600, 1000)))
                    rad = rad.stack(perturbation=('altitude',))

                brdfwf = self._iskengine.GetProperty('brdfwf')[1].reshape((len(self.wavelengths), -1))
                rad['wf_brdf'] = (['wavelength', 'los'], brdfwf)
                return rad
        else:
            return rad

    @wrap_skif_functionfail
    def cell_optical_depths(self):
        """
        Calculates the optical depths along each line of sight at each wavelength to each cell.

        Returns
        -------
        List[List[Dict]]
            A two-dimensional list of shape (num_wavelength, num_los) where each element is a dictionary containing the keys
            'cell_start_optical_depth' which has the optical depth at the start of the cell and 'cell_start_distance' which
            has the distance from the observer to the start of the cell along the line of sight vector.
        """
        old_los_single_scattering = copy(self.disable_los_single_scattering)
        old_scatter_order = copy(self.num_orders_of_scatter)

        self._options['storeopticaldepth'] = 1
        self.disable_los_single_scattering = True
        self.num_orders_of_scatter = 1

        self.calculate_radiance()

        output = []

        num_wav = len(self._wavelengths)
        num_ray = len(self._geometry.lines_of_sight)

        for wav_idx in range(num_wav):
            wav_output = []
            for ray_idx in range(num_ray):
                cell_start_optical_depth = self._iskengine.GetProperty('loscellstartopticaldepth[{}]'.format(wav_idx * num_ray + ray_idx))[1]
                cell_start_distances = self._iskengine.GetProperty('loscellstartdistance[{}]'.format(wav_idx * num_ray + ray_idx))[1]

                wav_output.append({'cell_start_optical_depth': cell_start_optical_depth,
                                   'cell_start_distance': cell_start_distances})
            output.append(wav_output)

        self._options['storeopticaldepth'] = False
        self.disable_los_single_scattering = old_los_single_scattering
        self.num_orders_of_scatter = old_scatter_order

        return output

    @property
    def num_orders_of_scatter(self):
        """
        Controls the number of orders of scatter used inside the model. The number of orders of scatter controls the
        accuracy of the model.  When set to 1, only light directly scattered from the sun is accounted for, with no
        multiple scattering.

        The model runs substantially faster when set
        to 1 scatter order, therefore a common use case is to do development/testing with this property set to 1.  There
        is very little speed difference between 2 orders of scatter and 50 orders of scatter, so usually it is
        recommended to set this to either 1 or 50.

        **Default:** 50
        """
        return self.options.get('numordersofscatter', 50)

    @property
    def refraction(self):
        """
        Determines whether or not refraction effects should be included in the radiative transfer calculation.
        If set to true, then refraction effects will be considered where the model supports it.  Currently this is
        only available for the observer line of sight rays (solar rays and multiple scattering will not account
        for refraction).

        Note that refraction support is experimental and may have some side effects.  Currently known limitations
        are that the surface_elevation and top_of_atmosphere_altitudeproperty may not work properly with refraction
        enabled.

        Enabling refraction also changes the ray tracer slightly, meaning that calculations with and without refraction
        enabled may differ slightly due to integration accuracies/quadrature points, etc.  It is recommended to
        experiment with the grid_spacing property to ensure that the calculation is done with enough precision.

        **Default:** false
        """
        return bool(self.options.get('userefraction', False))

    @refraction.setter
    def refraction(self, use_refraction: bool):
        self.options['userefraction'] = int(use_refraction)

    @num_orders_of_scatter.setter
    def num_orders_of_scatter(self, num_orders_of_scatter: int):
        if num_orders_of_scatter < 1 or not isinstance(num_orders_of_scatter, int):
            raise ValueError('num_orders_of_scatter must be a positive integer')
        self.options['numordersofscatter'] = num_orders_of_scatter

    @property
    def surface_elevation(self):
        """
        Sets the elevation of the surface in [m] for the model.  This setting can be used to account for changes in the Earth's
        topography.

        **Default:** 0 m
        """
        return self.options.get('surfaceheight', 0)

    @surface_elevation.setter
    def surface_elevation(self, surface_height: float):
        self.options['surfaceheight'] = surface_height

    @property
    def top_of_atmosphere_altitude(self):
        """
        Sets the altitude of the top of the atmosphere.  Above this altitude it is assumed there is no atmosphere.
        The default value of 100000 m is typically suitable for stratospheric applications but for the mesosphere it
        is sometimes necessary to increase this value.

        **Default:** 100000 m
        """
        return self.options.get('toaheight', 100000)

    @top_of_atmosphere_altitude.setter
    def top_of_atmosphere_altitude(self, toa_altitude: float):
        self.options['toaheight'] = toa_altitude

    @property
    def num_threads(self):
        """
        Controls the number of threads to use when multithreading the radiative transfer calculation.  The default value
        of 0 indicates that the number of threads used will be equal to the number of available logical cores on
        the target machine.

        Setting this value to a lower number than the number of available cores can be useful when running a radiative
        transfer calculation in the background and the computer is too slow to multitask.

        **Default:** 0
        """
        return self.options.get('numthreads', 0)

    @num_threads.setter
    def num_threads(self, num_threads: int):
        if num_threads < 0 or not isinstance(num_threads, int):
            raise ValueError('num_threads must be an integer >= 0')
        self.options['numthreads'] = num_threads

    @property
    def grid_spacing(self):
        """
        Sets many options dealing with general grids and discretizations in the model.  In general, if you are interested
        in things on the scale of 1000 m, then set this to 1000 m.  For most applications setting this property is
        sufficient, however for advanced users it may be necessary to set some grids individually.


        **Default:** 1000 m
        """
        return self.options.get('raytracingshells', 1000)

    @grid_spacing.setter
    def grid_spacing(self, grid: float):
        if grid < 0:
            raise ValueError('Grid Spacing should be positive')

        self.options['raytracingshells'] = grid
        self.options['solarraytracingshells'] = grid

        optical_heights = np.arange(self.surface_elevation, self.top_of_atmosphere_altitude, grid / 2)
        optical_heights = np.concatenate((optical_heights, [self.top_of_atmosphere_altitude]))

        self.options['manualopticalheights'] = optical_heights

        diffuse_heights = np.arange(self.surface_elevation + 1e-6, self.top_of_atmosphere_altitude, grid * 5)

        self.options['manualdiffuseheights'] = diffuse_heights

    @property
    def include_emissions(self):
        """
        Enables support for atmospheric emissions inside the model.  If set to True, any :py:class:`Emission <sasktran.Emission>`
        objects that are included in the :py:class:`Atmosphere <sasktran.Atmosphere>` will be used in the radiative
        transfer calculation.  Note that the default value is False, ignoring emissions.

        **Default:** False
        """
        return bool(self.options.get('useemissions', False))

    @include_emissions.setter
    def include_emissions(self, use_emissions: bool):
        self.options['useemissions'] = int(use_emissions)

    @property
    def atmosphere_dimensions(self):
        """
        Sets the number of dimensions for the atmosphere in the radiative transfer calculation.


        ========== ======================================
        Input      Setting
        ========== ======================================
        1          One-dimensional in altitude
        2          Two-dimensional in altitude and angle
                   along the line of sight
        3          Fully three-dimensional
        ========== ======================================

        **Default:** 1
        """
        opt_type = self.options.get('OpticalTableType', 0)

        if opt_type == 0:
            return 1
        if opt_type == 1:
            return 3
        if opt_type == 2:
            return 2

    @atmosphere_dimensions.setter
    def atmosphere_dimensions(self, atmosphere_dim: int):
        if atmosphere_dim < 0 or atmosphere_dim > 3:
            raise ValueError('atmosphere_dim must be either 1, 2, or 3')
        if atmosphere_dim == 1:
            self.options['opticaltabletype'] = 0
        if atmosphere_dim == 2:
            self.options['opticaltabletype'] = 2
        if atmosphere_dim == 3:
            self.options['opticaltabletype'] = 1

    @property
    def polarization(self):
        """
        The polarization mode of the model.  Can either be set to 'scalar' or 'vector'.

        **Default:** 'scalar'
        """
        internal_type = self.options.get('polarizationtype', 0)

        if internal_type == 0:
            return 'scalar'
        else:
            return 'vector'

    @polarization.setter
    def polarization(self, mode):
        if mode.lower() == 'scalar':
            self.options['polarizationtype'] = 0
        elif mode.lower() == 'vector':
            self.options['polarizationtype'] = 99
            self.options['polarizationhigherorderfraction'] = 1e-5

    @property
    def num_diffuse_profiles(self):
        """
        Advanced option that controls the number of diffuse profiles in the model.  If worried about the accuracy
        of the multiple scatter signal, this should be the first option to try changing.  A diffuse profile
        represents a (latitude, longitude) where the multiple scatter signal is calculated.  The default value of 1
        indicates that the multiple scatter signal is only calculated at the reference point, or the average tangent
        point.  If you increase the value to 5, then there will be 5 locations along the line of sight where the full
        multiple scatter calculation is done.  Outside of these areas the multiple scatter signal is interpolated. Odd
        numbers of diffuse profiles are preferred (but not necessary) since this facilitates one profile being placed
        at the tangent point.

        To determine if this setting has an effect on your calculation it is recommended to first compare the results
        between 1, 5, and 11 diffuse profiles.

        The amount of diffuse profiles to obtain an accurate depends heavily on a variety of things, but most heavily
        the solar zenith and solar scattering angles.  When the sun is high in the sky (low solar zenith angle) 1--3
        diffuse profiles is usually sufficient.  However, when looking across the solar terminator (sza ~90, ssa < 90)
        it might be necessary to use upwards of 50 diffuse profiles to obtain a signal accurate to better than 0.5%.
        More information on the number of required diffuse profiles can be found in
        Zawada, D. J., Dueck, S. R., Rieger, L. A., Bourassa, A. E., Lloyd, N. D., and Degenstein, D. A.: High-resolution and Monte Carlo additions to the SASKTRAN radiative transfer model, Atmos. Meas. Tech., 8, 2609-2623, https://doi.org/10.5194/amt-8-2609-2015, 2015.
        """
        return self.options.get('numdiffuseprofilesinplane', 1)

    @num_diffuse_profiles.setter
    def num_diffuse_profiles(self, num_diffuse: int):
        if num_diffuse < 1 or not isinstance(num_diffuse, int):
            raise ValueError('num_diffuse must be an integer greater than 0')

        self.options['numdiffuseprofilesinplane'] = num_diffuse

    @property
    def disable_los_single_scattering(self):
        """
        Experimental option.  If set to true, then single scattering along the line of sight is not calculated.
        Note that single scatter from the ground is still included if the line of sight intersects the Earth.
        Default: False.

        """
        current = self.options.get('disablelosscattering', 0)

        if current == 0:
            return False
        else:
            return True

    @disable_los_single_scattering.setter
    def disable_los_single_scattering(self, disable: bool):
        if disable:
            self._options['disablelosscattering'] = 1
            self._options['opticalphasetableresolution'] = 100
        else:
            self._options['disablelosscattering'] = 0

    @property
    def use_solartable_for_singlescatter(self):
        """
        Experimental option. If set to true, then the solar table will be used for single scattering calculations
        instead of directly tracing solar rays.
        Default: False.
        """
        current = self.options.get('usesolartableforsinglescattering', 0)

        if current == 0:
            return False
        else:
            return True

    @use_solartable_for_singlescatter.setter
    def use_solartable_for_singlescatter(self, use: bool):
        if use:
            self._options['usesolartableforsinglescattering'] = 1
        else:
            self._options['usesolartableforsinglescattering'] = 0

    @property
    def prefill_solartransmission_table(self):
        """
        Experimental option. If set to true, then the solar table will be prefilled instead of filling dynamically.
        Default: False.
        """
        current = self.options.get('prefillsolartransmission', 0)

        if current == 0:
            return False
        else:
            return True

    @prefill_solartransmission_table.setter
    def prefill_solartransmission_table(self, use: bool):
        if use:
            self._options['prefillsolartransmission'] = 1
        else:
            self._options['prefillsolartransmission'] = 0

    def configure_for_cloud(self, grid_spacing_m: float, cloud_altitudes: np.array, cloud_diffuse_spacing_m: float = 100,
                            max_optical_depth_of_cell: float = 0.01, min_extinction_ratio_of_cell: float = 1,
                            max_adaptive_optical_depth: float = 2):
        """
        Configures the HR model for use with a cloud.

        Parameters
        ----------
        grid_spacing_m: float
            Desired grid spacing of the model for all parameters exluding the cloud.
        cloud_altitudes: np.array
            Altitudes that the cloud is non-zero at. [m]
        cloud_diffuse_spacing_m: float, Optional
            Spacing in m that diffuse points should be placed inside the cloud. Lower values will result in a more
            accurate multiple scattering calculatino within the cloud.  Default: 100 m
        max_optical_depth_of_cell: float, Optional
            Maximumum optical depth that is allowed for the adaptive integration step within HR.  Lower values will
            result in a more accurate calculation.  Default: 0.01
        min_extinction_ratio_of_cell: float, Optional
            The minimum ratio of extinction that must be hit for the adaptive integration to split a cell in two, between 0 and 1.
            Higher values will result in a more accurate calculation. It is recommended to leave this at 1 unless you know what you
            are doing.  Default: 1
        max_adaptive_optical_depth: float, Optional
            The maximum optical depth along the ray where splitting where occur.  Default: 2
        """
        if grid_spacing_m < 0:
            raise ValueError('Grid Spacing should be positive')

        optical_heights = np.arange(self.surface_elevation, self.top_of_atmosphere_altitude, grid_spacing_m / 2)
        optical_heights = np.concatenate((optical_heights, [self.top_of_atmosphere_altitude]))

        optical_heights = np.sort(np.unique(np.concatenate((optical_heights, cloud_altitudes))))

        self.options['manualopticalheights'] = optical_heights

        solar_shells = np.arange(self.surface_elevation, self.top_of_atmosphere_altitude, grid_spacing_m)
        solar_shells = np.concatenate((solar_shells, [self.top_of_atmosphere_altitude]))
        solar_shells = np.sort(np.unique(np.concatenate((solar_shells, cloud_altitudes))))

        self.options['manualraytracingshells'] = solar_shells
        self.options['manualsolarraytracingshells'] = solar_shells

        diffuse_heights = np.arange(self.surface_elevation + 1e-6, self.top_of_atmosphere_altitude, grid_spacing_m * 5)
        diff_heights_cloud = np.arange(np.nanmin(cloud_altitudes), np.nanmax(cloud_altitudes), cloud_diffuse_spacing_m)

        diffuse_heights = np.sort(np.unique(np.concatenate((diffuse_heights, diff_heights_cloud))))

        self.options['manualdiffuseheights'] = diffuse_heights

        self.options['maxopticaldepthofcell'] = max_optical_depth_of_cell
        self.options['minextinctionratioofcell'] = min_extinction_ratio_of_cell
        self.options['maxadaptiveopticaldepthofray'] = max_adaptive_optical_depth


class EngineHRSSApprox(EngineHR):
    def __init__(self, geometry: sk.Geometry = None, atmosphere: sk.Atmosphere = None,
                 wavelengths: list = None, options: dict = None,
                 ms_wavelengths: np.array = None):
        """
        Computes the radiance by performing an approximation technique where a single scatter calculation is performed
        at all wavelengths, and then a multiple scatter fraction is done at select wavelengths.  The fraction of
        multiple scatter is then interpolated to approximate the full radiance at all wavelengths.


        Parameters
        ----------
        geometry: sk.Geometry
            Geometry for the calculation
        atmosphere: sk.Atmosphere
            Atmosphere for the calculation
        wavelengths: list
            Wavelengths for the calculation
        options: dict
        ms_wavelengths: np.array
        """
        super().__init__(geometry, atmosphere, wavelengths, options)

        self.ms_wavelengths = ms_wavelengths

    def calculate_radiance(self, output_format='numpy', full_stokes_vector=False, stokes_orientation='geographic'):
        if output_format.lower() != 'xarray':
            raise ValueError('EngineHRSSApprox currently only supports the xarray output format')

        all_wavelengths = copy(self.wavelengths)
        ms_wavelengths = self.ms_wavelengths

        old_scatter_order = copy(self.num_orders_of_scatter)

        # If we are calculating the wf we have to turn it off the SS calculation and reset it for the MS calculation
        old_wf_species = None
        if self.atmosphere.wf_species is not None:
            old_wf_species = copy(self.atmosphere.wf_species)
            self.atmosphere.wf_species = None
            old_wf_mode = self.options.get('calcwf', 0)
            self.options['calcwf'] = 0

        self.num_orders_of_scatter = 1
        ss_rad = super().calculate_radiance(output_format, full_stokes_vector, stokes_orientation)

        self.num_orders_of_scatter = old_scatter_order
        self.wavelengths = ms_wavelengths

        if old_wf_species is not None:
            self.atmosphere.wf_species = old_wf_species
            self.options['calcwf'] = old_wf_mode

        ts_rad = super().calculate_radiance(output_format, full_stokes_vector, stokes_orientation)

        self.wavelengths = all_wavelengths

        ss_frac = ss_rad['radiance'] / ts_rad['radiance']

        ss_frac = ss_frac.interp(wavelength=all_wavelengths)

        ss_rad['radiance'] /= ss_frac

        if old_wf_species is not None:
            ss_rad['wf_' + str(old_wf_species)] = ts_rad['wf_' + str(old_wf_species)].interp(wavelength=all_wavelengths)
            ss_rad['wf_' + str(old_wf_species)] *= (ss_rad['radiance'] / ts_rad['radiance'].interp(wavelength=all_wavelengths))

        return ss_rad

    @property
    def ms_wavelengths(self):
        return self._ms_wavelengths

    @ms_wavelengths.setter
    def ms_wavelengths(self, wavelengths: np.array):
        self._ms_wavelengths = wavelengths


class EngineMC(Engine):
    """
    Examples
    --------
    >>> import sasktran as sk
    >>> # configure your geometry
    >>> geometry = sk.VerticalImage()
    >>> geometry.from_sza_saa(60, 60, 0, 0, [10, 20, 30, 40], 54372, 0)
    >>> # configure your atmosphere
    >>> atmosphere = sk.Atmosphere()
    >>> atmosphere['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
    >>> atmosphere.atmospheric_state = sk.MSIS90()
    >>> atmosphere['o3'] = sk.Species(sk.O3DBM(), sk.Labow())
    >>> atmosphere['no2'] = sk.Species(sk.NO2Vandaele1998(), sk.Pratmo())
    >>> atmosphere.brdf = sk.Lambertian(0.3)
    >>> # configure your options
    >>> options = {'setmcdebugmode': 1234, 'setnumphotonsperlos': 500, 'setsolartabletype': 0}
    >>> # make your engine
    >>> engine = sk.EngineMC(geometry=geometry, atmosphere=atmosphere, wavelengths=[350, 400, 450], options=options)
    >>> rad = engine.calculate_radiance()
    >>> print(rad) #doctest: +ELLIPSIS
    [[ 0.11689...  0.10768...  0.04666...  0.01477...]
     [ 0.10643...  0.08429...  0.02972...  0.00802...]
     [ 0.10589...  0.06424...  0.01797...  0.00499...]]
    """
    def __init__(self, geometry: sk.Geometry = None, atmosphere: sk.Atmosphere = None,
                 wavelengths: list = None, options: dict = None):
        super().__init__('MC', geometry, atmosphere, wavelengths, options)
        self._simultaneous = False
        self._primary_wl = -1.0
        self._opt_wavelengths = 1.0
        self._solar_irr = 1.0
        self._max_raman_orders = []
        self._min_frac_raman = 0.1
        self._amfspecies = None

    def _add_lines_of_sight(self):
        super()._add_lines_of_sight()

        if self._geometry.reference_point is not None:
            self._iskengine.SetProperty('setreferencepoint', self._geometry.reference_point)

    def _add_options(self):
        if self.air_mass_factor == 2 and 'amfspecies' not in self.options:
            if self.air_mass_factor_species not in self.atmosphere.species:
                raise ValueError('Requested air mass factor species was not found in the atmosphere')
            supported = self.atmosphere.species[self.air_mass_factor_species].climatology.supported_species()
            if len(supported) < 1:
                raise ValueError('Requested air mass factor species has no supported species')
            self.options['amfspecies'] = supported[0]

        inelastic = False
        for _, species in self._atmosphere.species.items():
            if species.optical_property._name == 'INELASTICRAYLEIGH':
                inelastic = True

        if inelastic:
            if len(self.max_raman_orders) > 0:
                self.options['maxramanorders'] = self.max_raman_orders
                if self.secondary_output > 0:
                    self.options['scattertype'] = 4
                else:
                    self.options['scattertype'] = 2
                if isinstance(self.min_fraction_higher_raman_order, float):
                    N = len(self.max_raman_orders) + 1
                    self.options['minfractionhigherramanorders'] = [self.min_fraction_higher_raman_order / N] * N
                else:
                    if len(self.max_raman_orders) + 1 != len(self.min_fraction_higher_raman_order):
                        raise ValueError('min_fraction_higher_raman_order must be one element longer than max_raman_orders')
                    self.options['minfractionhigherramanorders'] = self.min_fraction_higher_raman_order
            else:
                if self.secondary_output > 0:
                    self.options['scattertype'] = 3
                else:
                    self.options['scattertype'] = 1

            if isinstance(self.optical_property_wavelengths, float):
                minwl = np.min(self.wavelengths) - 10.0
                maxwl = np.max(self.wavelengths) + 10.0
                dwl = self.optical_property_wavelengths
                self.options['opticalpropertieswavelengths'] = np.arange(minwl, maxwl, dwl)
            else:
                self.options['opticalpropertieswavelengths'] = self.optical_property_wavelengths
        else:
            self.options['scattertype'] = 0
            if self.simultaneous_wavelength:
                self.options['opticalpropertieswavelengths'] = np.sort(self.wavelengths)

        if 'opticalpropertieswavelengths' in self.options:
            if isinstance(self._solar_irr, float):
                self.options['solarspectrum'] = np.ones_like(self.options.get('opticalpropertieswavelengths'))
            else:
                self.options['solarspectrum'] = self.solar_irradiance

        if self.simultaneous_wavelength:
            self.options['wavelengthtype'] = 1
            self.options['radiancewavelengths'] = self.wavelengths
            # default to midpoint
            primary = self.primary_wavelength if self.primary_wavelength > 0 else .5 * (min(self.wavelengths) + max(self.wavelengths))
            # find the closest match
            distance = [abs(wl - primary) for wl in self.wavelengths]
            # tie goes to the longer wavelength (should see further into the atmosphere)
            primary = max([wl for d, wl in zip(distance, self.wavelengths) if d == min(distance)])
            self.options['setprimarywavelength'] = primary
        else:
            self.options['wavelengthtype'] = 0

        super()._add_options()

    def calculate_radiance(self, output_format='numpy', full_stokes_vector=False, stokes_orientation='geographic'):
        data = super().calculate_radiance(output_format, full_stokes_vector, stokes_orientation)
        variance = self.skif_object().GetProperty('variance')[1]
        nwl = len(self.wavelengths)
        nlos = len(self.geometry.lines_of_sight)
        if output_format.lower() == 'numpy':
            output_names = ['radiance', 'radiance_variance']
            output_data = [data, variance.reshape(nwl, nlos)]
            if self.secondary_output > 0:
                sec = np.reshape(self.skif_object().GetProperty('secondarymeasurement')[1], (nwl, nlos))
                secvar = np.reshape(self.skif_object().GetProperty('secondaryvariance')[1], (nwl, nlos))
                name = {1: 'elastic_radiance', 2: 'ring_spectrum', 3: 'filling_in'}[self.secondary_output]
                output_names += [name, f'{name}_variance']
                output_data += [sec, secvar]
            if self.air_mass_factor > 0:
                namf = len(self.air_mass_factor_shells) - 1
                amf = np.reshape(self.skif_object().GetProperty('airmassfactor')[1], (nwl, nlos, namf))
                amfvar = np.reshape(self.skif_object().GetProperty('airmassfactorvariance')[1], (nwl, nlos, namf))
                output_names += ['air_mass_factor', 'air_mass_factor_variance']
                output_data += [amf, amfvar]
            output = namedtuple('EngineOutput', output_names)
            return output(*output_data)
        if output_format.lower() == 'xarray':
            data['radiance_variance'] = (['wavelength', 'los'], variance.reshape(data.dims['wavelength'], data.dims['los']))
            if self.secondary_output > 0:
                sec = np.reshape(self.skif_object().GetProperty('secondarymeasurement')[1], (nwl, nlos))
                secvar = np.reshape(self.skif_object().GetProperty('secondaryvariance')[1], (nwl, nlos))
                name = {1: 'elastic_radiance', 2: 'ring_spectrum', 3: 'filling_in'}[self.secondary_output]
                data[name] = (['wavelength', 'los'], sec)
                data[f'{name}_variance'] = (['wavelength', 'los'], secvar)
            if self.air_mass_factor > 0:
                namf = len(self.air_mass_factor_shells) - 1
                amf = np.reshape(self.skif_object().GetProperty('airmassfactor')[1], (nwl, nlos, namf))
                amfvar = np.reshape(self.skif_object().GetProperty('airmassfactorvariance')[1], (nwl, nlos, namf))
                data['air_mass_factor'] = (['wavelength', 'los', 'amflayer'], amf)
                data['air_mass_factor_variance'] = (['wavelength', 'los', 'amflayer'], amfvar)
                data.coords['lower_shell'] = ('amflayer', self.air_mass_factor_shells[:-1])
                data.coords['upper_shell'] = ('amflayer', self.air_mass_factor_shells[1:])

        return data

    @property
    def num_orders_of_scatter(self):
        """
        Controls the number of orders of scatter used inside the model. The number of orders of scatter controls the
        accuracy of the model.  When set to 1, only light directly scattered from the sun is accounted for, with no
        multiple scattering.

        **Default:** 50
        """
        return self.options.get('numordersofscatter', 50)

    @num_orders_of_scatter.setter
    def num_orders_of_scatter(self, num_orders_of_scatter: int):
        if num_orders_of_scatter < 1 or not isinstance(num_orders_of_scatter, int):
            raise ValueError('num_orders_of_scatter must be a positive integer')
        self.options['numordersofscatter'] = num_orders_of_scatter

    @property
    def max_photons_per_los(self):
        """
        Sets the maximum number of photon paths to simulate for every line of sight. Simulation stops when this maximum
        or the target precision has been reached.

        **Default:** 10000
        """
        return self.options.get('setnumphotonsperlos', 10000)

    @max_photons_per_los.setter
    def max_photons_per_los(self, num_photons: int):
        self.options['setnumphotonsperlos'] = num_photons

    @property
    def target_std(self):
        """
        Sets the target for the relative precision of the simulation. Precision (stdev/result) varies as N^{-1/2}.
        Simulation stops when this target or the maximum number of photon paths has been reached.

        **Default:** 0.01
        """
        return self.options.get('settargetstd', 0.01)

    @target_std.setter
    def target_std(self, target: float):
        self.options['settargetstd'] = target

    @property
    def min_relative_path_weight(self):
        """
        MC simulation of a single photon path stops at some order if each higher-order scatter could add to the
        measured radiance by less than w*{amountAlreadyMeasured}. This can save a lot of time, as only the odd
        "important" high-order path has to be simulated, but it does introduce a small bias. The importance of each
        order of scattering decreases approximately geometrically (except if there are bright clouds, etc), so only
        the first truncated photon contributes significantly to this bias. A value of 0.0 indicates no truncation.

        Choosing  w = 1 / ( 3 * numPhotonsPerLOS^2 )  makes the negative bias approximately one third of the magnitude
        of random noise in the experiment. Very often a photon could make at most a contribution much less than w, so
        the bias is actually much smaller than might be expected.

        **Default:** 0.0
        """
        return self.options.get('setminrelpathweight', 0.0)

    @min_relative_path_weight.setter
    def min_relative_path_weight(self, weight: float):
        self.options['setminrelpathweight'] = weight

    @property
    def surface_elevation(self):
        """
        Sets the elevation of the surface in [m] for the model.  This setting can be used to account for changes in the Earth's
        topography.

        **Default:** 0 m
        """
        return self.options.get('surfaceheight', 0)

    @surface_elevation.setter
    def surface_elevation(self, surface_height: float):
        self.options['surfaceheight'] = surface_height

    @property
    def top_of_atmosphere_altitude(self):
        """
        Sets the altitude of the top of the atmosphere.  Above this altitude it is assumed there is no atmosphere.
        The default value of 100000 m is typically suitable for stratospheric applications but for the mesosphere it
        is sometimes necessary to increase this value.

        **Default:** 100000 m
        """
        return self.options.get('toaheight', 100000)

    @top_of_atmosphere_altitude.setter
    def top_of_atmosphere_altitude(self, toa_altitude: float):
        self.options['toaheight'] = toa_altitude

    @property
    def num_ray_tracing_shells(self):
        """
        Sets the number of evenly spaced atmospheric shells used for ray tracing.

        **Default:** 201
        """
        return self.options.get('setnumraytracingshells', 201)

    @num_ray_tracing_shells.setter
    def num_ray_tracing_shells(self, num_shells: int):
        self.options['setnumraytracingshells'] = num_shells

    @property
    def num_optical_property_altitudes(self):
        """
        Sets the number of evenly spaced altitudes for caching optical properties.

        **Default:** 201
        """
        return self.options.get('setnumoptpropalts', 201)

    @num_optical_property_altitudes.setter
    def num_optical_property_altitudes(self, num_altitudes: int):
        self.options['setnumoptpropalts'] = num_altitudes

    @property
    def solar_table_altitude_delta(self):
        """
        Sets the altitude spacing for the solar transmission table. Only relevant if solar_table_type is 1 or 2.

        **Default:** 500 m
        """
        return self.options.get('setsolartablealtitudedelta', 500.0)

    @solar_table_altitude_delta.setter
    def solar_table_altitude_delta(self, delta: float):
        self.options['setsolartablealtitudedelta'] = delta

    @property
    def solar_table_type(self):
        """
        Sets the method for calculating solar transmissions.

        ========== ======================================
        Input      Setting
        ========== ======================================
        0          No table; solar transmissions
                   calculated on the fly
        1          Two-dimensional in altitude and solar
                   zenith angle
                   along the line of sight
        2          Three-dimensional in altitude, solar
                   zenith angle, and solar longitude
        3          No solar contribution.
        ========== ======================================

        **Default:** 1
        """
        return self.options.get('setsolartabletype', 1)

    @solar_table_type.setter
    def solar_table_type(self, table_type: int):
        if table_type not in [0, 1, 2, 3]:
            raise ValueError('solar_table_type must be 0, 1, 2, or 3')
        self.options['setsolartabletype'] = table_type

    @property
    def ray_tracer_type(self):
        """
        Set the ray tracer type for rays from the sun into the atmosphere, rays for the observer line of sight, and
        rays traced during multiple scattering.

        Solar rays should be straight -- we don't have a good algorithm for tracing curved rays from the sun to a point.
        MS rays should probably be straight -- using curved rays has almost no effect and is a big performance hit.
        Using curved MS rays might be desireable for some people.
        LOS rays can be either straight or curved -- the performance hit for curved rays is pretty small here.

        ========== ======================================
        Input      Setting
        ========== ======================================
        0          Straight ray tracing: straight rays
                   intersecting spherical shells
        1          Curved ray tracing: considering
                   refraction
        2          Generic ray tracing: straight rays
                   intersecting arbitrary shapes
        3          Curve the observer line of sight only
        ========== ======================================

        **Default:** 2
        """
        los = self.options.get('setlosraytracertype', 2)
        ms = self.options.get('setmsraytracertype', 2)

        if los == 0 and ms == 0:
            return 0
        elif los == 1 and ms == 1:
            return 1
        elif los == 2 and ms == 2:
            return 2
        elif los == 1 and ms == 0:
            return 3
        else:
            raise ValueError('unknown ray_tracer_type')

    @ray_tracer_type.setter
    def ray_tracer_type(self, tracer_type: int):
        if tracer_type not in [0, 1, 2, 3]:
            raise ValueError('ray_tracer_type must be 0, 1, or 2, 3')
        self.options['setlosraytracertype'] = tracer_type
        self.options['setmsraytracertype'] = tracer_type
        self.options['setsolarraytracertype'] = tracer_type
        if tracer_type == 1:
            self.options['setsolarraytracertype'] = 2  # curved solar transmission not recommended
        if tracer_type == 3:
            self.options['setlosraytracertype'] = 1
            self.options['setmsraytracertype'] = 0
            self.options['setsolarraytracertype'] = 0

    @property
    def optical_property_integration_type(self):
        """
        Set the optical property table integrator (quadrature) type. Adaptive will split the integration if the
        optical depth exceeds max_optical_depth_of_cell.

        ========== ======================================
        Input      Setting
        ========== ======================================
        0          Straight
        1          Adaptive
        ========== ======================================

        **Default:** 0
        """
        return self.options.get('setoptpropinttype', 0)

    @optical_property_integration_type.setter
    def optical_property_integration_type(self, int_type: int):
        if int_type not in [0, 1]:
            raise ValueError('optical_property_integration_type must be 0 or 1')
        self.options['setoptpropinttype'] = int_type

    @property
    def max_optical_depth_of_cell(self):
        """
        Set the maximum optical depth of a cell before the adaptive integrator will split it. Only relevant if
        optical_property_integration_type is 1.

        **Default:** 0.1
        """
        return self.options.get('maxopticaldepthofcell', 0.1)

    @max_optical_depth_of_cell.setter
    def max_optical_depth_of_cell(self, optical_depth: float):
        if optical_depth <= 0:
            raise ValueError('max_optical_depth_of_cell must be positive')
        self.options['maxopticaldepthofcell'] = optical_depth

    @property
    def optical_table_type(self):
        """
        Sets the type of table used to hold the optical properties of the atmosphere.

        ========== ======================================
        Input      Setting
        ========== ======================================
        0          One-dimensional (altitude)
        1          Three-dimensional (Delaunay sphere)
        ========== ======================================

        **Default:** 0
        """
        return self.options.get('setopticaltabletype', 0)

    @optical_table_type.setter
    def optical_table_type(self, table_type: int):
        self.options['setopticaltabletype'] = table_type

    @property
    def debug_mode(self):
        """
        A non-zero value will fix the random number generator seed to that value and disable multi-threading.

        **Default:** 0
        """
        return self.options.get('setmcdebugmode', 0)

    @debug_mode.setter
    def debug_mode(self, mode: int):
        self.options['setmcdebugmode'] = mode

    @property
    def scatter_position_resolution(self):
        """
        Sets the distance to determine the simulated scatter positions within.

        **Default:** 50 m
        """
        return self.options.get('setscatterpositionres', 50.0)

    @scatter_position_resolution.setter
    def scatter_position_resolution(self, distance: float):
        self.options['setscatterpositionres'] = distance

    @property
    def min_fraction_higher_order(self):
        """
        Sets the fraction of samples that will automatically propagate to the highest order of scatter, regardless of
        contribution to the variance. Remaining samples will be truncated to minimize the variance with the minimum
        number of simulated scatter events. Set to 1.0 to disable this optimization.

        **Default:** 0.1
        """
        return self.options.get('setminfractionhigherorder', 0.1)

    @min_fraction_higher_order.setter
    def min_fraction_higher_order(self, fraction: float):
        if fraction < 0.0 or fraction > 1.0:
            raise ValueError('min_fraction_higher_order must be between 0.0 and 1.0 inclusive')
        self.options['setminfractionhigherorder'] = fraction

    @property
    def simultaneous_wavelength(self):
        """
        If true, multiple wavelengths will be calculated simultaneously using common ray tracing. This should not be
        used if air mass factors are being calculated.

        **Default:** False
        """
        return self._simultaneous

    @simultaneous_wavelength.setter
    def simultaneous_wavelength(self, simultaneous: bool):
        self._simultaneous = simultaneous
        self._set_needs_reconfigure()

    @property
    def primary_wavelength(self):
        """
        Sets the wavelength that will be used to guide the ray tracing when simultaneous_wavelength is True. If the
        wavelength set here does not belong to self.wavelengths, the nearest match will be used. If the wavelength set
        here is not positive, the wavelength nearest to the midpoint of the range will be used.

        **Default:** -1.0
        """
        return self._primary_wl

    @primary_wavelength.setter
    def primary_wavelength(self, wavelength: float):
        self._primary_wl = wavelength
        self._set_needs_reconfigure()

    @property
    def air_mass_factor(self):
        """
        Sets the type of air mass factor calculation.

        ========== ======================================
        Input      Setting
        ========== ======================================
        0          No air mass factor calculation
        1          Air mass factors calculated via path
                   length
        2          Air mass factors calculated via path
                   optical depth
        ========== ======================================

        **Default:** 0
        """
        amf = self.options.get('secondaryoutput', 0)
        if amf == 1:
            return 1
        elif amf == 2:
            return 2
        else:
            return 0

    @air_mass_factor.setter
    def air_mass_factor(self, amf_type: int):
        if amf_type not in [0, 1, 2]:
            raise ValueError('air_mass_factor must be 0, 1, or 2')
        self.options['secondaryoutput'] = amf_type

    @property
    def air_mass_factor_species(self):
        """
        Sets the species used for air mass factor optical depth calculations.

        The first supported species listed by the climatology in the atmosphere object associated with this name will
        be used. If a supported species further down the list is required, set the 'amfspecies' option manually with
        a handle from sk.handles.standard_handles().

        Only relevant if air_mass_factor is set to 2.

        **Default:** None
        """
        return self._amfspecies

    @air_mass_factor_species.setter
    def air_mass_factor_species(self, amf_species: str):
        self._amfspecies = amf_species
        self._set_needs_reconfigure()

    @property
    def air_mass_factor_shells(self):
        """
        Sets the altitudes that define the layers to calculate air mass factors in. Only relevant if air_mass_factor is
        not 0.

        **Default:** 500 m intervals from 0 m to 100 000 m
        """
        return self.options.get('manualamfshells', np.linspace(0, 1e5, 101))

    @air_mass_factor_shells.setter
    def air_mass_factor_shells(self, shells: np.ndarray):
        self.options['manualamfshells'] = shells

    @property
    def secondary_output(self):
        """
        Select secondary quantities to return alongside inelastic radiance. Only relevant if an inelastic species is
        included in the atmosphere.

        ========== ======================================
        Input      Setting
        ========== ======================================
        0          No secondary output
        1          Elastic radiance
        2          Ring spectrum
        3          Filling-in parameter
        ========== ======================================

        **Default:** 0
        """
        secondaryoutput = self.options.get('secondaryoutput', 0)
        if secondaryoutput > 2:
            secondaryoutput -= 2
        else:
            secondaryoutput = 0
        return secondaryoutput

    @secondary_output.setter
    def secondary_output(self, output: int):
        if output not in [0, 1, 2, 3]:
            raise ValueError('inealstic_radiance must be 0, 1, 2, 3, or 4')
        if output > 0:
            output += 2
        self.options['secondaryoutput'] = output

    @property
    def max_raman_orders(self):
        """
        Defines the maximum order of Raman scatter (equal to the length of the list), as well as the maximum
        total order of scatter for each order of Raman scatter (the first element corresponding to the first Raman
        order, and so on).

        Setting this turns the optimized inelastic mode on; if it is left as an empty list, Raman scattering events
        will occur at the same rate they would physically. In optimized inelastic mode, Raman scattering events are
        forced in order to minimize the variance of the radiance, or the variance of the secondary output if
        secondary_output is not 0.

        Only relevant if an inelastic species is included in the atmosphere.

        **Default:** []
        """
        return self._max_raman_orders

    @max_raman_orders.setter
    def max_raman_orders(self, orders: list):
        if isinstance(orders, int):
            orders = [orders]
        for i, o in enumerate(orders):
            if o < i:
                raise ValueError('no element of max_raman_orders can be smaller than its corresponding Raman order')
        self._max_raman_orders = orders
        self._set_needs_reconfigure()

    @property
    def min_fraction_higher_raman_order(self):
        """
        Defines the fraction of all samples that will, for each order of Raman scatter, always be simulated up to the
        maximum total order defined by max_raman_orders. The length of min_fraction_higher_order must be one greater
        than the length of max_raman_order, with the extra element at the front of the list represent the 0th Raman
        order.

        Another option is to specify a scalar here, which will be divided evenly among the Raman orders (this will
        likely not be optimal).

        The purpose of these minimum fractions is to reduce the small bias introduced by using variance knowledge
        to determine what samples will be taken next.

        Only relevant if an inelastic species is included in the atmosphere and max_raman_orders has length
        greater than 0.

        **Default:** 0.1
        """
        return self._min_frac_raman

    @min_fraction_higher_raman_order.setter
    def min_fraction_higher_raman_order(self, fractions: list):
        if not isinstance(fractions, float):
            for f in fractions:
                if f < 0.0:
                    raise ValueError('no element of min_fraction_higher_raman_order can be negative')
            if sum(fractions) > 1.0:
                raise ValueError('the sum of min_fraction_higher_raman_order must not be greater than 1.0')
        self._min_frac_raman = fractions
        self._set_needs_reconfigure()

    @property
    def optical_property_wavelengths(self):
        """
        The wavelengths where optical properties will be cached. This is only relevant if inelastic radiance is being
        calculated. If inelastic scatters take a photon outside the range defined here, a warning will be thrown and
        optical properties will be extrapolated. If a single number is set here, it will be interpreted as a resolution
        and the wavelengths will be automatically set over a suitable range at this resolution.

        **Default:** 1.0
        """
        return self._opt_wavelengths

    @optical_property_wavelengths.setter
    def optical_property_wavelengths(self, wavelengths_or_resolution: np.ndarray):
        self._opt_wavelengths = wavelengths_or_resolution
        self._set_needs_reconfigure()

    @property
    def solar_irradiance(self):
        """
        The solar spectrum used for inelastic radiance calculation, on the wavelength grid specified by
        optical_property_wavelengths. If a scalar is given, the solar spectrum will be constant.

        **Default:** 1.0
        """
        return self._solar_irr

    @solar_irradiance.setter
    def solar_irradiance(self, irradiance: np.ndarray):
        self._solar_irr = irradiance
        self._set_needs_reconfigure()