import sasktranif.sasktranif as skif
import numpy as np
import sasktran as sk
from sasktran.exceptions import wrap_skif_functionfail
from collections import namedtuple
from sasktran.util import to_iter
from sasktran.engine import Engine


class EngineTIR(Engine):
    """
    Examples
    --------
    >>> import sasktran as sk
    >>> import numpy as np
    >>> from sasktran.tir.engine import EngineTIR
    >>> from sasktran.tir.opticalproperty import HITRANChemicalTIR
    >>> # Select wavenumbers to do calculation at
    >>> wavenum = np.linspace(775.0, 775.003, 4)
    >>> # Wavelengths passed to engine must in increasing order
    >>> wavelen = 1.0e7 / wavenum[::-1]
    >>> # Create limb view geometry for a 20 km tangent altitude observed from a 40 km balloon altitude
    >>> geometry = sk.VerticalImage()
    >>> geometry.from_sza_saa(60, 60, 0, 0, [10, 20, 30, 40], 54372, 0)
    >>> # Create an atmosphere with absorbing/emitting ozone
    >>> atmosphere = sk.Atmosphere()
    >>> o3_opt = HITRANChemicalTIR('O3', micro_window_margin=50)
    >>> atmosphere['ozone'] = sk.Species(o3_opt, sk.Labow())
    >>> # Create the engine and perform the calculation
    >>> engine = EngineTIR(geometry=geometry, atmosphere=atmosphere, wavelengths=wavelen)
    >>> rad = engine.calculate_radiance()
    >>> print(rad)
    [[8.73341577e+12 1.09506519e+13 1.32697534e+13 4.41957608e+12]
     [1.47059359e+13 1.67048571e+13 2.02328109e+13 1.51375632e+13]
     [1.06787128e+13 1.28978370e+13 1.59440769e+13 7.41160317e+12]
     [5.22908338e+12 7.13817250e+12 7.17447283e+12 9.33290088e+11]]
    """

    def __init__(self, geometry: sk.Geometry = None, atmosphere: sk.Atmosphere = None,
                 wavelengths: list = None, options: dict = None):
        super().__init__('tir', geometry, atmosphere, wavelengths, options)
        self._model_parameters = dict()
        self._wf_shape = None
        self._can_use_cached_cross_sections = False
        self._use_cached_cross_sections = False
        self._cached_species = []
        self._prev_wf_species = []
        self._wf_species_and_temperature = None

    def _initialize_model(self):
        super()._initialize_model()

        self._model_parameters['referencepoint'] = self._iskengine.GetProperty('ReferencePoint')[1]

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

    @wrap_skif_functionfail
    def _add_lines_of_sight(self):
        # overridden to ignore the sun
        for los in to_iter(self._geometry.lines_of_sight):
            self._iskengine.AddLineOfSight(los.mjd, los.observer, los.look_vector)
        self._can_use_cached_cross_sections = False

    @wrap_skif_functionfail
    def _add_atmosphere_wf(self):
        self._wf_species_and_temperature = self.atmosphere.wf_species
        if self.do_temperature_wf:
            if self.atmosphere.wf_species is None:
                self._wf_species_and_temperature = 'temperature'
            else:
                self._wf_species_and_temperature = np.append(self._wf_species_and_temperature, 'temperature')

        if self._wf_species_and_temperature is not None:
            wf_species_handles = []
            for a_species in np.atleast_1d(self._wf_species_and_temperature):
                if a_species == 'temperature':
                    continue
                if not (a_species in self.atmosphere.species):
                    err_msg = "Weighting function species '{}' was not in the species list".format(a_species)
                    raise ValueError(err_msg)
                wf_species_handles.append(self.atmosphere.species[a_species].species)
            if len(wf_species_handles) > 0:  # will be 0 length when only temperature wf is computed
                self.options['wfspecies'] = " ".join(wf_species_handles)
            if not ('wfheights' in self.options):
                self.options['wfheights'] = np.linspace(500, 99500, 100)
            if not ('wfwidths' in self.options):
                self.options['wfwidths'] = np.array([1000] * 100)
            if hasattr(self._wf_species_and_temperature, '__iter__') and not isinstance(
                    self._wf_species_and_temperature, str):
                self._wf_shape = (
                    len(self.wavelengths),
                    len(self.geometry.lines_of_sight),
                    len(self._wf_species_and_temperature),
                    len(self.options['wfheights'])
                )
            else:
                self._wf_shape = (
                    len(self.wavelengths),
                    len(self.geometry.lines_of_sight),
                    len(self.options['wfheights'])
                )

        for a_species in np.atleast_1d(self.atmosphere.wf_species):
            if not (a_species in self._prev_wf_species):
                # weighting function species have changed so can't use cached cross sections
                self._can_use_cached_cross_sections = False

    @property
    def model_parameters(self):
        return self._model_parameters

    @Engine.wavelengths.setter
    def wavelengths(self, value):
        self._wavelengths = np.atleast_1d(value)
        self._can_use_cached_cross_sections = False

    @wrap_skif_functionfail
    def calculate_radiance(self, output_format='numpy', full_stokes_vector=False):
        if full_stokes_vector:
            raise ValueError('EngineTIR does not support polarized calculations')

        if set(self._cached_species) != set(self.atmosphere.species.keys()):
            self._can_use_cached_cross_sections = False

        if self._use_cached_cross_sections and self._can_use_cached_cross_sections:
            return self._calculate_radiance_from_cache(output_format)
        self._iskengine.SetProperty('usecache', False)
        self.options['wavelengths'] = self.wavelengths
        rad = super().calculate_radiance(output_format)
        self._can_use_cached_cross_sections = True
        self._cached_species = self.atmosphere.species.keys()
        self._prev_wf_species = self.atmosphere.wf_species
        if self._prev_wf_species is None:
            self._prev_wf_species = []
        if self._wf_species_and_temperature is not None:
            wf = self._iskengine.GetWeightingFunctions()[1].reshape(self._wf_shape)
            if output_format.lower() == 'numpy':
                output = namedtuple('EngineOutput', ['radiance', 'weighting_function'])
                return output(radiance=rad, weighting_function=wf)
            elif output_format.lower() == 'xarray':
                if hasattr(self._wf_species_and_temperature, '__iter__') and not isinstance(
                        self._wf_species_and_temperature, str):
                    for idx, species in enumerate(self._wf_species_and_temperature):
                        rad['wf_' + species] = (['wavelength', 'los', 'perturbation'], wf[:, :, idx, :])
                else:
                    rad['wf_' + self._wf_species_and_temperature] = (['wavelength', 'los', 'perturbation'], wf)
                rad.coords['altitude'] = ('altitude', self.options.get('wfheights', np.linspace(500, 99500, 100)))
                rad = rad.stack(perturbation=('altitude',))
                return rad
        else:
            return rad

    @wrap_skif_functionfail
    def _calculate_radiance_from_cache(self, output_format='numpy'):
        self._iskengine.SetProperty('usecache', True)
        if self._atmosphere.has_changed:
            self._add_atmosphere_wf()
            self._add_atmosphere()
            self._needs_reconfigure = False
            self._geometry._added_to_engine()
            self._atmosphere._added_to_engine()
        data = self._iskengine.CalculateRadiance()[1]

        if output_format.lower() == 'numpy':
            rad = data
        elif output_format.lower() == 'xarray':
            try:
                import xarray as xr
            except ModuleNotFoundError:
                raise ValueError('Output Format "xarray" requires the package xarray to be installed')

            los_vectors = np.vstack([l.look_vector for l in self.geometry.lines_of_sight])
            obs_positions = np.vstack([l.observer for l in self.geometry.lines_of_sight])
            mjds = np.vstack([l.mjd for l in self.geometry.lines_of_sight]).flatten()

            rad = xr.Dataset({'radiance': (['wavelength', 'los'], data),
                              'mjd': (['los'], mjds),
                              'los_vector': (['los', 'xyz'], los_vectors),
                              'observer_position': (['los', 'xyz'], obs_positions)},
                             coords={'wavelength': self.wavelengths,
                                     'xyz': ['x', 'y', 'z']})
        else:
            raise ValueError('Supported values for output_format are "numpy" and "xarray"')

        if self._wf_species_and_temperature is not None:
            wf = self._iskengine.GetWeightingFunctions()[1].reshape(self._wf_shape)
            if output_format.lower() == 'numpy':
                output = namedtuple('EngineOutput', ['radiance', 'weighting_function'])
                return output(radiance=rad, weighting_function=wf)
            elif output_format.lower() == 'xarray':
                if hasattr(self._wf_species_and_temperature, '__iter__') and not isinstance(
                        self._wf_species_and_temperature, str):
                    for idx, species in enumerate(self._wf_species_and_temperature):
                        rad['wf_' + species] = (['wavelength', 'los', 'perturbation'], wf[:, :, idx, :])
                else:
                    rad['wf_' + self._wf_species_and_temperature] = (['wavelength', 'los', 'perturbation'], wf)
                rad.coords['altitude'] = ('altitude', self.options.get('wfheights', np.linspace(500, 99500, 100)))
                rad = rad.stack(perturbation=('altitude',))
                return rad
        else:
            return rad

    @property
    def refraction(self):
        """
        If set to true, then refractive ray tracing is performed for the observer line of sight rays.

        **Default:** false
        """
        return bool(self.options.get('userefraction', False))

    @refraction.setter
    def refraction(self, use_refraction: bool):
        self._can_use_cached_cross_sections = False
        self.options['userefraction'] = int(use_refraction)

    @property
    def surface_elevation(self):
        """
        Sets the surface elevation in meters.

        **Default:** 0 m
        """
        return self.options.get('surfaceheight', 0)

    @surface_elevation.setter
    def surface_elevation(self, surface_height: float):
        self._can_use_cached_cross_sections = False
        self.options['surfaceheight'] = surface_height

    @property
    def top_of_atmosphere_altitude(self):
        """
        Sets the altitude of the top of the atmosphere. Above this altitude it is assumed there is no atmosphere.
        Due to the assumption of local thermodynamic equilibrium within the Thermal InfraRed engine, increasing
        the value above 100000 m is not recommended because this simplification can not be used at higher altitudes.

        **Default:** 100000 m
        """
        return self.options.get('toaheight', 100000)

    @top_of_atmosphere_altitude.setter
    def top_of_atmosphere_altitude(self, toa_altitude: float):
        self._can_use_cached_cross_sections = False
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

    # TODO: add optical properties and emission table spacing setters
    @property
    def grid_spacing(self):
        """
        Sets the vertical spacing of homogeneous atmosphere layers in [m].

        **Default:** 1000 m
        """
        return self.options.get('raytracingshells', 1000)

    @grid_spacing.setter
    def grid_spacing(self, grid: float):
        self._can_use_cached_cross_sections = False
        if grid < 0:
            raise ValueError('Grid Spacing should be positive')

        self.options['raytracingshells'] = grid
        self.options['opticalpropertiesheightres'] = grid

    @property
    def atmosphere_dimensions(self):
        """
        Sets the number of dimensions for the atmosphere in the radiative transfer calculation.


        ========== ======================================
        Input      Setting
        ========== ======================================
        1          One-dimensional in altitude
        2          Two-dimensional in altitude and along
                   along the line of sight
        ========== ======================================

        **Default:** 1
        """
        return self.options.get('opticaltabledimensions', 1)

    @atmosphere_dimensions.setter
    def atmosphere_dimensions(self, atmosphere_dim: int):
        self._can_use_cached_cross_sections = False
        if atmosphere_dim < 0 or atmosphere_dim > 2:
            raise ValueError('atmosphere_dim must be either 1 or 2')
        self.options['opticaltabledimensions'] = atmosphere_dim

    @property
    def adaptive_integration(self):
        """
        A boolean property that, if True, enables dynamic cell splitting during the radiative transfer integration.
        During the ray tracing process, the line of sight is split into cells based on the intersection of the line of
        sight with geometric shapes used to split the atmosphere into homogeneous layers. When the optical depth of a
        cell is computed, the cell is split if one or both of the following conditions are satisfied:

            1. The optical depth of the cell exceeds the value specified by the property ``max_optical_depth_of_cell``
            2. The ratio of the smaller of the two values of extinction at the endpoints of the cell to the larger
               value is less than the property ``min_extinction_ratio_of_cell``

        **Default:** True
        """
        return self.options.get('useadaptiveintegration', True)

    @adaptive_integration.setter
    def adaptive_integration(self, value):
        self.options['useadaptiveintegration'] = value

    @property
    def max_optical_depth_of_cell(self):
        """
        Set the maximum optical depth of a single cell is allowed to have before it is split into two cells, when
        adaptive integration is enabled. If the property ``adaptive_integration`` is set to False, this property has no
        effect.

        **Default:** 0.1
        """
        return self.options.get('maxopticaldepthofcell', 0.1)

    @max_optical_depth_of_cell.setter
    def max_optical_depth_of_cell(self, value):
        self.options['maxopticaldepthofcell'] = value

    @property
    def min_extinction_ratio_of_cell(self):
        """
        Set the minimum ratio between the extinction at the endpoints of single cell that is allowed before the cell
        is split, when adaptive integration is enabled. If the property ``adaptive_integration`` is set to False, this
        property has no effect. The ratio is calculated as the smaller of the two extinction values divided by the
        larger.

        **Default** 0.9
        """
        return self.options.get('minextinctionratioofcell', 0.9)

    @min_extinction_ratio_of_cell.setter
    def min_extinction_ratio_of_cell(self, value):
        self.options['minextinctionratioofcell'] = value

    @property
    def linear_extinction(self):
        """
        If this property is set to True, the engine calculates the optical depth of a cell by allowing the extinction
        to vary linearly with height, between the values at the endpoints of the cell. If this property is set to False,
        the extinction of a cell is constant, calculated as the average of the values at the endpoints.

        **Default:** True
        """
        return self.options.get('uselinearextinction', True)

    @linear_extinction.setter
    def linear_extinction(self, value):
        self.options['uselinearextinction'] = value

    @property
    def ground_emissivity(self):
        """
        Sets the scalar emissivity of the surface of the Earth. If the line of sight ends at the ground, the emission
        from the Earth's surface must be taken into account. In the TIR engine this emission is modelled as a
        black body with constant emissivity. For limb-viewing geometries this setting will have no effect as there are
        no scattering effects in the TIR engine.

        **Default:** 1
        """
        return self.options.get('groundemissivity', 1)

    @ground_emissivity.setter
    def ground_emissivity(self, value):
        self.options['groundemissivity'] = value

    @property
    def wf_heights(self):
        """
        Heights, in meters, to compute analytic weighting functions at.

        **Default:** ``np.linspace(500, 99500, 100)``  # 500, 1500, 2500, ..., 99500
        """
        return self.options.get('wfheights', np.linspace(500, 99500, 100))

    @wf_heights.setter
    def wf_heights(self, value):
        self.options['wfheights'] = value

    @property
    def wf_widths(self):
        """
        Weighting function widths. The weighting functions determine the effect of perturbing a species population at
        each height specified by ``wf_heights`` on the radiance. The property wf_widths species the vertical distance
        from each value in ``wf_heights`` to the height where the perturbation is zero. Therefore, ``wf_heights`` and ``wf_widths``
        must have identical lengths.

        i.e. for a weighting function calculation at ``wf_heights[i]``, the perturbation is maximum at ``wf_heights[i]`` and
        decreases linearly (with height) to 0 at ``(wf_heights[i] - wf_widths[i])`` and ``(wf_heights[i] + wf_widths[i])``.

        **Default:** ``np.array([1000] * 100)``  # sets a width of 1000 m at every height in wf_heights
        """
        return self.options.get('wfwidths', np.array([1000] * 100))

    @wf_widths.setter
    def wf_widths(self, value):
        self.options['wfwidths'] = value

    @property
    def use_cached_cross_sections(self):
        """
        Designed to enable faster iterative radiance calculations, this setting instructs the TIR engine to compute
        cross sections using cached results. This can only be done if:

            1. ``calculate_radiance()`` has been called at least once
            2. The only change to the engine object since the last call to ``calculate_radiance()`` is the number density
               climatology of one or more species in the atmosphere property. The optical properties used for these
               species must the same; if a new optical property is specified, the resulting radiance will be as though
               the old optical property were used.
            3. The species whose climatologies have been changed must also be set as the ``wf_species`` of the atmosphere
               property, i.e. only the species for which weighting functions are computed may have their climatologies
               changed

        **Default:** True
        """
        return self._use_cached_cross_sections

    @use_cached_cross_sections.setter
    def use_cached_cross_sections(self, value):
        self._use_cached_cross_sections = value

    @property
    def do_temperature_wf(self):
        """
        Set this property to True to enable temperature weighting function calculations.

        **Default:** False
        """
        return self.options.get('calctemperaturewf', False)

    @do_temperature_wf.setter
    def do_temperature_wf(self, value):
        self.options['calctemperaturewf'] = value

    @property
    def do_vmr_wf(self):
        """
        Controls the units of gas species weighting functions. The default setting configures weighting functions
        to be given in (radiance / number density) where number density has units of cm^-3.
        Setting this property to True tells the engine to calculate weighting functions with units of
        (radiance / VMR).

        **Default:** False
        """
        return self.options.get('usevmrwfunit', False)

    @do_vmr_wf.setter
    def do_vmr_wf(self, value):
        self.options['usevmrwfunit'] = value
