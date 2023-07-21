import numpy as np
import sasktran as sk
import logging
from typing import Union
from sasktran.exceptions import wrap_skif_functionfail
from collections import namedtuple
from sasktran.util import LowercaseKeyDictWithCallback, to_iter


class EngineDO(sk.Engine):
    """
    Examples
    --------
    >>> import sasktran as sk
    >>> # Configure the geometry
    >>> geometry = sk.NadirGeometry()
    >>> tempo = sk.Geodetic()
    >>> tempo.from_lat_lon_alt(0, -100, 35786000)
    >>> lats = np.linspace(20, 50, 4)
    >>> lons = np.linspace(-110, -90, 3)
    >>> lats, lons = np.meshgrid(lats, lons, indexing='ij')
    >>> geometry.from_lat_lon(mjd=57906.63472, observer=tempo, lats=lats, lons=lons)
    >>> # Configure the atmosphere
    >>> atmosphere = sk.ConstituentAtmosphere()
    >>> atmosphere['rayleigh'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
    >>> atmosphere['o3'] = sk.Species(sk.O3DBM(), sk.Labow())
    >>> atmosphere['no2'] = sk.Species(sk.NO2Vandaele1998(), sk.Pratmo())
    >>> atmosphere.brdf = sk.Kokhanovsky()
    >>> # Calcualte radiance
    >>> engine = sk.EngineDO(geometry=geometry, atmosphere=atmosphere, wavelengths=[350, 450])
    >>> rad = engine.calculate_radiance('numpy')
    >>>
    >>> # print radiances
    >>> print(rad[0].reshape((4, 3))) # radiance at 350nm
    [[ 0.15726175  0.19468353  0.22767423]
     [ 0.1666762   0.20027906  0.22962238]
     [ 0.17116222  0.1995542   0.22411801]
     [ 0.17099479  0.1932349   0.21235274]]
    >>> print(rad[1].reshape((4, 3))) # radiance at 450nm
    [[ 0.15884236  0.19815234  0.23287579]
     [ 0.16822132  0.20341264  0.23415734]
     [ 0.17226753  0.20184478  0.22736965]
     [ 0.1716201   0.19441591  0.21402527]]
    """
    def __init__(self, geometry: sk.Geometry = None, atmosphere: sk.Atmosphere = None,
                 wavelengths: list = None, options: dict = None):
        super().__init__('do', geometry, atmosphere, wavelengths, options)
        self._model_parameters = dict()
        self._wf_shape = None
        self._options = LowercaseKeyDictWithCallback(options) if options is not None else LowercaseKeyDictWithCallback()
        self._options.add_modified_callback(self._set_needs_reconfigure)
        self._num_stokes = 1
        self._options['averagereferencepoint'] = self.options.get('averagereferencepoint', 1)

    def _initialize_model(self):
        super()._initialize_model()
        self._iskengine.SetPolarizationMode(self._num_stokes)
        self._model_parameters['sun'] = self._iskengine.GetProperty('sun')[1]
        self.model_parameters['numstreams'] = int(self._iskengine.GetProperty('numstreams')[1])
        self.model_parameters['numlayers'] = int(self._iskengine.GetProperty('numlayers')[1])
        self.model_parameters['numbrdfexpansionterms'] = int(self._iskengine.GetProperty('numbrdfexpansionterms')[1])

    def _wf_altitudes(self):
        # Logic to determine the weighting function altitudes
        if 'wfaltitudes' in self.options:
            return self.options['wfaltitudes']
        return self.options.get('altitudegrid', np.arange(0, 100001, 500))

    @wrap_skif_functionfail
    def _add_atmosphere_wf(self):
        if self.atmosphere.wf_species is not None:
            wf_species_handles = []
            for a_species in np.atleast_1d(self.atmosphere.wf_species):
                # We have special modes for brdf
                if a_species == 'brdf':
                    wf_species_handles.append(a_species)
                else:
                    if not (a_species in self.atmosphere.species):
                        err_msg = "Weighting function species '{}' was not in the species list".format(a_species)
                        raise ValueError(err_msg)
                    wf_species_handles.append(self.atmosphere.species[a_species].species)
            self.options['wfclimatologies'] = " ".join(wf_species_handles)
            if not ('wfaltitudes' in self.options):
                self.options['wfaltitudes'] = self._wf_altitudes()
            if hasattr(self.atmosphere.wf_species, '__iter__') and not isinstance(self.atmosphere.wf_species, str):
                self._brdf_wf_included = False
                if 'brdf' in self.atmosphere.wf_species:
                    self._brdf_wf_included = True
                    num_wf_species = len(self.atmosphere.wf_species) - 1
                else:
                    num_wf_species = len(self.atmosphere.wf_species)
            else:
                num_wf_species = 1

            if num_wf_species > 1:
                self._wf_shape = (
                    len(self.wavelengths),
                    len(self.geometry.lines_of_sight),
                    num_wf_species,
                    len(self.options['wfaltitudes'])
                )
            else:
                self._wf_shape = (
                    len(self.wavelengths),
                    len(self.geometry.lines_of_sight),
                    len(self.options['wfaltitudes'])
                )
            if self._num_stokes > 1:
                self._wf_shape += (self._num_stokes,)

    @property
    def model_parameters(self):
        return self._model_parameters

    @wrap_skif_functionfail
    def _add_lines_of_sight(self):
        for los in to_iter(self._geometry.lines_of_sight):
            self._iskengine.AddLineOfSight(los.mjd, los.observer, los.look_vector)

        if self._geometry.sun is not None:
            self.options['sun'] = self._geometry.sun

    @wrap_skif_functionfail
    def calculate_radiance(self, output_format='xarray'):
        if self._num_stokes > 1:
            rad = super().calculate_radiance(output_format, full_stokes_vector=True)
            if output_format.lower() == 'xarray':
                stokes = np.hstack((rad['I'].values, rad['Q'].values, rad['U'].values, rad['V'].values)).reshape(len(self.wavelengths), 4, len(self.geometry.lines_of_sight))
                rad['radiance'] = (['wavelength', 'los', 'stokes'], stokes[:, :self.num_stokes, :].transpose([0, 2, 1]))
                rad['stokes'] = ['I', 'Q', 'U', 'V'][:self.num_stokes]
        else:
            rad = super().calculate_radiance(output_format, full_stokes_vector=False)

        if self._atmosphere.wf_species is not None:
            include_brdf = False
            if hasattr(self.atmosphere.wf_species, '__iter__'):
                if 'brdf' in self.atmosphere.wf_species:
                    include_brdf = True
            if self._num_stokes > 1:
                raw_wf = self._iskengine.GetWeightingFunctions()[1].reshape(len(self.wavelengths), len(self.geometry.lines_of_sight), -1, self._num_stokes)
            else:
                raw_wf = self._iskengine.GetWeightingFunctions()[1].reshape(len(self.wavelengths), len(self.geometry.lines_of_sight), -1)
            if include_brdf:
                wf_brdf = raw_wf[:, :, -1]
                if raw_wf.shape[2] > 1:
                    wf = raw_wf[:, :, :-1].reshape(self._wf_shape)
                else:
                    wf = None
            else:
                wf = raw_wf.reshape(self._wf_shape)

            if output_format.lower() == 'numpy':
                output = namedtuple('EngineOutput', ['radiance', 'weighting_function'])
                return output(radiance=rad, weighting_function=wf)
            elif output_format.lower() == 'xarray':
                if self.num_stokes > 1:
                    wf_dims = ['wavelength', 'los', 'perturbation', 'stokes']
                    brdf_wf_dims = ['wavelength', 'los', 'stokes']
                else:
                    wf_dims = ['wavelength', 'los', 'perturbation']
                    brdf_wf_dims = ['wavelength', 'los']

                if hasattr(self.atmosphere.wf_species, '__iter__') and not isinstance(self.atmosphere.wf_species, str):
                    if self._brdf_wf_included:
                        rad['wf_brdf'] = (brdf_wf_dims, wf_brdf)
                        num_species = len(self.atmosphere.wf_species) - 1
                    else:
                        num_species = len(self.atmosphere.wf_species)
                    wf_species = [species for species in self.atmosphere.wf_species if species != 'brdf']
                    for idx, species in enumerate(wf_species):
                        if num_species > 1:
                            rad['wf_' + species] = (wf_dims, wf[:, :, idx])
                        else:
                            rad['wf_' + species] = (wf_dims, wf)
                else:
                    if self.atmosphere.wf_species == 'brdf':
                        rad['wf_' + self.atmosphere.wf_species] = (brdf_wf_dims, wf_brdf)
                    else:
                        rad['wf_' + self.atmosphere.wf_species] = (wf_dims, wf)
                rad.coords['altitude'] = ('altitude', self.options.get('wfaltitudes', np.linspace(500, 99500, 100)))
                rad = rad.stack(perturbation=('altitude',))

                if self.options.get('diagnostics', 0) == 1:
                    od = np.zeros((len(self.wavelengths), len(self.geometry.lines_of_sight), self.num_layers))
                    ssa = np.zeros((len(self.wavelengths), len(self.geometry.lines_of_sight), self.num_layers))
                    heights = np.zeros((len(self.wavelengths), len(self.geometry.lines_of_sight), self.num_layers + 1))
                    observer_od = np.zeros((len(self.wavelengths), len(self.geometry.lines_of_sight)))
                    for w in range(len(self.wavelengths)):
                        for l in range(len(self.geometry.lines_of_sight)):
                            ok, od[w, l, :] = self.skif_object().GetProperty('daopticaldepths[%i, %i]' % (w, l))
                            ok, ssa[w, l, :] = self.skif_object().GetProperty('dassa[%i, %i]' % (w, l))
                            ok, heights[w, l, :] = self.skif_object().GetProperty('daboundaryaltitudes[%i, %i]' % (w, l))
                            ok, observer_od[w, l] = self.skif_object().GetProperty('observerod[%i, %i]' % (w, l))

                    rad['layer_optical_depth'] = (['wavelength', 'los', 'layer'], od)
                    rad['layer_ssa'] = (['wavelength', 'los', 'layer'], ssa)
                    rad['layer_boundary_heights'] = (['wavelength', 'los', 'layer_b'], heights)
                    rad['observer_od'] = (['wavelength', 'los'], observer_od)

                los_transmission = self.skif_object().GetProperty('lostransmission')[1].reshape((len(self.wavelengths), -1))
                rad['los_transmission'] = (['wavelength', 'los'], los_transmission)

                if self.options.get('outputopticaldepth', 0) == 1:
                    layer_od = self.skif_object().GetProperty('layeropticaldepth')[1].reshape((len(self.wavelengths), -1))
                    rad['layer_optical_depth'] = (['wavelength', 'layer'], layer_od)

                    layer_ssa = self.skif_object().GetProperty('layerssa')[1].reshape((len(self.wavelengths), -1))
                    rad['layer_ssa'] = (['wavelength', 'layer'], layer_ssa)

                return rad
        else:
            return rad

    @wrap_skif_functionfail
    def get_reference_points(self):
        reference_points = np.ndarray((len(self._geometry.lines_of_sight), 4))
        for i, los in enumerate(self._geometry.lines_of_sight):
            ok, ref_pt = self._iskengine.GetProperty('referencepoint[{}]'.format(i))
            reference_points[i, :] = ref_pt
        return reference_points

    @property
    def num_streams(self):
        """
        The number of streams used inside the calculation.  Together with the number of layers, the number of streams
        is one of the primary settings that controls the accuracy of the model.
        Default: 16
        """
        return self.model_parameters.get('numstreams', 16)

    @num_streams.setter
    def num_streams(self, n: int):
        self.options['numstreams'] = n

    @property
    def num_layers(self):
        return self.model_parameters.get('numlayers', 50)

    @num_layers.setter
    def num_layers(self, n: int):
        """
        The number of layers used inside the calculation.  Together with the number of streams, the number of layers
        is one of the primary settings that controls the accuracy of the model.
        Default: 50
        """
        self.options['numlayers'] = n

    @property
    def num_brdf_quadrature_terms(self):
        """
        The number of quadrature terms used to calculate the Legendre coefficients for the BRDF.  Only has an effect
        if a Non-Lambertian surface is being used.
        Default: 64
        """
        return self.model_parameters.get('numbrdfexpansionterms')

    @num_brdf_quadrature_terms.setter
    def num_brdf_quadrature_terms(self, n: int):
        self.options['numbrdfexpansionterms'] = n

    @property
    def num_stokes(self):
        """
        Number of stokes parameters to include in the calculation and output.  Set to 1 for scalar calculations, and
        to 3 for vector calculations.
        Default: 1
        """
        return self._num_stokes

    @num_stokes.setter
    def num_stokes(self, n: int):
        self._num_stokes = n

    @property
    def alt_grid(self):
        """
        Internal altitude grid inside the model.  ALl climatologies/optical properties are sampled on this grid, and
        then integrated over it to find the layer quantities.
        Default: Linearly spaced from 0 to 100 km with spacing of 0.5 km
        """
        return self.options.get('altitudegrid', np.arange(0, 100001, 500))

    @alt_grid.setter
    def alt_grid(self, alts: np.array):
        self.options['altitudegrid'] = alts

    @property
    def layer_construction(self):
        """
        Method used to place the layer boundaries in altitude.  Should be one of `uniform_pressure` for
        uniform pressure layers, `uniform_optical_depth` for layers of equal optical depth, `uniform_height` for layers
        spaced linearly in height, `match_altitude_grid` to match layers with the internal altitude grid, or a numpy
        array of altitudes for manually specified layer boundaries.
        Default: uniform_pressure
        """
        return None

    @layer_construction.setter
    def layer_construction(self, layers: Union[str, np.array]):
        if type(layers) == str:
            if layers.lower() == 'uniform_pressure':
                self._options['layerconstructionmethod'] = 0
            elif layers.lower() == 'uniform_height':
                self._options['layerconstructionmethod'] = 2
            elif layers.lower() == 'match_altitude_grid':
                self._options['layerconstructionmethod'] = 3
                self.num_layers = len(self.alt_grid) - 1
        else:
            self.options['manuallayeraltitudes'] = layers
            self.num_layers = len(layers) - 1

    @property
    def num_spherical_sza(self):
        """
        Number of solar zenith angles to calculate the plane parallel solution at when operating in spherical mode.
        Default: 2
        """
        return self.options.get('numsphericalsza', 2)

    @num_spherical_sza.setter
    def num_spherical_sza(self, sza: int):
        self.options['numsphericalsza'] = sza

    @property
    def viewing_mode(self):
        """
        Geometry viewing mode of the model for the line of sight rays, either `spherical` or `plane_parallel`.
        Default: plane_parallel
        """
        mode = self.options.get('uselosspherical', 0)

        if mode == 0:
            return 'plane_parallel'
        else:
            return 'spherical'

    @viewing_mode.setter
    def viewing_mode(self, mode: str):
        if mode.lower() == 'spherical':
            self.options['uselosspherical'] = 1
        else:
            self.options['uselosspherical'] = 0

    @property
    def num_threads(self):
        """
        Number of threads to use in the calculation when running on multiple wavelengths.  The value of 0 indicates
        use every logical core present on the machine.  It is recommended to set this value to the number of physical
        cores on the machine for maximum performance.
        Default: 0
        """
        return self.options.get('numthreads', 0)

    @num_threads.setter
    def num_threads(self, threads: int):
        self.options['numthreads'] = threads
