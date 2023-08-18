import sasktranif.sasktranif as skif
import numpy as np
import sasktran as sk
import logging
from copy import copy
from sasktran.exceptions import wrap_skif_functionfail
from collections import namedtuple
from sasktran.util import DictWithCallback, LowercaseKeyDictWithCallback, to_iter
from sasktran.stokesvector import StokesVector


class EngineCO(sk.Engine):
    """
    """
    def __init__(self, geometry: sk.Geometry = None, atmosphere: sk.Atmosphere = None,
                 wavelengths: list = None, options: dict = None):
        super().__init__('co', geometry, atmosphere, wavelengths, options)

        self._model_parameters = dict()
        self._wf_shape = None
        self._nstokes = 1

    def _initialize_model(self):
        super()._initialize_model()

        self._model_parameters['referencepoint'] = self._iskengine.GetProperty('ReferencePoint')[1]
        # self._model_parameters['sun'] = self._iskengine.GetProperty('sun')[1]

    def _add_atmosphere(self):
        super()._add_atmosphere()

    def _add_lines_of_sight(self):
        super()._add_lines_of_sight()

        # if self._geometry.reference_point is not None:
        #     self._iskengine.SetProperty('setreferencepoint', self._geometry.reference_point)

    @property
    def model_parameters(self):
        return self._model_parameters

    @wrap_skif_functionfail
    def _add_atmosphere_wf(self):
        if self.atmosphere.wf_species is not None:
            wf_species_handles = []
            for a_species in np.atleast_1d(self.atmosphere.wf_species):
                if not (a_species in self.atmosphere.species):
                    err_msg = "Weighting function species '{}' was not in the species list".format(a_species)
                    if err_msg is not None:
                        raise ValueError(err_msg)
                else:
                    wf_species_handles.append(self.atmosphere.species[a_species].species)
            self.options['wfspecies'] = " ".join(wf_species_handles)

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
                    1,
                    -1
                )

    @wrap_skif_functionfail
    def calculate_radiance(self, output_format='numpy', full_stokes_vector=False, stokes_orientation='geographic'):
        self._stokes_stack = self._nstokes
        if full_stokes_vector:
            raise ValueError('EngineCO does not support full_stokes_vector=True in calculate radiance, instead set engine.nstokes=3')

        rad = super().calculate_radiance(output_format=output_format, full_stokes_vector=full_stokes_vector,
                                         stokes_orientation=stokes_orientation)

        if self._atmosphere.wf_species is not None:
            wf = self._iskengine.GetWeightingFunctions()[1].reshape(self._wf_shape)

            if hasattr(self.atmosphere.wf_species, '__iter__') and not isinstance(self.atmosphere.wf_species, str):
                for species_idx, species in enumerate(self.atmosphere.wf_species):
                    rad['wf_' + species] = (['wavelength', 'los', 'perturbation'], wf[:, :, species_idx, :])
            else:
                rad['wf_' + self.atmosphere.wf_species] = (['wavelength', 'los', 'perturbation'], wf[:, :, 0, :])

        return rad

    @property
    def nstokes(self) -> int:
        return self._nstokes

    @nstokes.setter
    def nstokes(self, ns: int):
        self._nstokes = ns
        self._options['nstokes'] = self._nstokes
