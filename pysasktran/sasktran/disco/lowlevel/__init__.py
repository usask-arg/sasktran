import sys
import os
import platform
import numpy as np
from pathlib import Path
import ctypes
import xarray as xr


def shared_library_path() -> Path:
    """
    Returns
    -------
    Path
        Full path to the shared library object for the isk_disco object
    """
    filename = sys.modules["sasktran_core"].__file__
    dirname = os.path.dirname(filename)

    if platform.system() == 'Windows':
        dllname = dirname + os.sep + '_sasktran_core_internals.dll'
    elif platform.system() == 'Darwin':
        dllname = dirname + os.sep + 'lib_sasktran_core_internals.dylib'
    else:
        dllname = dirname + os.sep + 'lib_sasktran_core_internals.so'

    path = Path(dllname)
    if not path.exists():
        raise EnvironmentError('Could Not Find sasktran shared library, tried {}'.format(path.as_posix()))

    return path


class _Atmosphere(ctypes.Structure):
    """
    Internal C structure representing the lowlevel C atmosphere object
    """
    _fields_ = [("od", ctypes.POINTER(ctypes.c_double)),
                ("ssa", ctypes.POINTER(ctypes.c_double)),
                ("f", ctypes.POINTER(ctypes.c_double)),
                ("a1", ctypes.POINTER(ctypes.c_double)),
                ("a2", ctypes.POINTER(ctypes.c_double)),
                ("a3", ctypes.POINTER(ctypes.c_double)),
                ("a4", ctypes.POINTER(ctypes.c_double)),
                ("b1", ctypes.POINTER(ctypes.c_double)),
                ("b2", ctypes.POINTER(ctypes.c_double)),
                ("ss_phase", ctypes.POINTER(ctypes.c_double)),
                ("albedo", ctypes.POINTER(ctypes.c_double)),
                ("layerboundaryaltitude", ctypes.POINTER(ctypes.c_double)),
                ("earthradius", ctypes.c_double)
                ]


class Atmosphere(object):
    def __init__(self, nstr: int, nlyr: int, nwavel: int, nlos: int=0, nstokes: int=1):
        """
        Lowlevel atmosphere object, represents the atmosphere as a number of stacked homogeneous layers.
        Each layer is specified by it's total optical depth, single scattering albedo, and legendre coefficients.

        We also require the layer boundary altitudes, which are used to determine the solar attenuation factors
        in a spherical atmosphere.  Similarly the earth radius is needed.

        Additionally the single scatter phase function can be specified so that the engine uses the TMS single
        scattering correction.

        Layers must be specified from top of the atmosphere down to the surface.  i.e., layer index 0 is at the
        top of the atmosphere.

        The DO model can be run in parallel, where multiple "wavelengths" are specified at the same time.

        Parameters
        ----------
        nstr: int
            Number of streams used in the model.  This is the full number of streams, so upwelling and downwelling
            will each have nstr/2 directions.
        nlyr: int
            Number of layers used in the model
        nwavel: int
            Number of "wavelengths"
        nlos: int
            Number of lines of sight (only needed if specifying the single scatter phase function for each line of sight)
            Optional.
        nstokes: int
            Number of stokes parameters for the calculation (only needed if specifying the single scatter phase function)
            Optional.
        """
        self._od = np.zeros((nlyr, nwavel), dtype=np.float64, order='F')
        self._ssa = np.zeros((nlyr, nwavel), dtype=np.float64, order='F')
        self._f = np.zeros((nlyr, nwavel), dtype=np.float64, order='F')
        self._a1 = np.zeros((nstr, nlyr, nwavel), dtype=np.float64, order='F')
        self._a2 = np.zeros((nstr, nlyr, nwavel), dtype=np.float64, order='F')
        self._a3 = np.zeros((nstr, nlyr, nwavel), dtype=np.float64, order='F')
        self._a4 = np.zeros((nstr, nlyr, nwavel), dtype=np.float64, order='F')
        self._b1 = np.zeros((nstr, nlyr, nwavel), dtype=np.float64, order='F')
        self._b2 = np.zeros((nstr, nlyr, nwavel), dtype=np.float64, order='F')

        self._ss_phase = np.zeros((nstokes, nlos, nlyr, nwavel), dtype=np.float64, order='F')

        self._albedo = np.zeros((nwavel), dtype=np.float64, order='F')

        self._layerboundaryaltitudes = np.zeros((nlyr), dtype=np.float64, order='F')

        self._earthradius = 6372000.0

    def c_atmosphere(self) -> _Atmosphere:
        return _Atmosphere(np.ctypeslib.as_ctypes(self._od.flatten('F')),
                           np.ctypeslib.as_ctypes(self._ssa.flatten('F')),
                           np.ctypeslib.as_ctypes(self._f.flatten('F')),
                           np.ctypeslib.as_ctypes(self._a1.flatten('F')),
                           np.ctypeslib.as_ctypes(self._a2.flatten('F')),
                           np.ctypeslib.as_ctypes(self._a3.flatten('F')),
                           np.ctypeslib.as_ctypes(self._a4.flatten('F')),
                           np.ctypeslib.as_ctypes(self._b1.flatten('F')),
                           np.ctypeslib.as_ctypes(self._b2.flatten('F')),
                           np.ctypeslib.as_ctypes(self._ss_phase.flatten('F')),
                           np.ctypeslib.as_ctypes(self._albedo),
                           np.ctypeslib.as_ctypes(self._layerboundaryaltitudes),
                           self._earthradius
                           )

    @property
    def od(self):
        """
        Returns
        -------
        np.ndarray
            Layer Optical Depths, shape (nlayer, nwavel)
        """
        return self._od

    @property
    def ssa(self):
        """
        Returns
        -------
        np.ndarray
            Layer single scatter albedo, shape (nlayer, nwavel)
        """
        return self._ssa

    @property
    def f(self):
        """
        Returns
        -------
        np.ndarray
            Layer delta scaling coefficient, shape (nlayer, nwavel)
        """
        return self._f

    @property
    def a1(self):
        """
        Returns
        -------
        np.ndarray
            Layer a1 Legendre coefficient, shape (nstr, nlayer, nwavel)
        """
        return self._a1

    @property
    def a2(self):
        """
        Returns
        -------
        np.ndarray
            Layer a2 Legendre coefficient, shape (nstr, nlayer, nwavel)
        """
        return self._a2

    @property
    def a3(self):
        """
        Returns
        -------
        np.ndarray
            Layer a3 Legendre coefficient, shape (nstr, nlayer, nwavel)
        """
        return self._a3

    @property
    def a4(self):
        """
        Returns
        -------
        np.ndarray
            Layer a4 Legendre coefficient, shape (nstr, nlayer, nwavel)
        """
        return self._a4

    @property
    def b1(self):
        """
        Returns
        -------
        np.ndarray
            Layer b1 Legendre coefficient, shape (nstr, nlayer, nwavel)
        """
        return self._b1

    @property
    def b2(self):
        """
        Returns
        -------
        np.ndarray
            Layer b2 Legendre coefficient, shape (nstr, nlayer, nwavel)
        """
        return self._b2

    @property
    def ss_phase(self):
        """
        Returns
        -------
        np.ndarray
            Single scatter phase function, shape (nstokes, nlos, nlayer, nwavel)
        """
        return self._ss_phase

    @property
    def albedo(self):
        """
        Returns
        -------
        np.ndarray
            Surface albedo, shape (nwavel)
        """
        return self._albedo

    @property
    def layer_boundary_altitudes(self):
        """
        Returns
        -------
        np.ndarray
            Boundary altitudes of the layer, shape (nlyr+1).
        """
        return self._layerboundaryaltitudes

    @property
    def earth_radius(self):
        return self._earthradius

    @earth_radius.setter
    def earth_radius(self, radius: float):
        self._earthradius = radius


class _WeightingFunctions(ctypes.Structure):
    _fields_ = [("d_od", ctypes.POINTER(ctypes.c_double)),
                ("d_ssa", ctypes.POINTER(ctypes.c_double)),
                ("d_f", ctypes.POINTER(ctypes.c_double)),
                ("d_a1", ctypes.POINTER(ctypes.c_double)),
                ("d_a2", ctypes.POINTER(ctypes.c_double)),
                ("d_a3", ctypes.POINTER(ctypes.c_double)),
                ("d_a4", ctypes.POINTER(ctypes.c_double)),
                ("d_b1", ctypes.POINTER(ctypes.c_double)),
                ("d_b2", ctypes.POINTER(ctypes.c_double)),
                ("d_ss_phase", ctypes.POINTER(ctypes.c_double)),
                ("d_albedo", ctypes.POINTER(ctypes.c_double)),
                ("d_layerindex", ctypes.POINTER(ctypes.c_int)),
                ("numderiv", ctypes.c_int)
                ]


class WeightingFunctions(object):
    def __init__(self, nstr: int, nlos: int, nstokes: int, nwavel: int, nderiv: int):
        """
        Weighting function inputs for the model.  Each layer

        Parameters
        ----------
        nstr
        nlos
        nstokes
        nwavel
        nderiv
        """
        self._d_od = np.zeros((nderiv, nwavel), dtype=np.float64, order='F')
        self._d_ssa = np.zeros((nderiv, nwavel), dtype=np.float64, order='F')
        self._d_f = np.zeros((nderiv, nwavel), dtype=np.float64, order='F')
        self._d_a1 = np.zeros((nderiv, nstr, nwavel), dtype=np.float64, order='F')
        self._d_a2 = np.zeros((nderiv, nstr, nwavel), dtype=np.float64, order='F')
        self._d_a3 = np.zeros((nderiv, nstr, nwavel), dtype=np.float64, order='F')
        self._d_a4 = np.zeros((nderiv, nstr, nwavel), dtype=np.float64, order='F')
        self._d_b1 = np.zeros((nderiv, nstr, nwavel), dtype=np.float64, order='F')
        self._d_b2 = np.zeros((nderiv, nstr, nwavel), dtype=np.float64, order='F')

        self._d_ss_phase = np.zeros((nstokes, nderiv, nlos, nwavel), dtype=np.float64, order='F')

        self._d_albedo = np.zeros((nderiv, nwavel), dtype=np.float64, order='F')

        self._d_layerindex = np.zeros((nderiv), dtype=np.int32, order='F')

        self._nderiv = nderiv

    def c_weighting_functions(self):
        return _WeightingFunctions(np.ctypeslib.as_ctypes(self._d_od.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_ssa.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_f.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_a1.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_a2.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_a3.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_a4.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_b1.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_b2.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_ss_phase.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_albedo.flatten('F')),
                                   np.ctypeslib.as_ctypes(self._d_layerindex),
                                   self._nderiv
                                   )

    @property
    def nderiv(self):
        return self._nderiv

    @property
    def d_od(self):
        return self._d_od

    @property
    def d_ssa(self):
        return self._d_ssa

    @property
    def d_f(self):
        return self._d_f

    @property
    def d_a1(self):
        return self._d_a1

    @property
    def d_albedo(self):
        return self._d_albedo

    @property
    def d_ss_phase(self):
        return self._d_ss_phase

    @property
    def d_layerindex(self):
        return self._d_layerindex


class _Config(ctypes.Structure):
    _fields_ = [('nstr', ctypes.c_int),
                ('nwavel', ctypes.c_int),
                ('nlyr', ctypes.c_int),
                ('nstokes', ctypes.c_int),
                ('nthreads', ctypes.c_int),
                ('useexactsinglescatter', ctypes.c_bool),
                ('usepseudospherical', ctypes.c_bool),
                ('numazimuthexpansion', ctypes.c_int)]


class Config(object):
    def __init__(self, nstr: int, nwavel: int, nlyr: int, nstokes: int, nthreads: int, exact_ss: bool = False,
                 use_pseudo_spherical: bool = True,
                 num_azimuth_expansion: int = 0
                 ):
        self._nstr = nstr
        self._nwavel = nwavel
        self._nlyr = nlyr
        self._nstokes = nstokes
        self._nthreads = nthreads
        self._exact_ss = exact_ss
        self._use_pseudo_spherical = use_pseudo_spherical
        self._num_azimuth_expansion = num_azimuth_expansion

    def c_config(self):
        return _Config(self._nstr, self._nwavel, self._nlyr, self._nstokes, self._nthreads, self._exact_ss,
                       self._use_pseudo_spherical,
                       self._num_azimuth_expansion
                       )

    @property
    def nstr(self):
        return self._nstr

    @property
    def nwavel(self):
        return self._nwavel

    @property
    def nlyr(self):
        return self._nlyr

    @property
    def nstokes(self):
        return self._nstokes


class _ViewingGeometry(ctypes.Structure):
    _fields_ = [('cos_vza', ctypes.POINTER(ctypes.c_double)),
                ('cos_sza', ctypes.c_double),
                ('saa', ctypes.POINTER(ctypes.c_double)),
                ('viewingaltitude', ctypes.POINTER(ctypes.c_double)),
                ('nlos', ctypes.c_int)]


class ViewingGeometry(object):
    def __init__(self, nlos: int):
        self._cos_vza = np.zeros((nlos), dtype=np.float64, order='F')
        self._cos_sza = 0.0
        self._saa = np.zeros((nlos), dtype=np.float64, order='F')
        self._viewingaltitude = np.ones((nlos), dtype=np.float64, order='F') * -1
        self._nlos = nlos

    def c_viewing_geometry(self):
        return _ViewingGeometry(np.ctypeslib.as_ctypes(self._cos_vza),
                                self._cos_sza,
                                np.ctypeslib.as_ctypes(self._saa),
                                np.ctypeslib.as_ctypes(self._viewingaltitude),
                                self._nlos
                                )

    @property
    def cos_vza(self):
        return self._cos_vza

    @property
    def cos_sza(self):
        return self._cos_sza

    @cos_sza.setter
    def cos_sza(self, sza: float):
        self._cos_sza = sza

    @property
    def saa(self):
        return self._saa

    @property
    def viewingaltitude(self):
        return self._viewingaltitude

    @property
    def nlos(self):
        return self._nlos


class _Output(ctypes.Structure):
    _fields_ = [('radiance', ctypes.POINTER(ctypes.c_double)),
                ('d_radiance', ctypes.POINTER(ctypes.c_double))]


class Output(object):
    def __init__(self, nstokes: int, nlos: int, nwavel: int, nderiv: int):
        # For the output we have to have internal buffers that are 1d that we can then map to multiple dimensions
        # afterwards
        self._radiance = np.zeros((nstokes * nlos * nwavel), dtype=np.float64, order='F')
        self._d_radiance = np.zeros((nderiv * nstokes * nlos * nwavel), dtype=np.float64, order='F')

        self._nstokes = nstokes
        self._nlos = nlos
        self._nwavel = nwavel
        self._nderiv = nderiv

    def c_output(self):
        return _Output(np.ctypeslib.as_ctypes(self._radiance),
                       np.ctypeslib.as_ctypes(self._d_radiance))

    def xarray(self):
        return xr.Dataset({'radiance': (['stokes', 'los', 'wavelength'], self._radiance.reshape((self._nstokes, self._nlos, self._nwavel), order='F')),
                           'd_radiance': (['wf', 'stokes', 'los', 'wavelength'], self._d_radiance.reshape((self._nderiv, self._nstokes, self._nlos, self._nwavel), order='F'))
                           })


def calculate(atmosphere: Atmosphere, config: Config, weightingfunctions: WeightingFunctions, geometry: ViewingGeometry):
    dll = ctypes.CDLL(shared_library_path().as_posix())

    atm_pointer = ctypes.POINTER(_Atmosphere)
    conf_pointer = ctypes.POINTER(_Config)
    wf_pointer = ctypes.POINTER(_WeightingFunctions)
    geo_pointer = ctypes.POINTER(_ViewingGeometry)
    out_pointer = ctypes.POINTER(_Output)

    dll.calculate.argtypes = [atm_pointer, conf_pointer, wf_pointer, geo_pointer, out_pointer]

    output = Output(config.nstokes, geometry.nlos, config.nwavel, weightingfunctions.nderiv)

    c_atm = atmosphere.c_atmosphere()
    c_conf = config.c_config()
    c_wf = weightingfunctions.c_weighting_functions()
    c_geo = geometry.c_viewing_geometry()
    c_out = output.c_output()

    dll.calculate(c_atm,
                  c_conf,
                  c_wf,
                  c_geo,
                  c_out)

    return output


if __name__ == "__main__":
    nstr = 4
    nlyr = 2
    nwavel = 100000
    nderiv = 0
    nstokes = 1
    nlos = 1

    atmosphere = Atmosphere(nstr, nlyr, nwavel)
    config = Config(nstr, nwavel, nlyr, nstokes, 0)
    weightingfunctions = WeightingFunctions(nstr, nlyr, nwavel, nderiv)
    viewing_geometry = ViewingGeometry(nlos)

    atmosphere.od[0, :] = 0.2
    atmosphere.od[1, :] = 0.2

    atmosphere.ssa[0, :] = 0.8
    atmosphere.ssa[1, :] = 0.7

    atmosphere.a1[0, 0, :] = 1
    atmosphere.a1[0, 1, :] = 1
    atmosphere.a1[2, 0, :] = 0.5
    atmosphere.a1[2, 1, :] = 0.5

    atmosphere.layer_boundary_altitudes[0] = 100000
    atmosphere.layer_boundary_altitudes[1] = 10000

    viewing_geometry.cos_sza = 0.8
    viewing_geometry.cos_vza[0] = 0.7
    viewing_geometry.saa[0] = 0

    output = calculate(atmosphere, config, weightingfunctions, viewing_geometry)

    x = 5