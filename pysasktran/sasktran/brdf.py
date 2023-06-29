import sasktranif.sasktranif as skif
from sasktran.exceptions import wrap_skif_functionfail
from typing import List
import numpy as np
import logging


class BRDF(object):
    """
    Class which implements a Bidirectional Reflectance Distribution Function.  This is a light wrapper around the
    sasktranif ISKBrdf object.  It is recommended to use the specialized BRDF classes rather than this one.

    A BRDF is defined as the ratio of the radiance scattered from a surface into a given direction to
    the collimated power incident on a unit area of the surface. Much of the literature on BRDFs (Spurr 2002, Roujean
    1992, Rahman 1993, Wanner 1995, Cox and Munk 1954, Kokhanovsky 2012) define reflectance functions through a
    quantity that is scaled up by pi compared to this definition.

    Sasktran defines the relative azimuth angle such that 0 degrees corresponds to forward scattering. Much of the
    literature on BRDFs (Spurr 2002, Roujean 1992, Rahman 1993, Wanner 1995, Cox and Munk 1954, Hapke 2012) define it
    the opposite way.

    Parameters
    ----------
    name: str
        Name of the BRDF object to create

    Examples
    --------
    >>> from sasktran import BRDF
    >>> brdf = BRDF('lambertian')
    >>> brdf.skif_object().SetProperty('Albedo', 0.3)
    True
    >>> brdf.reflectance(340, 0, 0, 54372, 0.5, 0.5, 0)
    0.0954929658551372
    """

    @wrap_skif_functionfail
    def __init__(self, name: str):
        self._name = name
        self._iskbrdf = skif.ISKBrdf(name.upper())

    def __getstate__(self):
        state = {k: v for k, v in self.__dict__.items() if k != '_iskbrdf'}
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._iskbrdf = skif.ISKBrdf(self._name.upper())
        self._update_skif_object()

    def __repr__(self):
        repr = ('BRDF Object: {}'.format(self._name))

        return repr

    def skif_object(self, **kwargs):
        """
        Returns the internel SasktranIF object
        """
        return self._iskbrdf

    @wrap_skif_functionfail
    def reflectance(self, wavel_nm: float, latitude: float, longitude: float, mjd: float, mu_in: float, mu_out: float,
                    cos_dphi: float):
        """
        Calculates the value of the BRDF for a given wavelength, location, time and geometry.

        Parameters
        ----------
        wavel_nm: float
            Wavelength in nanometers.
        latitude: float
            Latitude of surface location.
        longitude: float
            Longitude of surface location.
        mjd: float
            Modified Julian Date
        mu_in: float
            Cosine of the incoming zenith angle. Must be between 0 and 1
        mu_out: float
            Cosine of the outgoing zenith angle. Must be between 0 and 1
        cos_dphi: float
            Cosine of the azimuth angle between incoming and outgoing directions. Must be between -1 (backscattering)
            and +1 (forward scattering).
        """
        pt = [latitude, longitude, 0, mjd]
        return self._iskbrdf.BRDF(wavel_nm, pt, mu_in, mu_out, cos_dphi)[1]

    def _update_skif_object(self):
        pass


class Lambertian(BRDF):
    """
    Class which implements a Lambertian Bidirectional Reflectance Distribution Function.

    Parameters
    ----------
    albedo: float
        Ratio of incoming to outgoing radiance.  The reflectance distribution is equal to albedo/pi

    Examples
    --------
    >>> from sasktran.brdf import Lambertian
    >>> brdf = Lambertian(albedo=0.3)
    >>> brdf.reflectance(340, 0, 0, 54372, 0.5, 0.5, 0)
    0.0954929658551372
    >>> brdf.albedo = 1
    >>> brdf.reflectance(340, 0, 0, 54372, 0.5, 0.5, 0)
    0.3183098861837907
    """

    @wrap_skif_functionfail
    def __init__(self, albedo):
        super().__init__('lambertian')
        self.albedo = albedo

    @property
    def albedo(self):
        return self._albedo

    @albedo.setter
    def albedo(self, value):
        self._albedo = value
        self._update_skif_object()

    @wrap_skif_functionfail
    def _update_skif_object(self):
        self._iskbrdf.SetProperty('Albedo', self._albedo)


class Kokhanovsky(BRDF):
    """
    Implementation of the snow BRDF from Kokhanovsky 2012.  The parameters L and M are optional, the default values
    correspond to 1020 nm over greenland.

    Equation 1 is implemented with the following changes:

    * the reflectance is divided by pi

    Parameters
    ----------
    L: float, optional
        Parameter approximately equal to 13 times the average snow grain optical diameter. Default 3.6e6

    M: float, optional
        Parameter proportional to the mass concentration of pollutants in the snow. Default 5.5e-8

    Examples
    --------
    >>> from sasktran.brdf import Kokhanovsky
    >>> brdf = Kokhanovsky()
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    0.25185744490707646

    References
    ----------
    A. A. Kokhanovsky and F. M. Breon, "Validation of an Analytical Snow BRDF Model
    Using PARASOL Multi-Angular and Multispectral Observations,"
    IEEE Geoscience and Remote Sensing Letters, vol. 9, no. 5, pp. 928-932, Sept. 2012.
    doi: 10.1109/LGRS.2012.2185775
    URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6166850&isnumber=6205680
    """

    def __init__(self, L: float=3.6e6, M: float=5.5e-8):
        self._L = L
        self._M = M
        super().__init__('SNOW_KOKHANOVSKY2012')

        self.L = L
        self.M = M

    @property
    def L(self):
        return self._L

    @L.setter
    def L(self, value):
        self._L = value
        self._update_skif_object()

    @property
    def M(self):
        return self._M

    @M.setter
    def M(self, value):
        self._M = value
        self._update_skif_object()

    @wrap_skif_functionfail
    def _update_skif_object(self):
        self._iskbrdf.SetProperty('BRDFParameters', [self._L, self._M])


class Roujean(BRDF):
    """
    Implementation of the BRDF from Roujean 1992 for heterogenous surfaces. Pre-defined values are given for 11
    different surfaces in the visible and near infrared bands.

    Equation 10 is implemented with values of k0, k1, and k2 from Table 1, but with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi
    * the reflectance is also divided by 100 (k0, k1, k2 are percentages)

    Parameters
    ----------
    surface: str
        One of the following: 'plowed field', 'annual grass', 'hard wheat', 'steppe', 'corn', 'orchard grass',
        'irrigated wheat', 'pineforest', 'deciduous forest', 'soybean', or 'grass lawn'
    spectral_region: str, optional
        One of the following: 'vis' for the visible spectrum or 'nir' for the near-infrared spectrum. Default 'vis'

    Examples
    --------
    >>> from sasktran.brdf import Roujean
    >>> brdf = Roujean(surface='corn', spectral_region='nir')
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    0.09605247619742084

    References
    ----------
    Roujean, J.-L., M. Leroy, and P.-Y. Deschamps (1992), A bidirectional
    reflectance model of the Earth's surface for the correction of remote
    sensing data, J. Geophys. Res., 97(D18), 20455–20468, doi:10.1029/92JD01411.
    """
    _code_map = {'PLOWED FIELD': 1,
                 'ANNUAL GRASS': 2,
                 'HARD WHEAT': 3,
                 'STEPPE': 4,
                 'CORN': 5,
                 'ORCHARD GRASS': 6,
                 'IRRIGATED WHEAT': 7,
                 'PINEFOREST': 8,
                 'DECIDUOUS FOREST': 9,
                 'SOYBEAN': 10,
                 'GRASS LAWN': 11}

    def __init__(self, surface, spectral_region='vis'):
        super().__init__('Roujean')
        self._surface = surface
        self._spectral_region = spectral_region
        self._update_skif_object()

    @property
    def surface(self):
        return self._surface

    @surface.setter
    def surface(self, value):
        self._surface = value
        self._update_skif_object()

    @property
    def spectral_region(self):
        return self._spectral_region

    @spectral_region.setter
    def spectral_region(self, value):
        self._spectral_region = value
        self._update_skif_object()

    @wrap_skif_functionfail
    def _update_skif_object(self):
        code_id = self._code_map.get(self._surface.upper())
        if code_id is None:
            raise ValueError('Surface is not valid')

        if self._spectral_region.lower() == 'nir':
            code_id += 20

        self._iskbrdf.SetProperty('SetPredefinedParameters', code_id)


class CoxMunk(BRDF):
    """
    Implementation of the water BRDF based on Cox and Munk 1954 and Spurr 2002.

    Equation A.18 from Spurr is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

    Parameters
    ----------
    wind_speed: float, optional
        Wind speed in m/s. Default 5.0
    water_index: float, optional
        Index of refraction of water. Default 1.334

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.CoxMunk(wind_speed=5.0, water_index=1.334)
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 1)
    0.21213348733518164

    References
    ----------
    Cox C, Munk W. The measurement of the roughness of the sea surface
    from photographs of the Sun glitter. J Opt Soc Ann 1954; 44:838–50.

    Robert J. D. Spurr, "A new approach to the retrieval of surface properties
    from earthshine measurements," Journal of Quantitative Spectroscopy
    & Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
    doi:10.1016/S0022-4073(02)00283-2
    url: http://www.sciencedirect.com/science/article/pii/S0022407302002832

    """
    def __init__(self, wind_speed: float=5.0, water_index: float=1.334):
        super().__init__('cox_munk')
        self._wind = wind_speed
        self._water = water_index
        self._update_skif_object()

    @property
    def wind_speed(self):
        return self._wind

    @wind_speed.setter
    def wind_speed(self, value):
        self._wind = value
        self._update_skif_object()

    @property
    def water_index(self):
        return self._water

    @water_index.setter
    def water_index(self, value):
        self._water = value
        self._update_skif_object()

    def _update_skif_object(self):
        self._iskbrdf.SetProperty('BRDFParameters', [self._wind, self._water])


class Rahman(BRDF):
    """
    Implementation of the Semiempirical Rahman BRDF based on Rahman 1993 and Spurr 2002.

    Equation 2 from Rahman is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

    Parameters
    ----------
    rho0: float
        Empirical parameter related to the overall intensity of reflection (not proportional to single-scatter albedo)
    theta: float
        Empirical parameter that controls the relative amount of backward (-1<=theta<=0) or forward (0<=theta<=1)
        scattering
    k: float
        Emprical parameter that indicates the level of anisotropy of the surface

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.Rahman(rho0=0.1, theta=-0.5, k=1.5)
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    0.01505098240847888

    References
    ----------
    Rahman H, Pinty B, Verstraete M. Coupled surface atmosphere reflectance
    (CSAR) model: 2. Semi-empirical surface model usable with NOAA AVHRR
    data. J Geophys Res 1993; 98:20,791–801.

    Robert J. D. Spurr, "A new approach to the retrieval of surface properties
    from earthshine measurements," Journal of Quantitative Spectroscopy
    & Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
    doi:10.1016/S0022-4073(02)00283-2
    url: http://www.sciencedirect.com/science/article/pii/S0022407302002832
    """
    def __init__(self, rho0: float, theta: float, k: float):
        super().__init__('rahman')
        self._rho0 = rho0
        self._k = k
        self._theta = theta
        self._update_skif_object()

    @property
    def rho0(self):
        return self._rho0

    @rho0.setter
    def rho0(self, value):
        self._rho0 = value
        self._update_skif_object()

    @property
    def theta(self):
        return self._theta

    @theta.setter
    def theta(self, value):
        self._theta = value
        self._update_skif_object()

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, value):
        self._k = value
        self._update_skif_object()

    def _update_skif_object(self):
        self._iskbrdf.SetProperty('BRDFParameters', [self._rho0, self._theta, self._k])


class Hapke(BRDF):
    """
    Implementation of a vegetation BRDF based on Hapke 2012 and Spurr 2002.

    Equation A.10 of Spurr is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

    Parameters
    ----------
    omega: float
        Single scatter albedo
    delta: float
        Angular width at half maximum of the hotspot / opposition peak
    b0: float
        Amplitude of the hotspot / opposition peak

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.Hapke(omega=0.6, delta=0.06, b0=1.0)
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    0.0814968276380363

    References
    ----------
    Hapke B. Theory of reflectance and emittance spectroscopy. Cambridge:
    Cambridge University Press, 2012.

    Robert J. D. Spurr, "A new approach to the retrieval of surface properties
    from earthshine measurements," Journal of Quantitative Spectroscopy
    & Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
    doi:10.1016/S0022-4073(02)00283-2
    url: http://www.sciencedirect.com/science/article/pii/S0022407302002832
    """
    def __init__(self, omega: float, delta: float, b0: float):
        super().__init__('hapke')
        self._omega = omega
        self._delta = delta
        self._b0 = b0
        self._update_skif_object()

    @property
    def omega(self):
        return self._omega

    @omega.setter
    def omega(self, value):
        self._omega = value
        self._update_skif_object()

    @property
    def delta(self):
        return self._delta

    @delta.setter
    def delta(self, value):
        self._delta = value
        self._update_skif_object()

    @property
    def b0(self):
        return self._b0

    @b0.setter
    def b0(self, value):
        self._b0 = value
        self._update_skif_object()

    def _update_skif_object(self):
        self._iskbrdf.SetProperty('BRDFParameters', [self._omega, self._delta, self._b0])


class MODIS(BRDF):
    """
    Implementation of the Ross-thick Li-sparse-reciprocal BRDF used in MODIS' Albedo/BRDF product. It is a linear
    combination of an isotropic (Lambertian) term, a volume scatterer (Ross-thick kernel) and a geometric scatterer
    (Li-sparse-reciprocal kernel).

    Parameters
    ----------
    f_iso: float
        Weight of the isotropic (Lambertian) kernel
    f_vol: float
        Weight of the volume (Ross-thick) kernel
    f_geo: float
        Weight of the geometric (Li-sparse-reciprocal) kernel

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.MODIS(f_iso=0.2, f_vol=0.2, f_geo=0.1)
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    0.031577469878097834

    References
    ---------
    A. H. Strahler, J. P. Muller, "MODIS BRDF/Albedo Product: Algorithm
    Theoretical Basis Document Version 5.0," April 1999.
    URL: https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/atbd_mod09_v5.pdf
    """
    def __init__(self, f_iso, f_vol, f_geo):
        super().__init__('modis')
        self._iso = f_iso
        self._vol = f_vol
        self._geo = f_geo
        self._update_skif_object()

    @property
    def f_isotropic(self):
        return self._iso

    @f_isotropic.setter
    def f_isotropic(self, value):
        self._iso = value
        self._update_skif_object()

    @property
    def f_volume(self):
        return self._vol

    @f_volume.setter
    def f_volume(self, value):
        self._vol = value
        self._update_skif_object()

    @property
    def v_geometric(self):
        return self._geo

    @v_geometric.setter
    def v_geometric(self, value):
        self._geo = value
        self._update_skif_object()

    def _update_skif_object(self):
        self._iskbrdf.SetProperty('BRDFParameters', [self._iso, self._vol, self._geo])


class LinearCombination(BRDF):
    """
    Implements a BRDF that is a linear combination of an arbitrary number of kernels.

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.LinearCombination()

    >>> kernel1 = sk.Lambertian(albedo=1.0)
    >>> kernel2 = sk.RossThickKernel()
    >>> kernel3 = sk.LiSparseReciprocalKernel(crown_shape=1.0, relative_height=2.0)

    >>> brdf.kernels = [kernel1, kernel2, kernel3]
    >>> brdf.kernel_weights = [0.2, 0.2, 0.1]

    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    0.031577469878097834
    """
    def __init__(self):
        super().__init__('linear_combination')
        self._kernels = []
        self._weights = []

    def __getstate__(self):
        state = self.__dict__
        state['_kernels'] = [kernel.__getstate__() for kernel in state['_kernels']]

    def __setstate__(self, state):
        super().__setstate__(state)             # should fill self.__dict__
        self.kernels = self._kernels            # trigger kernel.setter
        self.kernel_weights = self._weights     # trigger kernel_weights.setter

    def add_kernel(self, kernel: BRDF):
        """
        Add a BRDF kernel to the linear combination.

        Parameters
        ----------
        kernel: sk.BRDF
            The BRDF kernel to add.
        """
        self._kernels.append(kernel)
        self._iskbrdf.SetProperty('AddKernel', kernel.skif_object())

    def remove_kernel(self, index: int=-1):
        """
        Remove the BRDF kernel at the position given by index. Accepts negative indices as Python does.

        Parameters
        ----------
        index: int
            Position of the kernel to remove. Invalid indices will remove nothing and produce a warning. Default -1
        """
        if index >= -len(self._kernels) and index < len(self._kernels):
            del self._kernels[index]
            if index >= 0:
                self._iskbrdf.SetProperty('RemoveKernel', index)
            else:
                self._iskbrdf.SetProperty('RemoveKernel', index + len(self._kernels))
        else:
            logging.warning('LinearCombinationBRDF.remove_kernel(): the kernel was not removed because {} is an '
                            'invalid index'.format(index))

    @property
    def kernel_weights(self):
        """
        List of kernel weights. If the number of weights given does not match the number of kernels present, the
        given weights will not be set.
        """
        return self._weights

    @kernel_weights.setter
    def kernel_weights(self, value: List[float]):
        if len(value) == len(self._kernels):
            self._weights = value
            self._iskbrdf.SetProperty('KernelWeights', value)
        else:
            logging.warning('LinearCombinationBRDF.kernel_weights(): the kernel weights were not set because {} '
                            'weights were given for {} kernels'.format(len(value), len(self._kernels)))

    @property
    def kernels(self):
        """
        Full list of BRDF kernels.
        """
        return self._kernels

    @kernels.setter
    def kernels(self, value: List[BRDF]):
        for i in range(len(self._kernels)):
            self.remove_kernel(0)
        for kernel in value:
            self.add_kernel(kernel)
        self._kernels = value


class RoujeanKernel(BRDF):
    """
    Implements the Roujean kernel from Roujean 1992 and Spurr 2002. This kernel is part of the above Roujean BRDF.
    It is recommended to use kernels as part of a LinearCombination and not as standalone BRDFs.

    Equation 2 from Roujean is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.RoujeanKernel()
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    -0.4471903004981418

    References
    ----------
    Robert J. D. Spurr, "A new approach to the retrieval of surface properties
    from earthshine measurements," Journal of Quantitative Spectroscopy
    & Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
    doi:10.1016/S0022-4073(02)00283-2
    url: http://www.sciencedirect.com/science/article/pii/S0022407302002832

    Roujean, J.-L., M. Leroy, and P.-Y. Deschamps (1992), A bidirectional
    reflectance model of the Earth's surface for the correction of remote
    sensing data, J. Geophys. Res., 97(D18), 20455–20468, doi:10.1029/92JD01411.

    """
    def __init__(self):
        super().__init__('roujean_kernel')


class LiKernel(BRDF):
    """
    Base class for Li kernels from Wanner 1995 and Spurr 2002. It is recommended to use derived classes
    (LiSparseKernel, LiDenseKernel or LiSparseReciprocalKernel) instead of this class. It is also recommended to use
    kernels as part of a LinearCombination and not as standalone BRDFs.

    Tree crowns are modelled as spheroids of height 2b and width 2r, with centers a height of h above the ground.

    Parameters
    ----------
    crown_shape: float
        The ratio b/r describing the shape of the crown.
    relative_height: float
        The ratio h/b describing the height of the crown above the ground.

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.LiKernel('li_sparse_kernel', crown_shape=1.0, relative_height=2.0)
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    -0.8753521870054244

    References
    ----------
    W. Wanner and A. H. Strahler, "On the derivation of kernels for kernel-driven
    modes of bidirectional reflectance," Journal of Geophysical Research, vol.
    100, no. d10, pp. 21077-21089, Oct. 20, 1995.
    doi: 10.1029/95JD02371
    URL: http://onlinelibrary.wiley.com/doi/10.1029/95JD02371/abstract

    Robert J. D. Spurr, "A new approach to the retrieval of surface properties
    from earthshine measurements," Journal of Quantitative Spectroscopy
    & Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
    doi:10.1016/S0022-4073(02)00283-2
    url: http://www.sciencedirect.com/science/article/pii/S0022407302002832
    """
    def __init__(self, name: str, crown_shape: float, relative_height: float):
        super().__init__(name)
        self._crown_shape = crown_shape
        self._relative_height = relative_height
        self._update_skif_object()

    @property
    def crown_shape(self):
        return self._crown_shape

    @crown_shape.setter
    def crown_shape(self, value):
        self._crown_shape = value
        self._update_skif_object()

    @property
    def relative_height(self):
        return self._relative_height

    @relative_height.setter
    def relative_height(self, value):
        self._relative_height = value
        self._update_skif_object()

    def _update_skif_object(self):
        self._iskbrdf.SetProperty('BRDFParameters', [self._crown_shape, self._relative_height])


class LiSparseKernel(LiKernel):
    """
    The Li-sparse kernel models sparse canopy where mutual shadowing can be ignored. It is recommended to use
    kernels as part of a LinearCombination and not as standalone BRDFs.

    Equation 32 from Wanner is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

    Parameters
    ----------
    crown_shape: float, optional
        See sk.LiSparse - Default 1.0
    relative_height: float, optional
        See sk.LiSparse - Default 2.0

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.LiSparseKernel(crown_shape=1.0, relative_height=2.0)
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    -0.8753521870054244

    References
    ----------
    W. Wanner and A. H. Strahler, "On the derivation of kernels for kernel-driven
    modes of bidirectional reflectance," Journal of Geophysical Research, vol.
    100, no. d10, pp. 21077-21089, Oct. 20, 1995.
    doi: 10.1029/95JD02371
    URL: http://onlinelibrary.wiley.com/doi/10.1029/95JD02371/abstract

    Robert J. D. Spurr, "A new approach to the retrieval of surface properties
    from earthshine measurements," Journal of Quantitative Spectroscopy
    & Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
    doi:10.1016/S0022-4073(02)00283-2
    url: http://www.sciencedirect.com/science/article/pii/S0022407302002832
    """
    def __init__(self, crown_shape: float=1.0, relative_height: float=2.0):
        super().__init__('li_sparse_kernel', crown_shape, relative_height)


class LiSparseReciprocalKernel(LiKernel):
    """
    The Li-sparse-reciprocal kernel is a variation on the Li-sparse kernel which allows the incident/viewing angles
    to be swapped without changing the BRDF's value (reciprocal) - See the MODIS Albedo/BRDF ATBD. It is recommended
    to use kernels as part of a LinearCombination and not as standalone BRDFs.

    Equation 39 from Strahler is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

    Parameters
    ----------
    crown_shape: float, optional
        See sk.LiSparse - Default 1.0
    relative_height: float, optional
        See sk.LiSparse - Default 2.0

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.LiSparseReciprocalKernel(crown_shape=1.0, relative_height=2.0)
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    -0.477464829275686

    References
    ----------
    A. H. Strahler, J. P. Muller, "MODIS BRDF/Albedo Product: Algorithm
    Theoretical Basis Document Version 5.0," April 1999.
    URL: https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/atbd_mod09_v5.pdf
    """
    def __init__(self, crown_shape=1.0, relative_height=2.0):
        super().__init__('li_sparse_reciprocal_kernel', crown_shape, relative_height)


class LiDenseKernel(LiKernel):
    """
    The Li-dense kernel models dense canopy where mutual shadowing cannot be ignored. It is recommended to use
    kernels as part of a LinearCombination and not as standalone BRDFs.

    Equation 47 from Wanner is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

    Parameters
    ----------
    crown_shape: float, optional
        See sk.LiSparse - Default 2.5
    relative_height: float, optional
        See sk.LiSparse - Default 2.0

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.LiDenseKernel(crown_shape=2.5, relative_height=2.0)
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    -0.46940635114445084

    References
    ----------
    W. Wanner and A. H. Strahler, "On the derivation of kernels for kernel-driven
    modes of bidirectional reflectance," Journal of Geophysical Research, vol.
    100, no. d10, pp. 21077-21089, Oct. 20, 1995.
    doi: 10.1029/95JD02371
    URL: http://onlinelibrary.wiley.com/doi/10.1029/95JD02371/abstract

    Robert J. D. Spurr, "A new approach to the retrieval of surface properties
    from earthshine measurements," Journal of Quantitative Spectroscopy
    & Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
    doi:10.1016/S0022-4073(02)00283-2
    url: http://www.sciencedirect.com/science/article/pii/S0022407302002832
    """
    def __init__(self, crown_shape=2.5, relative_height=2.0):
        super().__init__('li_dense_kernel', crown_shape, relative_height)


class RossThinKernel(BRDF):
    """
    The Ross-thin kernel models canopies with a small leaf area index (LAI) - see Wanner 1995 and Spurr 2002. It is
    recommended to use kernels as part of a LinearCombination and not as standalone BRDFs.

    Equation 13 from Wanner is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.RossThinKernel()
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    0.8132395113781661

    References
    ----------
    Robert J. D. Spurr, "A new approach to the retrieval of surface properties
    from earthshine measurements," Journal of Quantitative Spectroscopy
    & Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
    doi:10.1016/S0022-4073(02)00283-2
    url: http://www.sciencedirect.com/science/article/pii/S0022407302002832

    W. Wanner and A. H. Strahler, "On the derivation of kernels for kernel-driven
    modles of bidirectional reflectance," Journal of Geophysical Research, vol.
    100, no. d10, pp. 21077-21089, Oct. 20, 1995.
    doi: 10.1029/95JD02371
    URL: http://onlinelibrary.wiley.com/doi/10.1029/95JD02371/abstract
    """
    def __init__(self):
        super().__init__('ross_thin_kernel')


class RossThickKernel(BRDF):
    """
    The Ross-thick kernel models canopies with a large leaf area index (LAI) - see Wanner 1995 and Spurr 2002. It is
    recommended to use kernels as part of a LinearCombination and not as standalone BRDFs.

    Equation 7 from Wanner is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

    Examples
    --------
    >>> import sasktran as sk
    >>> brdf = sk.RossThickKernel()
    >>> brdf.reflectance(1050, 0, 0, 54372, 0.5, 0.5, 0)
    0.07830987784454152

    References
    ----------
    Robert J. D. Spurr, "A new approach to the retrieval of surface properties
    from earthshine measurements," Journal of Quantitative Spectroscopy
    & Radiative Transfer, vol. 83, pp. 15-46, Oct. 9, 2002.
    doi:10.1016/S0022-4073(02)00283-2
    url: http://www.sciencedirect.com/science/article/pii/S0022407302002832

    W. Wanner and A. H. Strahler, "On the derivation of kernels for kernel-driven
    modles of bidirectional reflectance," Journal of Geophysical Research, vol.
    100, no. d10, pp. 21077-21089, Oct. 20, 1995.
    doi: 10.1029/95JD02371
    URL: http://onlinelibrary.wiley.com/doi/10.1029/95JD02371/abstract
    """
    def __init__(self):
        super().__init__('ross_thick_kernel')


class LatLonBRDF(BRDF):
    """
    A custom BRDF that is constructed as the combination of other BRDF's on a uniform latitude/longitude grid.

    Parameters
    ----------
    longitudes: np.array
        Longitude grid that the BRDF is specified on in degrees.  Longitudes must be specified from 0 to 360
        and not from -180 to 180.
    latitudes: np.array
        Latitude grid that the BRDF is specified on in degrees
    brdfs: np.ndarray
        2-D Array of dimension (len(longitudes), len(latitudes)) of sk.BRDF objects

    Examples
    --------
    >>> import sasktran as sk
    >>> latitudes = [-10, 10]
    >>> longitudes = [40, 60]
    >>> brdfs = np.empty((len(longitudes), len(latitudes)), dtype=object)
    >>> brdfs[0, 0] = sk.Lambertian(0.1)
    >>> brdfs[0, 1] = sk.Lambertian(0.5)
    >>> brdfs[1, 0] = sk.Lambertian(0.8)
    >>> brdfs[1, 1] = sk.Lambertian(1.0)
    >>> brdf = sk.LatLonBRDF(longitudes, latitudes, brdfs)
    >>> brdf.reflectance(wavel_nm=750, latitude=0, longitude=50, mjd=54372, mu_in=0, mu_out=0, cos_dphi=0)
    0.1909859317102744
    """
    def __init__(self, longitudes: np.array, latitudes: np.array, brdfs: np.ndarray):
        super().__init__('USERDEFINED_LATLON')
        self._longitudes = longitudes
        self._latitudes = latitudes
        self._brdfs = brdfs

        self._update_skif_object()

    @wrap_skif_functionfail
    def _update_skif_object(self):
        # Need to set these properties first as they reallocate the internal grid of BRDFs
        self.skif_object().SetProperty('latitudes', self._latitudes)
        self.skif_object().SetProperty('longitudes', self._longitudes)

        index = 0
        for idy in range(self._brdfs.shape[1]):
            for idx in range(self._brdfs.shape[0]):
                # Set the index of the BRDF we are setting and the BRDF
                self.skif_object().SetProperty('brdfindex', index)
                self.skif_object().SetProperty('brdf', self._brdfs[idx, idy].skif_object())

                index += 1

        # Set a dummy property which triggers initializing of the BRDFs + the grids
        self.skif_object().SetProperty('initialize', 0)


class SpectrallyVaryingBRDF(BRDF):
    """
    A custom BRDF that is composed of a collection of other BRDFs on a wavelength grid to create a spectrally
    varying BRDF.

    Parameters
    ----------
    wavelengths: np.array
        Wavelength grid that the BRDF is specified on in nm
    brdfs: np.ndarray
        1-D Array of dimension (len(wavelengths)) of sk.BRDF objects

    Examples
    --------
    >>> import sasktran as sk
    >>> import numpy as np
    >>> wavelengths = np.array([350, 750])
    >>> brdfs = np.empty((len(wavelengths)), dtype=object)
    >>> brdfs[0] = sk.Lambertian(0.1)
    >>> brdfs[1] = sk.Lambertian(0.5)
    >>> brdf = sk.SpectrallyVaryingBRDF(wavelengths, brdfs)
    >>> brdf.reflectance(wavel_nm=350, latitude=0, longitude=50, mjd=54372, mu_in=0, mu_out=0, cos_dphi=0)
    0.03183098861837907
    >>> brdf.reflectance(wavel_nm=550, latitude=0, longitude=50, mjd=54372, mu_in=0, mu_out=0, cos_dphi=0)
    0.09549296585513721
    """
    def __init__(self, wavelengths: np.array, brdfs: np.ndarray):
        super().__init__('SPECTRAL_VARYING')
        self._wavelengths = wavelengths
        self._brdfs = brdfs

        self._update_skif_object()

    @wrap_skif_functionfail
    def _update_skif_object(self):
        # Need to set these properties first as they reallocate the internal grid of BRDFs
        self.skif_object().SetProperty('wavelengths', self._wavelengths)

        index = 0
        for idx in range(self._brdfs.shape[0]):
            # Set the index of the BRDF we are setting and the BRDF
            self.skif_object().SetProperty('brdfindex', index)
            self.skif_object().SetProperty('brdf', self._brdfs[idx].skif_object())

            index += 1

        # Set a dummy property which triggers initializing of the BRDFs + the grids
        self.skif_object().SetProperty('initialize', 0)
