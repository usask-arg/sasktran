from .config import update_registry_from_config
from .geometry import Geometry, VerticalImage, NadirGeometry
from .atmosphere import Atmosphere

# Wrappers for core sasktranif objects
from .climatology import Climatology, ClimatologyUserDefined, ClimatologyUserDefined2D, ClimatologyUserDefined3D, Labow, Pratmo, MSIS90, ECMWF, GloSSAC
from .opticalproperty import OpticalProperty, MieAerosol, HITRANChemical, OpticalPropertyConvolved, UserDefinedAbsorption, UserDefinedAbsorptionPressure
from .opticalproperty import O3DBM, O3OSIRISRes, NO2Vandaele1998, Rayleigh, SimpleRayleigh, InelasticRayleigh, NO2OSIRISRes, BaumIceCrystal, SO2Vandaele2009
from .opticalproperty import O2O2Fally2000, O2O2HITRAN2016, O2O2Thalman2013, UserDefinedScatterConstantHeight
from .engine import Engine, EngineHR, EngineHRSSApprox, EngineMC
from .engine_occ import EngineOCC
from .engine_so import EngineSO
from .geodetic import Geodetic
from .solarspectrum import SolarSpectrum
from .brdf import BRDF, Lambertian, Roujean, Kokhanovsky, CoxMunk, Rahman, Hapke, MODIS
from .brdf import LinearCombination, RoujeanKernel, LiKernel, LiDenseKernel, LiSparseKernel, LiSparseReciprocalKernel
from .brdf import RossThinKernel, RossThickKernel, LatLonBRDF, SpectrallyVaryingBRDF
from .emission import Emission, EmissionTable, EmissionThermal
from .emission_hitranphotochemical import HITRANPhotoChemical_O2_ABand, HITRANPhotoChemical_O2_SingletDelta, HITRANPhotoChemical
from .aband.model import ABandEmission
from .mie import Mie, MieWiscombe
from .tir.engine import EngineTIR
from .disco.engine import EngineDO

# Utility classes overtop of sasktranif objects
from .lineofsight import LineOfSight
from .species import Species, SpeciesAerosol, SpeciesAerosolGloSSAC, SpeciesBaumIceCloud
from .stokesvector import StokesVector

from .handles import standard_handles
from .exceptions import SasktranError
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

update_registry_from_config()
