from dataclasses import dataclass
from typing import Optional, List
import numpy as np


@dataclass
class LayerQuantities:
    od: np.ndarray                     # [lyr, wavel]
    ssa: np.ndarray                    # [lyr, wavel]
    a1: Optional[np.ndarray] = None    # [str, lyr, wavel]
    a2: Optional[np.ndarray] = None    # [str, lyr, wavel]
    a3: Optional[np.ndarray] = None    # [str, lyr, wavel]
    a4: Optional[np.ndarray] = None    # [str, lyr, wavel]
    b1: Optional[np.ndarray] = None    # [str, lyr, wavel]
    b2: Optional[np.ndarray] = None    # [str, lyr, wavel]
    ss_phase: Optional[np.ndarray] = None  # [stokes, los, lyr, wavel]


@dataclass
class LayerQuantitiesDerivative:
    d_od: np.ndarray                    # [deriv, wavel]
    d_ssa: np.ndarray                   # [deriv, wavel]
    d_albedo: np.ndarray                # [deriv, wavel]
    layer_index: np.ndarray             # [deriv]
    d_a1: Optional[np.ndarray] = None   # [deriv, str, wavel]
    d_a2: Optional[np.ndarray] = None   # [deriv, str, wavel]
    d_a3: Optional[np.ndarray] = None   # [deriv, str, wavel]
    d_a4: Optional[np.ndarray] = None   # [deriv, str, wavel]
    d_b1: Optional[np.ndarray] = None   # [deriv, str, wavel]
    d_b2: Optional[np.ndarray] = None   # [deriv, str, wavel]


@dataclass
class LayerDefinition:
    nlyr: int
    layer_boundary_altitudes: np.array   # [nlyr+1], altitude relative to earth radius
    earth_radius: float
    layer_boundary_pressures: Optional[np.array]   # [nlyr+1], only required if a species is defined on pressure levels
