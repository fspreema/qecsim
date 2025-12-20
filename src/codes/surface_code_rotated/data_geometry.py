from dataclasses import dataclass

from src.codes.surface_code_rotated.get_stab_pairings import SurfacePairings
from src.codes.surface_code_rotated.surface_geom import SurfaceGeometry


@dataclass
class MasterGeometry:
    geometry_std: SurfaceGeometry
    geometry_ybasis: SurfaceGeometry


@dataclass
class MasterPairings:
    pairings_std: SurfacePairings
    pairings_ybasis: SurfacePairings
    pairings_ymemory: SurfacePairings
    pairings_yswitch: SurfacePairings
    pairings_log_h: SurfacePairings
