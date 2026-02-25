from dataclasses import dataclass

from src.codes.surface_code_rotated.get_stab_pairings import SurfacePairings
from src.codes.surface_code_rotated.surface_geom import SurfaceGeometry


@dataclass
class MasterGeometry:
    geometry_std: SurfaceGeometry | None
    geometry_ybasis: SurfaceGeometry | None


@dataclass
class MasterPairings:
    pairings_std: SurfacePairings | None
    pairings_ybasis: SurfacePairings | None
    pairings_ymemory: SurfacePairings | None
    pairings_yswitch: SurfacePairings | None
    pairings_log_h: SurfacePairings | None
