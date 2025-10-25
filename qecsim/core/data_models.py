from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Mapping, Optional, Any

# Type aliases
Coord = complex
Label = str
Index = int
Pair = Tuple[Coord, Coord]


#--------------------------
# Shared Helper Dataclasses
#--------------------------

@dataclass
class NoiseModel:
    """Group all noise probabilites used by sub-builders"""
    before_round_depol: float = 0.0
    before_m_flip_prob: float = 0.0
    after_r_flip: float = 0.0
    after_c_depol_prob: float = 0.0


@dataclass
class CircuitResult:
    """
    Uniform return type for circuit-building functions.

    Keep this generic so circuit builders can use any circuit object
    """
    circuit: Any
    obs_indices: Optional[List[int]] = None

    def __iadd__(self, other: "CircuitResult | Any") -> "CircuitResult":
        """Allow `result += other` regardless of other’s type."""
        if isinstance(other, CircuitResult):
            # prefer to delegate to circuit's in-place add if available
            try:
                self.circuit += other.circuit
            except Exception:
                # fallback: try simple concatenation/extend where applicable
                try:
                    self.circuit = self.circuit + other.circuit  # type: ignore
                except Exception:
                    raise
        else:
            self.circuit += other  # type: ignore
        return self


#--------------------------
# Configs (surface vs lattice-surgery)
#--------------------------

@dataclass
class Config_Surface:
    """Configuration used by surface/rotated builders"""
    distance: int
    state_init: str
    obs: str
    rounds: int


@dataclass
class Config_LatticeSurgery:
    """Configuration used by lattice-surgery builders"""
    distance: int
    target_state_init: str
    control_state_init: str


#--------------------------
# Internal pick-up helper
#--------------------------

def _pick_up_indices(coords: Dict[Coord, Label], q2i: Mapping[Coord, Index], *labels: str) -> List[Index]:
    """
    Return q2i indices matching any of the provided labels in coords.
    """
    label_set = set(labels)
    return [q2i[q] for q, t in coords.items() if t in label_set]


#--------------------
# Surface-style Patch 
#--------------------

@dataclass
class Patch:
    """All geometry information is stored here i.e. data/x_stab indices"""
    coords: Dict[Coord, Label]
    data: List[Index]
    x_stab: List[Index]
    stab_switch_apply_h: List[Index]
    x_stab_memory: List[Index]
    z_stab_memory: List[Index]
    z_stab: List[Index]
    upper_h: List[Index]
    right_h: List[Index]

    @classmethod
    def from_coords(
        cls,
        coords: Dict[Coord, Label],
        q2i: Mapping[Coord, Index],
    ) -> "Patch":
        def pick_up(*labels: str) -> List[Index]:
            return _pick_up_indices(coords, q2i, *labels)

        return cls(
            coords=coords,
            data=pick_up("DATA"),
            x_stab=pick_up("X-STAB", "X-STAB-BOUND-U", "X-STAB-BOUND-B", "X-STAB-BOUND-R"),
            stab_switch_apply_h=pick_up("X-STAB", "Z-STAB-BOUND-U", "X-STAB-BOUND-B", "Z-STAB-BOUND-U-H"),
            x_stab_memory=pick_up("X-STAB", "Z-STAB-BOUND-U-H", "X-STAB-BOUND-B"),
            z_stab_memory=pick_up("Z-STAB", "Z-STAB-BOUND-L", "X-STAB-BOUND-R-H"),
            z_stab=pick_up("Z-STAB", "Z-STAB-BOUND-L", "Z-STAB-BOUND-R", "Z-STAB-BOUND-U"),
            upper_h=pick_up("Z-STAB-BOUND-U-H"),
            right_h=pick_up("X-STAB-BOUND-R-H"),
        )


#--------------------------
# Lattice-surgery-style patches
#--------------------------

@dataclass
class Patch_Ancilla:
    """All geometry information is stored here i.e. data x_stab indices"""
    coords: Dict[Coord, Label]
    data: List[Index]
    x_stab: List[Index]
    z_stab: List[Index]
    x_bdyB: List[Index]
    z_bdyR: List[Index]

    @classmethod
    def from_coords(
        cls,
        coords: Dict[Coord, Label],
        q2i: Mapping[Coord, Index],
    ) -> "Patch_Ancilla":
        def pick_up(*labels: str) -> List[Index]:
            return _pick_up_indices(coords, q2i, *labels)

        return cls(
            coords=coords,
            data=pick_up("DATA"),
            x_stab=pick_up("X-STAB", "X-STAB-BOUND-A-A", "X-STAB-BOUND-B-A"),
            z_stab=pick_up("Z-STAB", "Z-STAB-BOUND-L-A", "Z-STAB-BOUND-R-A"),
            x_bdyB=pick_up("X-STAB-BOUND-B-A"),
            z_bdyR=pick_up("Z-STAB-BOUND-R-A"),
        )


@dataclass
class Patch_Control:
    coords: Dict[Coord, Label]
    data: List[Index]
    x_stab: List[Index]
    z_stab: List[Index]

    @classmethod
    def from_coords(
        cls,
        coords: Dict[Coord, Label],
        q2i: Mapping[Coord, Index],
    ) -> "Patch_Control":
        def pick_up(*labels: str) -> List[Index]:
            return _pick_up_indices(coords, q2i, *labels)

        return cls(
            coords=coords,
            data=pick_up("DATA"),
            x_stab=pick_up("X-STAB", "X-STAB-BOUND-A-C", "X-STAB-BOUND-B-C"),
            z_stab=pick_up("Z-STAB", "Z-STAB-BOUND-L-C", "Z-STAB-BOUND-R-C"),
        )


@dataclass
class Patch_Target:
    coords: Dict[Coord, Label]
    data: List[Index]
    x_stab: List[Index]
    z_stab: List[Index]

    @classmethod
    def from_coords(
        cls,
        coords: Dict[Coord, Label],
        q2i: Mapping[Coord, Index],
    ) -> "Patch_Target":
        def pick_up(*labels: str) -> List[Index]:
            return _pick_up_indices(coords, q2i, *labels)

        return cls(
            coords=coords,
            data=pick_up("DATA"),
            x_stab=pick_up("X-STAB", "X-STAB-BOUND-A-T", "X-STAB-BOUND-B-T"),
            z_stab=pick_up("Z-STAB", "Z-STAB-BOUND-L-T", "Z-STAB-BOUND-R-T"),
        )


@dataclass
class Patch_Surgery:
    coords: Dict[Coord, Label]
    x_stab_m: List[Index]
    z_stab_m: List[Index]
    x_stab_b: List[Index]
    z_stab_l: List[Index]

    @classmethod
    def from_coords(
        cls,
        coords: Dict[Coord, Label],
        q2i: Mapping[Coord, Index],
    ) -> "Patch_Surgery":
        def pick_up(*labels: str) -> List[Index]:
            return _pick_up_indices(coords, q2i, *labels)

        return cls(
            coords=coords,
            x_stab_m=pick_up("X-STAB-SURGERY-M"),
            z_stab_m=pick_up("Z-STAB-SURGERY-M"),
            x_stab_b=pick_up("X-STAB-SURGERY-B"),
            z_stab_l=pick_up("Z-STAB-SURGERY-L"),
        )


#--------------------------
# Unified Contexts
#--------------------------

@dataclass
class Context:
    """
    Shared information across surface-style lattices
    
    -> Extra mappings are optional and default to empty dicts so code can selectively use the pieces it needs.
    """
    q2i: Dict[Coord, Index]
    i2q: Dict[Index, Coord]
    stab_to_data: Dict[Pair, str]
    stab_to_data_modified: Dict[Pair, str] = field(default_factory=dict)
    stab_to_data_modified2: Dict[Pair, str] = field(default_factory=dict)
    stab_to_data_modified3: Dict[Pair, str] = field(default_factory=dict)


@dataclass
class LatticeContext:
    """
    Shared information across lattice-surgery-style lattices
    """
    q2i: Dict[Coord, Index]
    i2q: Dict[Index, Coord]
    stab_to_data: Dict[Pair, str]
    stab_to_data_surgery_ac: Dict[Pair, str]
    stab_to_data_surgery_at: Dict[Pair, str]
    surgery_coords: Dict[Coord, Label]


# Backwards-compatibility helpers (thin aliases)
# These let you import the familiar names from this core module and then
# update the individual code modules to reference these symbols.
ConfigSurface = Config_Surface
ConfigLatticeSurgery = Config_LatticeSurgery

__all__ = [
    "Coord",
    "Label",
    "Index",
    "Pair",
    "NoiseModel",
    "CircuitResult",
    "ConfigSurface",
    "ConfigLatticeSurgery",
    "Patch",
    "Patch_Ancilla",
    "Patch_Control",
    "Patch_Target",
    "Patch_Surgery",
    "Context",
    "LatticeContext",
]
