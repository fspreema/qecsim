from dataclasses import dataclass, field
from typing import Mapping
import stim

Coord = complex
Label = str
Index = int
Pair  = tuple[Coord, Coord]

#--------------------------
# Helper Dataclasses
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
    """Uniform return type for circuit-building functions."""
    circuit: stim.Circuit
    obs_indices: list[int] | None = None

    def __iadd__(self, other: "CircuitResult | stim.Circuit") -> "CircuitResult":
        """Allow `result += other` regardless of other’s type."""
        if isinstance(other, CircuitResult):
            self.circuit += other.circuit
            # No automatic merge of obs indices
        else:
            self.circuit += other
        return self

@dataclass
class Config:
    """Non tunable init setting get passed as *cfg*"""
    distance: int
    state_init: str
    obs: str
    rounds: int

@dataclass
class Patch:
    """All geometry information is stored here i.e. data/x_stab indices"""
    coords: dict[Coord, Label]
    data: list[Index]
    x_stab: list[Index]
    stab_switch_apply_h: list[Index]
    x_stab_memory: list[Index]
    z_stab_memory: list[Index]
    z_stab: list[Index]
    upper_h: list[Index]
    right_h: list[Index]

    @classmethod
    def from_coords(
        cls,
        coords: dict[Coord, Label],
        q2i: Mapping[Coord, Index],
    ) -> "Patch":
        def pick_up(*labels: str) -> list[Index]:
            return [q2i[q] for q, t in coords.items() if t in labels]

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

@dataclass
class Context:
    """Shared information across all lattices -> used by all helper functions *lct*"""
    q2i: dict[Coord, Index]
    i2q: dict[Index, Coord]
    stab_to_data: dict[Pair, str]
    stab_to_data_modified: dict[Pair, str] = field(default_factory=dict)
    stab_to_data_modified2: dict[Pair, str] = field(default_factory=dict)
    stab_to_data_modified3: dict[Pair, str] = field(default_factory=dict)


