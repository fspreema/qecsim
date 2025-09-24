from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Mapping, Optional

import stim

Coord = complex
Label = str
Index = int
Pair  = Tuple[Coord, Coord]

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
    obs_indices: Optional[list[int]] = None

    def __iadd__(self, other: "CircuitResult | stim.Circuit") -> "CircuitResult":
        """Allow `result += other` regardless of other’s type."""
        if isinstance(other, CircuitResult):
            self.circuit += other.circuit
            # No automatic merge of obs indices; caller decides how to use them
        else:
            self.circuit += other
        return self

@dataclass
class Config:
    """Non tunable init setting get passed as *cfg*"""
    distance : int
    state_init : str
    obs : str
    rounds : int

@dataclass
class Patch:
    """All geometry information is stored here i.e. data x_stab indicees"""
    coords:  Dict[Coord, Label]
    data:    List[Index]
    x_stab:  List[Index]
    stab_switch_apply_h : List[Index]
    z_stab:  List[Index]
    upper_h:  List[Index]
    right_h:  List[Index]

    #Create Factory to build up the Lists from the coords and the q2i dict

    @classmethod
    def from_coords(
        cls,
        coords: Dict[Coord, Label],
        q2i : Mapping[Coord, Index],
    ) -> "Patch":
        
        def pick_up(*labels: str) -> List[Index]:
            return [q2i[q] for q, t in coords.items() if t in labels]
        
        return cls(
            coords = coords,
            data   = pick_up("DATA"),
            x_stab = pick_up("X-STAB", "X-STAB-BOUND-U", "X-STAB-BOUND-B", "X-STAB-BOUND-R"),
            stab_switch_apply_h = pick_up("X-STAB", "Z-STAB-BOUND-U", "X-STAB-BOUND-B", "Z-STAB-BOUND-U-H"),
            z_stab = pick_up("Z-STAB", "Z-STAB-BOUND-L", "Z-STAB-BOUND-R", "Z-STAB-BOUND-U"),
            upper_h = pick_up("Z-STAB-BOUND-U-H"),
            right_h = pick_up("X-STAB-BOUND-R-H"),
        )

@dataclass
class Context:
    """Shared information across all lattices -> used by all helper functions *lct*"""
    q2i: Dict[Coord, Index]
    i2q: Dict[Index, Coord]
    stab_to_data: Dict[Pair, str]
    stab_to_data_modified: Optional[Dict[Pair, str]] = field(default_factory=dict)
    stab_to_data_modified2: Optional[Dict[Pair, str]] = field(default_factory=dict)