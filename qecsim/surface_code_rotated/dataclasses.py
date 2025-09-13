from dataclasses import dataclass
from typing import Dict, List, Tuple, Mapping

Coord = complex
Label = str
Index = int
Pair  = Tuple[Coord, Coord]

#--------------------------
# Helper Dataclasses
#--------------------------

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
    z_stab:  List[Index]

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
            x_stab = pick_up("X-STAB", "X-STAB-BOUND-U", "X-STAB-BOUND-B"),
            z_stab = pick_up("Z-STAB", "Z-STAB-BOUND-L", "Z-STAB-BOUND-R"),
        )

@dataclass
class Context:
    """Shared information across all lattices -> used by all helper functions *lct*"""
    q2i: Dict[Coord, Index]
    i2q: Dict[Index, Coord]
    stab_to_data: Dict[Pair, str]
    stab_to_data_flipped: Dict[Pair, str]