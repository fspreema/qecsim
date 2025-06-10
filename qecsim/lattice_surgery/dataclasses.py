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
    target_state_init: str
    control_state_init: str

@dataclass
class Patch_Ancilla:
    """All geometry information is stored here i.e. data x_stab indicees"""
    coords:  Dict[Coord, Label]
    data:    List[Index]
    x_stab:  List[Index]
    z_stab:  List[Index]
    x_bdyB:  List[Index]
    z_bdyR:  List[Index]

    #Create Factory to build up the Lists from the coords and the q2i dict

    @classmethod
    def from_coords(
        cls,
        coords: Dict[Coord, Label],
        q2i : Mapping[Coord, Index],
    ) -> "Patch_Ancilla":
        
        def pick_up(*labels: str) -> List[Index]:
            return [q2i[q] for q, t in coords.items() if t in labels]
        
        return cls(
            coords = coords,
            data   = pick_up("DATA"),
            x_stab = pick_up("X-STAB", "X-STAB-BOUND-A-A", "X-STAB-BOUND-B-A"),
            z_stab = pick_up("Z-STAB", "Z-STAB-BOUND-L-A", "Z-STAB-BOUND-R-A"),
            x_bdyB = pick_up("X-STAB-BOUND-B-A"),
            z_bdyR = pick_up("Z-STAB-BOUND-R-A"),
        )
    
@dataclass
class Patch_Control:
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
    ) -> "Patch_Control":
        
        def pick_up(*labels: str) -> List[Index]:
            return [q2i[q] for q, t in coords.items() if t in labels]
        
        return cls(
            coords = coords,
            data   = pick_up("DATA"),
            x_stab = pick_up("X-STAB", "X-STAB-BOUND-A-C", "X-STAB-BOUND-B-C"),
            z_stab = pick_up("Z-STAB", "Z-STAB-BOUND-L-C", "Z-STAB-BOUND-R-C"),
        )
    
@dataclass
class Patch_Target:
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
    ) -> "Patch_Target":
        
        def pick_up(*labels: str) -> List[Index]:
            return [q2i[q] for q, t in coords.items() if t in labels]
        
        return cls(
            coords = coords,
            data   = pick_up("DATA"),
            x_stab = pick_up("X-STAB", "X-STAB-BOUND-A-T", "X-STAB-BOUND-B-T"),
            z_stab = pick_up("Z-STAB", "Z-STAB-BOUND-L-T", "Z-STAB-BOUND-R-T"),
        )
    
@dataclass
class Patch_Surgery:
    """All geometry information is stored here i.e. data x_stab indicees"""
    coords:  Dict[Coord, Label]
    x_stab_m:  List[Index]
    z_stab_m:  List[Index]
    x_stab_b:  List[Index]
    z_stab_l:  List[Index]

    #Create Factory to build up the Lists from the coords and the q2i dict

    @classmethod
    def from_coords(
        cls,
        coords: Dict[Coord, Label],
        q2i : Mapping[Coord, Index],
    ) -> "Patch_Surgery":
        
        def pick_up(*labels: str) -> List[Index]:
            return [q2i[q] for q, t in coords.items() if t in labels]
        
        return cls(
            coords = coords,
            x_stab_m = pick_up("X-STAB-SURGERY-M"),
            z_stab_m = pick_up("Z-STAB-SURGERY-M"),
            x_stab_b = pick_up("X-STAB-SURGERY-B"),
            z_stab_l = pick_up("Z-STAB-SURGERY-L")
        )

@dataclass
class LatticeContext:
    """Shared information across all lattices -> used by all helper functions *lct*"""
    q2i: Dict[Coord, Index]
    i2q: Dict[Index, Coord]
    stab_to_data: Dict[Pair, str]
    stab_to_data_surgery_ac: Dict[Pair, str]
    stab_to_data_surgery_at: Dict[Pair, str]
    surgery_coords: Dict[Coord, Label]