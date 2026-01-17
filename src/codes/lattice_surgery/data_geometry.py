from dataclasses import dataclass

from src.codes.lattice_surgery.get_stab_pairings import LatticeSurgeryPairings


@dataclass
class MasterPairings:
    std_pairings: LatticeSurgeryPairings
    ac_merge_pairings: LatticeSurgeryPairings
    at_merge_pairings: LatticeSurgeryPairings
