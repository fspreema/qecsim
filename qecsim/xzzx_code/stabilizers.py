from typing import Dict, Tuple

Coord = complex
Label = str

__all__ = ["populate_stab_to_data"]

# ----------------------------
# Public Function
# ----------------------------

def populate_stab_to_data(patch: Dict[Coord, Label]) -> Dict[Tuple[Coord, Coord], str]:
    """
    Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice
    """

    stab_to_data: Dict[Tuple[Coord, Coord], Label] = {}

    _attach_interior_cx(patch, stab_to_data)
    _attach_boundary_cx(patch, stab_to_data)

    return stab_to_data

# ------------------------------
# Internal helper function
# ------------------------------

def _attach_interior_cx(patch: Dict[Coord, Label], stab_to_data: Dict[Tuple[Coord, Coord], str]):
    """
    Adds the 4-body CX Schedule for the *interior* stabilizers
    """

    for coords,string in patch.items():
        #Already in Correct Orientation for Measurement of CX
        if string == "STAB-Ver":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord3 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_cord4 = (coords.real + 1) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "1-CX"
            stab_to_data[new_cord2, coords] = "2-CZ"
            stab_to_data[new_cord3, coords] = "3-CZ"
            stab_to_data[new_cord4, coords] = "4-CX"
            
        elif string == "STAB-Hor":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord3 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_cord4 = (coords.real + 1) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "1-CX"
            stab_to_data[new_cord2, coords] = "2-CZ"
            stab_to_data[new_cord3, coords] = "3-CZ"
            stab_to_data[new_cord4, coords] = "4-CX"

def _attach_boundary_cx(patch: Dict[Coord, Label], stab_to_data: Dict[Tuple[Coord, Coord], str]):
    """
    Adds the 2-body CX Schedule for the *boundary* stabilizers
    """

    for coords,string in patch.items():

        if string == "STAB-BOUND-L-Hor":
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "3-CZ"
            stab_to_data[new_cord2, coords] = "4-CX"
  
        elif string == "STAB-BOUND-R-Hor":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "1-CX"
            stab_to_data[new_cord2, coords] = "2-CZ"
               
        elif string == "STAB-BOUND-A-Ver":
            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "3-CZ"
            stab_to_data[new_cord2, coords] = "4-CX"

        elif string == "STAB-BOUND-B-Ver":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
            stab_to_data[new_cord1, coords] = "1-CX"
            stab_to_data[new_cord2, coords] = "2-CZ" 
