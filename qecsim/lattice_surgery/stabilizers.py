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

    * X-Stabs: control ist the **Data** qubit -> (data, stab)
    * Z-Stabs: control is the **stab** qubit -> (stab, data)
    """

    stab_to_data: Dict[Tuple[Coord, Coord], str] = {}

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
        if string == "X-STAB":
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
            new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "1-CX"
            stab_to_data[new_cord2, coords] = "2-CX"
            stab_to_data[new_cord3, coords] = "3-CX"
            stab_to_data[new_cord4, coords] = "4-CX"

        elif string == "Z-STAB":
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
            new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "1-CX"
            stab_to_data[coords, new_cord2] = "2-CX"
            stab_to_data[coords, new_cord3] = "3-CX"
            stab_to_data[coords, new_cord4] = "4-CX"

def _attach_boundary_cx(patch: Dict[Coord, Label], stab_to_data: Dict[Tuple[Coord, Coord], str]):
    """
    Adds the 2-body CX Schedule for the *boundary* stabilizers
    """
    for coords,string in patch.items():

        if string == "Z-STAB-BOUND-L-A":
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "1-CX"
            stab_to_data[coords, new_cord2] = "2-CX"

        elif string == "Z-STAB-BOUND-R-A":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "3-CX"
            stab_to_data[coords, new_cord2] = "4-CX"
        
        elif string == "X-STAB-BOUND-A-A":
            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "4-CX"
            stab_to_data[new_cord2, coords] = "3-CX"

        elif string == "X-STAB-BOUND-B-A":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
            stab_to_data[new_cord1, coords] = "2-CX"
            stab_to_data[new_cord2, coords] = "1-CX"

        elif string == "Z-STAB-BOUND-L-T":
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "5-CX"
            stab_to_data[coords, new_cord2] = "6-CX"

        elif string == "Z-STAB-BOUND-R-T":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "5-CX"
            stab_to_data[coords, new_cord2] = "6-CX"
        
        elif string == "X-STAB-BOUND-A-T":
            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "5-CX"
            stab_to_data[new_cord2, coords] = "6-CX"

        elif string == "X-STAB-BOUND-B-T":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
            stab_to_data[new_cord1, coords] = "5-CX"
            stab_to_data[new_cord2, coords] = "6-CX"

        elif string == "Z-STAB-BOUND-L-C":
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "5-CX"
            stab_to_data[coords, new_cord2] = "6-CX"

        elif string == "Z-STAB-BOUND-R-C":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "5-CX"
            stab_to_data[coords, new_cord2] = "6-CX"
        
        elif string == "X-STAB-BOUND-A-C":
            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "5-CX"
            stab_to_data[new_cord2, coords] = "6-CX"

        elif string == "X-STAB-BOUND-B-C":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
            stab_to_data[new_cord1, coords] = "5-CX"
            stab_to_data[new_cord2, coords] = "6-CX"