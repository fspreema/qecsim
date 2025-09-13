from typing import Dict, Tuple

Coord = complex
Label = str

__all__ = ["populate_stab_to_data"]

# ----------------------------
# Public Function
# ----------------------------

def populate_stab_to_data(patch: Dict[Coord, Label], is_flipped : bool) -> Dict[Tuple[Coord, Coord], str]:
    """
    Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

    * X-Stabs: control ist the **Data** qubit -> (data, stab)
    * Z-Stabs: control is the **stab** qubit -> (stab, data)
    """

    stab_to_data: Dict[Tuple[Coord, Coord], Label] = {}

    _attach_interior_cx(patch, stab_to_data, is_flipped)
    _attach_boundary_cx(patch, stab_to_data, is_flipped)

    return stab_to_data


# ------------------------------
# Internal helper function
# ------------------------------

def _attach_interior_cx(patch: Dict[Coord, Label], stab_to_data: Dict[Tuple[Coord, Coord], str], flipped : bool):
    """
    Adds the 4-body CX Schedule for the *interior* stabilizers
    """

    if flipped == False:

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

    else:

        # Switching the Stabilizer roles -> Z stab has now the X Stab routine and vice versa

        for coords,string in patch.items():
            #Already in Correct Orientation for Measurement of CX
            if string == "Z-STAB":
                new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
                new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                stab_to_data[new_cord1, coords] = "1-CX"
                stab_to_data[new_cord2, coords] = "2-CX"
                stab_to_data[new_cord3, coords] = "3-CX"
                stab_to_data[new_cord4, coords] = "4-CX"

            elif string == "X-STAB":
                new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                stab_to_data[coords, new_cord1] = "1-CX"
                stab_to_data[coords, new_cord2] = "2-CX"
                stab_to_data[coords, new_cord3] = "3-CX"
                stab_to_data[coords, new_cord4] = "4-CX"

def _attach_boundary_cx(patch: Dict[Coord, Label], stab_to_data: Dict[Tuple[Coord, Coord], str], flipped : bool):
    """
    Adds the 2-body CX Schedule for the *boundary* stabilizers
    """

    if flipped == False:

        for coords,string in patch.items():

            if string == "Z-STAB-BOUND-L":
                new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                stab_to_data[coords, new_cord1] = "1-CX"
                stab_to_data[coords, new_cord2] = "2-CX"

            elif string == "Z-STAB-BOUND-R":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
                stab_to_data[coords, new_cord1] = "3-CX"
                stab_to_data[coords, new_cord2] = "4-CX"
                
            elif string == "X-STAB-BOUND-U":
                new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                stab_to_data[new_cord1, coords] = "4-CX"
                stab_to_data[new_cord2, coords] = "3-CX"

            elif string == "X-STAB-BOUND-B":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                stab_to_data[new_cord1, coords] = "2-CX"
                stab_to_data[new_cord2, coords] = "1-CX"

    else:

        for coords,string in patch.items():

            # Switching the stabilizers so the Z-Stab have switched CX and therefore act like a X-Stab and vice versa

            if string == "Z-STAB-BOUND-L":
                new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                stab_to_data[new_cord1, coords] = "2-CX"
                stab_to_data[new_cord2, coords] = "1-CX"

            elif string == "Z-STAB-BOUND-R":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
                stab_to_data[new_cord1, coords] = "4-CX"
                stab_to_data[new_cord2, coords] = "3-CX"
                
            elif string == "X-STAB-BOUND-U":
                new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                stab_to_data[coords, new_cord1] = "4-CX"
                stab_to_data[coords, new_cord2] = "3-CX"

            elif string == "X-STAB-BOUND-B":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                stab_to_data[coords, new_cord1] = "2-CX"
                stab_to_data[coords, new_cord2] = "1-CX"