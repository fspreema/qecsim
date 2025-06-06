from typing import Dict, Tuple

Coord = complex
Label = str

__all__ = ["populate_stab_to_data"]

# ----------------------------
# Public Function
# ----------------------------

def populate_stab_to_data(patch: Dict[Coord, Label], merging : bool, merging_type : str = False) -> Dict[Tuple[Coord, Coord], str]:
    """
    Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

    * X-Stabs: control ist the **Data** qubit -> (data, stab)
    * Z-Stabs: control is the **stab** qubit -> (stab, data)
    """

    stab_to_data: Dict[Tuple[Coord, Coord], Label] = {}

    _attach_interior_cx(patch, stab_to_data, merging, merging_type)
    _attach_boundary_cx(patch, stab_to_data, merging, merging_type)

    return stab_to_data


# ------------------------------
# Internal helper function
# ------------------------------

def _attach_interior_cx(patch: Dict[Coord, Label], stab_to_data: Dict[Tuple[Coord, Coord], str], merging : bool, merging_type : str):
    """
    Adds the 4-body CX Schedule for the *interior* stabilizers
    """

    if merging == False:

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

    elif merging == True:
        
        if merging_type == "AC":

            for coords,string in patch.items():

                if string in {"X-STAB", "X-STAB-BOUND-A-C"}:
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
                    new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "1-CX"
                    stab_to_data[new_cord2, coords] = "2-CX"
                    stab_to_data[new_cord3, coords] = "3-CX"
                    stab_to_data[new_cord4, coords] = "4-CX"

                elif string in {"Z-STAB", "Z-STAB-SURGERY-M"}:
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "1-CX"
                    stab_to_data[coords, new_cord2] = "2-CX"
                    stab_to_data[coords, new_cord3] = "3-CX"
                    stab_to_data[coords, new_cord4] = "4-CX"

        elif merging_type == "AT":

            for coords,string in patch.items():

                if string in {"X-STAB", "X-STAB-SURGERY-M"}:
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
                    new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "1-CX"
                    stab_to_data[new_cord2, coords] = "2-CX"
                    stab_to_data[new_cord3, coords] = "3-CX"
                    stab_to_data[new_cord4, coords] = "4-CX"

                elif string in {"Z-STAB", "Z-STAB-BOUND-L-T"}:
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "1-CX"
                    stab_to_data[coords, new_cord2] = "2-CX"
                    stab_to_data[coords, new_cord3] = "3-CX"
                    stab_to_data[coords, new_cord4] = "4-CX"

def _attach_boundary_cx(patch: Dict[Coord, Label], stab_to_data: Dict[Tuple[Coord, Coord], str], merging : bool, merging_type : str):
    """
    Adds the 2-body CX Schedule for the *boundary* stabilizers
    """

    if merging == False:

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

    elif merging == True:

        if merging_type == "AC":

            for coords,string in patch.items():

                if string in {"Z-STAB-BOUND-L-A", "Z-STAB-BOUND-L-C", "Z-STAB-SURGERY-L"}:
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "1-CX"
                    stab_to_data[coords, new_cord2] = "2-CX"

                elif string in {"Z-STAB-BOUND-R-A", "Z-STAB-BOUND-R-C"}:
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "3-CX"
                    stab_to_data[coords, new_cord2] = "4-CX"
                
                elif string == "X-STAB-BOUND-A-A":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "4-CX"
                    stab_to_data[new_cord2, coords] = "3-CX"

                elif string == "X-STAB-BOUND-B-C":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                    stab_to_data[new_cord1, coords] = "2-CX"
                    stab_to_data[new_cord2, coords] = "1-CX"

        elif merging_type == "AT":

            for coords,string in patch.items():

                if string == "Z-STAB-BOUND-L-A":
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "1-CX"
                    stab_to_data[coords, new_cord2] = "2-CX"

                elif string == "Z-STAB-BOUND-R-T":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "3-CX"
                    stab_to_data[coords, new_cord2] = "4-CX"
                
                elif string in {"X-STAB-BOUND-A-A", "X-STAB-BOUND-A-T"}:
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "4-CX"
                    stab_to_data[new_cord2, coords] = "3-CX"

                elif string in {"X-STAB-BOUND-B-A", "X-STAB-BOUND-B-T", "X-STAB-SURGERY-B"}:
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                    stab_to_data[new_cord1, coords] = "2-CX"
                    stab_to_data[new_cord2, coords] = "1-CX"