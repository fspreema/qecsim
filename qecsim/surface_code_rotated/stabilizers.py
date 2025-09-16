from typing import Dict, Tuple

Coord = complex
Label = str

__all__ = ["populate_stab_to_data"]

# ----------------------------
# Public Function
# ----------------------------

def populate_stab_to_data(patch: Dict[Coord, Label], is_flipped : bool = False, y_basis : bool = False, 
                          y_switch : bool = False, distance : int = 0) -> Dict[Tuple[Coord, Coord], str]:
    """
    Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

    * X-Stabs: control ist the **Data** qubit -> (data, stab)
    * Z-Stabs: control is the **stab** qubit -> (stab, data)

    * is_fipped -> Switching the role of X and Z stabilizers
    * y_basis -> Initlizing one z edge and one x edge. Each edge is build up out of two sides of the lattice
    * y_switch -> CX schedule after switch (H & SQRT X DEG gates) -> Implementation of XCY Schedule 
    """

    if not y_switch:
        stab_to_data: Dict[Tuple[Coord, Coord], Label] = {}
        _attach_interior_cx(patch, stab_to_data, is_flipped, y_switch, distance)
        _attach_boundary_cx(patch, stab_to_data, is_flipped, y_basis, y_switch, distance)

        return stab_to_data

    else:
        stab_to_data: Dict[Tuple[Coord, Coord], Label] = {}
        stab_to_data_xcy: Dict[Tuple[Coord, Coord], Label] = {}
        _attach_interior_cx(patch, stab_to_data, is_flipped, y_switch, distance, stab_to_data_xcy)
        _attach_boundary_cx(patch, stab_to_data, is_flipped, y_basis, y_switch, distance)

        return stab_to_data, stab_to_data_xcy

# ------------------------------
# Internal helper function
# ------------------------------

def _attach_interior_cx(patch: Dict[Coord, Label], stab_to_data: Dict[Tuple[Coord, Coord], str], flipped : bool, 
                        y_switch : bool, distance : int, stab_to_data_xcy: Dict[Tuple[Coord, Coord], str] = {}):
    """
    Adds the 4-body CX Schedule for the *interior* stabilizers
    """

    if flipped == False:

        if not y_switch:

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

        elif y_switch:

            # Filtering out the stabs needed for the two XCY Gate TICKS
            filtered_stabs_x : list[complex] = []
            filtered_stabs_z : list[complex] = []

            for i in range(distance - 1):
                cords_z = 4 + (i * 2) + 2j + (i * 2) * 1j
                cords_x = 2 + (i * 2) + 2j + (i * 2) * 1j
                filtered_stabs_z.append(cords_z)
                filtered_stabs_x.append(cords_x)

            # Filtering out Stabilizers without H applied -> Normal CX-Schedule
            stabs_norm_dict : dict = {}
            stabs_h_dict : dict = {}
            stabs_xcy_dict = {}

            for cords, qtype in patch.items():
                # Diagonal Cut
                if cords.real <= cords.imag:
                    stabs_norm_dict[cords] = qtype
                else:
                    stabs_h_dict[cords] = qtype

            # Implementing Solo XCY Gate (Other one in the CX Schedule)
            for coords,string in patch.items():

                if coords in filtered_stabs_z:
                    new_cord = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data_xcy[new_cord, coords] = "1-XCY"


            for coords,string in stabs_norm_dict.items():
                if string == "X-STAB":

                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
                    new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j

                    # Checking whether normal CX or the XCY gate
                    if coords not in filtered_stabs_x:
                        stab_to_data[new_cord1, coords] = "1-CX"
                        stab_to_data[new_cord2, coords] = "2-CX"
                        stab_to_data[new_cord3, coords] = "3-CX"
                        stab_to_data[new_cord4, coords] = "4-CX"

                    else:
                        stab_to_data[new_cord2, coords] = "2-CX"
                        stab_to_data[coords, new_cord2] = "2-XCY"
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

            """
            I HAVE NO CLUE WHY ONLY WEIGHT 3 instead of weight 4
            """

            # Implementing regular CX scheduele on H half -> only weight 3
            for coords,string in stabs_h_dict.items():
                #Implementing orientation of CX with sub schedule of XCY Gates
                if string == "X-STAB":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j

                    stab_to_data[coords, new_cord1] = "1-CX"
                    stab_to_data[coords, new_cord2] = "2-CX"
                    stab_to_data[new_cord3, coords] = "3-CX"
                    stab_to_data[new_cord1, coords] = "4-CX"


                elif string == "Z-STAB":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
                    
                    #Exclude the z stabs next to the diagonal:
                    if coords not in filtered_stabs_z:
                        stab_to_data[new_cord1, coords] = "1-CX"

                    stab_to_data[new_cord2, coords] = "2-CX"
                    stab_to_data[coords, new_cord3] = "3-CX"
                    stab_to_data[coords, new_cord1] = "4-CX"

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

def _attach_boundary_cx(patch: Dict[Coord, Label], stab_to_data: Dict[Tuple[Coord, Coord], str], flipped : bool, 
                        y_basis : bool, y_switch : bool, distance : int):
    """
    Adds the 2-body CX Schedule for the *boundary* stabilizers
    """

    if flipped == False:

        if not y_basis:

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

        if y_basis:

            if not y_switch:

                for coords,string in patch.items():

                    if string == "X-STAB-BOUND-L":
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[new_cord1, coords] = "2-CX"
                        stab_to_data[new_cord2, coords] = "1-CX"

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

                    elif string == "Z-STAB-BOUND-B":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                        stab_to_data[coords, new_cord1] = "1-CX"
                        stab_to_data[coords,new_cord2] = "2-CX"

            else:

                for coords,string in patch.items():

                    # As the H gates was applied we need to flip the corressponding stabilizer schedule
                    # ATTENTION -> ONLY FLIP THESE WHICH WHERE FLIPPED I.E. on only one diagonal half

                    if string == "X-STAB-BOUND-L":
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[new_cord1, coords] = "B1-CX"
                        stab_to_data[new_cord2, coords] = "B2-CX"

                    elif string in {"Z-STAB-BOUND-R", "Z-STAB-BOUND-R-H"}:
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j

                        # Check for lower boundary condition and exclude the cx which gets replaced by CYX
                        lower_boundary = [i for i in range(4, distance * 2, 4)][-1]
                        lower_coord = (distance * 2) + lower_boundary * 1j

                        if coords != lower_coord:
                            stab_to_data[new_cord1, coords] = "B1-CX"
                            stab_to_data[new_cord2, coords] = "B2-CX"
                        else:
                            stab_to_data[new_cord1, coords] = "B1-CX"
                        
                    elif string in {"X-STAB-BOUND-U", "X-STAB-BOUND-U-H"}:
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[coords, new_cord1] = "B1-CX"
                        stab_to_data[coords, new_cord2] = "B2-CX"

                    elif string == "Z-STAB-BOUND-B":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                        stab_to_data[coords, new_cord1] = "B2-CX"
                        stab_to_data[coords, new_cord2] = "B1-CX"
                    

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