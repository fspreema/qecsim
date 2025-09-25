Coord = complex
Label = str

__all__ = ["populate_stab_to_data"]

# ----------------------------
# Public Function
# ----------------------------

def populate_stab_to_data(patch: dict[Coord, Label], 
                          is_flipped: bool = False, 
                          y_basis: bool = False, 
                          y_switch: bool = False,
                          y_memory: bool = False, 
                          distance: int = 0
                          ) -> dict[tuple[Coord, Coord], str] | tuple[dict[tuple[Coord, Coord], str], dict[tuple[Coord, Coord], str]]:
    """
    Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

    * X-Stabs: control ist the **Data** qubit -> (data, stab)
    * Z-Stabs: control is the **stab** qubit -> (stab, data)

    * is_fipped -> Switching the role of X and Z stabilizers
    * y_basis -> Initlizing one z edge and one x edge. Each edge is build up out of two sides of the lattice
    * y_switch -> CX schedule after switch (H & SQRT X DEG gates) -> Implementation of XCY Schedule 
    """

    if not y_switch:
        stab_to_data: dict[tuple[Coord, Coord], Label] = {}
        _attach_interior_cx(patch, stab_to_data, is_flipped, y_switch, y_basis, distance)
        _attach_boundary_cx(patch, stab_to_data, is_flipped, y_switch, y_basis, y_memory, distance)

        return stab_to_data

    if y_switch:
        if not y_memory:
            stab_to_data: dict[tuple[Coord, Coord], Label] = {}
            stab_to_data_xcy: dict[tuple[Coord, Coord], Label] = {}
            _attach_interior_cx(patch, stab_to_data, is_flipped, y_switch, y_basis, distance, stab_to_data_xcy)
            _attach_boundary_cx(patch, stab_to_data, is_flipped, y_switch, y_basis, y_memory, distance)
        else:
            stab_to_data: dict[tuple[Coord, Coord], Label] = {}
            _attach_interior_cx(patch, stab_to_data, is_flipped, y_switch, y_basis, distance)
            _attach_boundary_cx(patch, stab_to_data, is_flipped, y_switch, y_basis, y_memory, distance)

        return stab_to_data, stab_to_data_xcy

# ------------------------------
# Internal helper function
# ------------------------------

def _neighbours(c: complex, dx: int, dy: int) -> complex:
    """Return neighbor coord at offsets (dx, dy)."""
    return (c.real + dx) + (c.imag + dy) * 1j

def _quad(c: complex) -> tuple[complex, complex, complex, complex]:
    """
    Return (E-Down, W-Down, E-Up, W-Up) neighbors of a stabilizer.
    -> Only the normal 4 weight stabilizers (excluding the Y-Case and Boundary)
    """
    return (
        _neighbours(c, +1, -1),  # q1
        _neighbours(c, -1, -1),  # q2
        _neighbours(c, +1, +1),  # q3
        _neighbours(c, -1, +1),  # q4
    )

def _assign_orders(
    table: dict[tuple[complex, complex], str],
    pairs: list[tuple[complex, complex]],
    orders: list[str],
    ) -> None:

    """Write (key -> order) entries into schedule dict."""
    for (a, b), order in zip(pairs, orders):
        table[(a, b)] = order


def _attach_interior_cx(patch: dict[Coord, Label], 
                        stab_to_data: dict[tuple[Coord, Coord], str], 
                        flipped: bool, 
                        y_switch: bool, 
                        y_basis: bool, 
                        distance: int, 
                        stab_to_data_xcy: dict[tuple[Coord, Coord], str] = {}):
    """
    Adds the 4-body CX Schedule for the *interior* stabilizers
    """

    # Adding global Order Lists -> List of qubit paris stays the same i.e. (q1, q1, q3, q4)
    ORDERS_X_NORMAL = ["1-CX", "2-CX", "3-CX", "4-CX"]
    ORDERS_X_YBASIS = ["4-CX", "3-CX", "2-CX", "1-CX"]
    ORDERS_Z_NORMAL = ["1-CX", "3-CX", "2-CX", "4-CX"]
    ORDERS_Z_YBASIS = ["4-CX", "2-CX", "3-CX", "1-CX"]

    if not flipped:
        if not y_switch:
            if not y_basis:

                for coords, qtype in patch.items():

                    # Defining Needed Neighbour Qubits
                    q1, q2, q3, q4 = _quad(coords)

                    if qtype == "X-STAB":
                        # Control is data -> (data, stab)
                        pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                        _assign_orders(stab_to_data, pairs, ORDERS_X_NORMAL)

                    elif qtype == "Z-STAB":
                        # Control is stab -> (stab, data)
                        pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                        _assign_orders(stab_to_data, pairs, ORDERS_Z_NORMAL)

            else:

                for coords, qtype in patch.items():

                    # Defining Needed Neighbour Qubits
                    q1, q2, q3, q4 = _quad(coords)

                    #Using Y Basis Order
                    if qtype == "X-STAB":
                        # Control is data -> (data, stab)
                        pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                        _assign_orders(stab_to_data, pairs, ORDERS_X_YBASIS)
                    elif qtype == "Z-STAB":
                        # Control is stab -> (stab, data)
                        pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                        _assign_orders(stab_to_data, pairs, ORDERS_Z_YBASIS)


        elif y_switch:

            # Filtering out the stabs needed for the two XCY Gate TICKS
            filtered_stabs_x : list[complex] = []
            filtered_stabs_z : list[complex] = []

            for i in range(distance - 1):
                stabs_z = 2 + (i * 2) + 2j + (i * 2) * 1j
                stabs_x = 4 + (i * 2) + 2j + (i * 2) * 1j
                filtered_stabs_z.append(stabs_z)
                filtered_stabs_x.append(stabs_x)

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
            for coords,qtype in patch.items():

                if coords in filtered_stabs_x:
                    new_cord = _neighbours(coords, -1, +1)
                    stab_to_data_xcy[coords, new_cord] = "1TICK"


            for coords, qtype in stabs_norm_dict.items():

                # Defining Needed Neighbour Qubits
                q1, q2, q3, q4 = _quad(coords)

                if qtype == "X-STAB":

                    pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                    _assign_orders(stab_to_data, pairs, ["2TICK", "3TICK", "4TICK", "5TICK"])

                elif qtype == "Z-STAB":
                    
                    # Checking whether normal CX or the XCY gate
                    if coords not in filtered_stabs_z:
                        pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                        _assign_orders(stab_to_data, pairs, ["2TICK", "4TICK", "3TICK", "5TICK"])

                    else:
                        _assign_orders(stab_to_data, [(q3, coords)], ["3.5TICK"])
                        _assign_orders(stab_to_data, [(coords, q2), (coords, q4)], ["4TICK", "5TICK"])

            """
            I HAVE NO CLUE WHY ONLY WEIGHT 3 instead of weight 4
            """

            # Implementing regular CX scheduele on H half -> only weight 3
            for coords,qtype in stabs_h_dict.items():
                #Implementing orientation of CX with sub schedule of XCY Gates
                if qtype == "X-STAB":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j

                    #Exclude the x stabs next to the diagonal:
                    if coords not in filtered_stabs_x:
                        stab_to_data[coords, new_cord1] = "2TICK"

                    stab_to_data[coords, new_cord2] = "3TICK"
                    stab_to_data[new_cord3, coords] = "4TICK"
                    stab_to_data[new_cord1, coords] = "5TICK"


                elif qtype == "Z-STAB":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
                    
                    stab_to_data[new_cord2, coords] = "3TICK"
                    stab_to_data[coords, new_cord3] = "4TICK"
                    stab_to_data[new_cord1, coords] = "2TICK"
                    # HOTFIX: DATASET WILL OVERWRITE 2TICK if they are the same
                    stab_to_data[coords , new_cord1, 2] = "5TICK"

    else:

        # Switching the Stabilizer roles -> Z stab has now the X Stab routine and vice versa
        for coords, qtype in patch.items():

            # Defining Needed Neighbour Qubits
            q1, q2, q3, q4 = _quad(coords)

            if qtype == "X-STAB":
                # Control is stab -> (stab, data)
                pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                _assign_orders(stab_to_data, pairs, ORDERS_Z_NORMAL)

            elif qtype == "Z-STAB":
                # Control is data -> (data, stab)
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                _assign_orders(stab_to_data, pairs, ORDERS_X_NORMAL)

def _assign_boundary_case(
    table: dict[tuple[complex, complex], str],
    label: str,
    coords: complex,
    case: dict,
    ) -> None:
    
    targets = [_neighbours(coords, dx, dy) for dx, dy in case["offsets"]]
    if case["pair"] == "stab->data":
        pairs = [(coords, t) for t in targets]
    else:
        pairs = [(t, coords) for t in targets]
    _assign_orders(table, pairs, case["orders"])


def _attach_boundary_cx(patch: dict[Coord, Label], 
                        stab_to_data: dict[tuple[Coord, Coord], str], 
                        flipped: bool, 
                        y_switch: bool,
                        y_basis: bool, 
                        y_memory: bool, 
                        distance : int):
    """
    Adds the 2-body CX Schedule for the *boundary* stabilizers
    """

    BOUNDARY_CASES_NORMAL = {
    "Z-STAB-BOUND-L": {"offsets": [(+1, -1), (+1, +1)], "orders": ["1-CX", "2-CX"], "pair": "stab->data"},
    "Z-STAB-BOUND-R": {"offsets": [(-1, -1), (-1, +1)], "orders": ["3-CX", "4-CX"], "pair": "stab->data"},
    "X-STAB-BOUND-U": {"offsets": [(-1, +1), (+1, +1)], "orders": ["4-CX", "3-CX"], "pair": "data->stab"},
    "X-STAB-BOUND-B": {"offsets": [(-1, -1), (+1, -1)], "orders": ["2-CX", "1-CX"], "pair": "data->stab"},
    }

    if not flipped:
        if not y_basis:

            for coords, qtype in patch.items():

                if qtype == "Z-STAB-BOUND-L":
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "1-CX"
                    stab_to_data[coords, new_cord2] = "2-CX"

                elif qtype == "Z-STAB-BOUND-R":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "3-CX"
                    stab_to_data[coords, new_cord2] = "4-CX"
                    
                elif qtype == "X-STAB-BOUND-U":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "4-CX"
                    stab_to_data[new_cord2, coords] = "3-CX"

                elif qtype == "X-STAB-BOUND-B":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                    stab_to_data[new_cord1, coords] = "2-CX"
                    stab_to_data[new_cord2, coords] = "1-CX"

        elif y_basis:

            if not y_switch and not y_memory:

                for coords, qtype in patch.items():

                    if qtype == "Z-STAB-BOUND-L":
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[coords, new_cord1] = "4-CX"
                        stab_to_data[coords, new_cord2] = "3-CX"

                    elif qtype == "X-STAB-BOUND-R":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[new_cord1, coords] = "3-CX"
                        stab_to_data[new_cord2, coords] = "1-CX"
                        
                    elif qtype == "Z-STAB-BOUND-U":
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[coords, new_cord1] = "1-CX"
                        stab_to_data[coords, new_cord2] = "3-CX"

                    elif qtype == "X-STAB-BOUND-B":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                        stab_to_data[new_cord1, coords] = "3-CX"
                        stab_to_data[new_cord2, coords] = "4-CX"

            elif y_switch and not y_memory:

                for coords, qtype in patch.items():
                    """
                    As the H gates was applied we need to flip the corressponding stabilizer schedule
                    
                    * ATTENTION -> ONLY FLIP THESE WHICH WHERE FLIPPED I.E. on only one diagonal half
                                -> Additional Boundary Stabs do CX both ways i.e. detecting z and x errors!
                    """
                    if qtype == "Z-STAB-BOUND-L":
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[coords, new_cord1] = "2TICK"
                        stab_to_data[coords, new_cord2] = "3TICK"

                    elif qtype == "X-STAB-BOUND-R":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j

                        # Check for lower boundary condition and exclude the cx which gets replaced by CYX
                        lower_boundary = [i for i in range(4, distance * 2, 4)][-1]
                        lower_coord = (distance * 2) + lower_boundary * 1j

                        if coords != lower_coord:
                            stab_to_data[coords, new_cord2,] = "2TICK"
                            
                        stab_to_data[coords, new_cord1] = "3TICK"

                    elif qtype == "X-STAB-BOUND-R-H":
                        # H Boundary has mixed cx direction
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[coords, new_cord1] = "4TICK"
                        stab_to_data[new_cord2, coords] = "2TICK"
                        stab_to_data[coords, new_cord2] = "5TICK"

                    elif qtype == "Z-STAB-BOUND-U-H":
                        # H Boundary has mixed direction
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j

                        if new_cord1 != 1 + 1j:
                            stab_to_data[coords, new_cord1] = "2TICK"

                        stab_to_data[new_cord2, coords] = "4TICK"
                        stab_to_data[new_cord1, coords] = "5TICK"

                        
                    elif qtype == "Z-STAB-BOUND-U":
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[new_cord1, coords] = "2TICK"
                        stab_to_data[new_cord2, coords] = "3TICK"

                    elif qtype == "X-STAB-BOUND-B":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                        stab_to_data[new_cord1, coords] = "3TICK"
                        stab_to_data[new_cord2, coords] = "2TICK"

            elif y_memory:

                for coords, qtype in patch.items():
                    """
                    In the memory Round we switch to the newly introduced boundary operators and deactivate the old ones
                    """

                    if qtype == "Z-STAB-BOUND-L":
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[coords, new_cord1] = "4-CX"
                        stab_to_data[coords, new_cord2] = "3-CX"


                    elif qtype == "X-STAB-BOUND-R-H":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[new_cord1, coords] = "2-CX"
                        stab_to_data[new_cord2, coords] = "1-CX"

                    elif qtype == "Z-STAB-BOUND-U-H":
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                        stab_to_data[new_cord1, coords] = "1-CX"
                        stab_to_data[new_cord2, coords] = "2-CX"


                    elif qtype == "X-STAB-BOUND-B":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                        stab_to_data[new_cord1, coords] = "3-CX"
                        stab_to_data[new_cord2, coords] = "4-CX"
                    

    else:

        for coords, qtype in patch.items():

            # Switching the stabilizers so the Z-Stab have switched CX and therefore act like a X-Stab and vice versa

            if qtype == "Z-STAB-BOUND-L":
                new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                stab_to_data[new_cord1, coords] = "2-CX"
                stab_to_data[new_cord2, coords] = "1-CX"

            elif qtype == "Z-STAB-BOUND-R":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
                stab_to_data[new_cord1, coords] = "4-CX"
                stab_to_data[new_cord2, coords] = "3-CX"
                
            elif qtype == "X-STAB-BOUND-U":
                new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
                stab_to_data[coords, new_cord1] = "3-CX"
                stab_to_data[coords, new_cord2] = "4-CX"

            elif qtype == "X-STAB-BOUND-B":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
                stab_to_data[coords, new_cord1] = "1-CX"
                stab_to_data[coords, new_cord2] = "2-CX"