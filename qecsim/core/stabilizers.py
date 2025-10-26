Coord = complex
Label = str
Pair = tuple[Coord, Coord]

__all__ = ["populate_stab_to_data"]


def _populate_xzzx(patch: dict[Coord, Label]) -> dict[Pair, str]:
    """
    Populate schedule for XZZX-style stabilizers
    """

    stab_to_data: dict[Pair, str] = {}

    def _attach_interior_cx():
        """
        Adds the 4-body CX Schedule for the *interior* stabilizers
        """

        for coords, string in patch.items():
            if string == "STAB-Ver":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                new_cord3 = (coords.real - 1) + (coords.imag + 1) * 1j
                new_cord4 = (coords.real + 1) + (coords.imag + 1) * 1j
                stab_to_data[(new_cord1, coords)] = "1-CX"
                stab_to_data[(new_cord2, coords)] = "2-CZ"
                stab_to_data[(new_cord3, coords)] = "3-CZ"
                stab_to_data[(new_cord4, coords)] = "4-CX"
            elif string == "STAB-Hor":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                new_cord3 = (coords.real - 1) + (coords.imag + 1) * 1j
                new_cord4 = (coords.real + 1) + (coords.imag + 1) * 1j
                stab_to_data[(new_cord1, coords)] = "1-CX"
                stab_to_data[(new_cord2, coords)] = "2-CZ"
                stab_to_data[(new_cord3, coords)] = "3-CZ"
                stab_to_data[(new_cord4, coords)] = "4-CX"

    def _attach_boundary_cx():
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        for coords, string in patch.items():
            if string == "STAB-BOUND-L-Hor":
                new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                stab_to_data[(new_cord1, coords)] = "3-CZ"
                stab_to_data[(new_cord2, coords)] = "4-CX"
            elif string == "STAB-BOUND-R-Hor":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                stab_to_data[(new_cord1, coords)] = "1-CX"
                stab_to_data[(new_cord2, coords)] = "2-CZ"
            elif string == "STAB-BOUND-A-Ver":
                new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                stab_to_data[(new_cord1, coords)] = "3-CZ"
                stab_to_data[(new_cord2, coords)] = "4-CX"
            elif string == "STAB-BOUND-B-Ver":
                new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                stab_to_data[(new_cord1, coords)] = "1-CX"
                stab_to_data[(new_cord2, coords)] = "2-CZ"

    # Add interior and boundary CXs
    _attach_interior_cx()
    _attach_boundary_cx()
    return stab_to_data


def _populate_lattice_surgery(
    patch: dict[Coord, Label],
    merging: bool,
    merging_type: str | None = None,
) -> dict[Pair, str]:
    """
    Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

    * X-Stabs: control ist the **Data** qubit -> (data, stab)
    * Z-Stabs: control is the **stab** qubit -> (stab, data)
    """
    stab_to_data: dict[Pair, str] = {}

    def attach_interior(merging: bool, merging_type: str | None):
        """
        Adds the 4-body CX Schedule for the *interior* stabilizers
        """

        if not merging:
            for coords, string in patch.items():
                if string == "X-STAB":
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
                    new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[(new_cord1, coords)] = "1-CX"
                    stab_to_data[(new_cord2, coords)] = "2-CX"
                    stab_to_data[(new_cord3, coords)] = "3-CX"
                    stab_to_data[(new_cord4, coords)] = "4-CX"
                elif string == "Z-STAB":
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[(coords, new_cord1)] = "1-CX"
                    stab_to_data[(coords, new_cord2)] = "2-CX"
                    stab_to_data[(coords, new_cord3)] = "3-CX"
                    stab_to_data[(coords, new_cord4)] = "4-CX"
        else:
            if merging_type == "AC":
                for coords, string in patch.items():
                    if string in {"X-STAB", "X-STAB-BOUND-A-C"}:
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
                        new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                        stab_to_data[(new_cord1, coords)] = "1-CX"
                        stab_to_data[(new_cord2, coords)] = "2-CX"
                        stab_to_data[(new_cord3, coords)] = "3-CX"
                        stab_to_data[(new_cord4, coords)] = "4-CX"
                    elif string in {"Z-STAB", "Z-STAB-SURGERY-M"}:
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                        new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                        stab_to_data[(coords, new_cord1)] = "1-CX"
                        stab_to_data[(coords, new_cord2)] = "2-CX"
                        stab_to_data[(coords, new_cord3)] = "3-CX"
                        stab_to_data[(coords, new_cord4)] = "4-CX"
            elif merging_type == "AT":
                for coords, string in patch.items():
                    if string in {"X-STAB", "X-STAB-SURGERY-M"}:
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
                        new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                        stab_to_data[(new_cord1, coords)] = "1-CX"
                        stab_to_data[(new_cord2, coords)] = "2-CX"
                        stab_to_data[(new_cord3, coords)] = "3-CX"
                        stab_to_data[(new_cord4, coords)] = "4-CX"
                    elif string in {"Z-STAB", "Z-STAB-BOUND-L-T"}:
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                        new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                        stab_to_data[(coords, new_cord1)] = "1-CX"
                        stab_to_data[(coords, new_cord2)] = "2-CX"
                        stab_to_data[(coords, new_cord3)] = "3-CX"
                        stab_to_data[(coords, new_cord4)] = "4-CX"

    def attach_boundary(merging: bool, merging_type: str | None):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        if not merging:
            for coords, string in patch.items():
                if string == "Z-STAB-BOUND-L-A":
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "1-CX"
                    stab_to_data[coords, new_cord2] = "2-CX"

                elif string == "Z-STAB-BOUND-R-A":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "3-CX"
                    stab_to_data[coords, new_cord2] = "4-CX"

                elif string == "X-STAB-BOUND-A-A":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "4-CX"
                    stab_to_data[new_cord2, coords] = "3-CX"

                elif string == "X-STAB-BOUND-B-A":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                    stab_to_data[new_cord1, coords] = "2-CX"
                    stab_to_data[new_cord2, coords] = "1-CX"

                elif string == "Z-STAB-BOUND-L-T":
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "5-CX"
                    stab_to_data[coords, new_cord2] = "6-CX"

                elif string == "Z-STAB-BOUND-R-T":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "5-CX"
                    stab_to_data[coords, new_cord2] = "6-CX"

                elif string == "X-STAB-BOUND-A-T":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "5-CX"
                    stab_to_data[new_cord2, coords] = "6-CX"

                elif string == "X-STAB-BOUND-B-T":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                    stab_to_data[new_cord1, coords] = "5-CX"
                    stab_to_data[new_cord2, coords] = "6-CX"

                elif string == "Z-STAB-BOUND-L-C":
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "5-CX"
                    stab_to_data[coords, new_cord2] = "6-CX"

                elif string == "Z-STAB-BOUND-R-C":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "5-CX"
                    stab_to_data[coords, new_cord2] = "6-CX"

                elif string == "X-STAB-BOUND-A-C":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "5-CX"
                    stab_to_data[new_cord2, coords] = "6-CX"

                elif string == "X-STAB-BOUND-B-C":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                    stab_to_data[new_cord1, coords] = "5-CX"
                    stab_to_data[new_cord2, coords] = "6-CX"

        else:
            if merging_type == "AC":
                for coords, string in patch.items():
                    if string in {"Z-STAB-BOUND-L-A", "Z-STAB-BOUND-L-C", "Z-STAB-SURGERY-L"}:
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                        stab_to_data[(coords, new_cord1)] = "1-CX"
                        stab_to_data[(coords, new_cord2)] = "2-CX"

                    elif string in {"Z-STAB-BOUND-R-A", "Z-STAB-BOUND-R-C"}:
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                        stab_to_data[(coords, new_cord1)] = "3-CX"
                        stab_to_data[(coords, new_cord2)] = "4-CX"

                    elif string == "X-STAB-BOUND-A-A":
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                        stab_to_data[(new_cord1, coords)] = "4-CX"
                        stab_to_data[(new_cord2, coords)] = "3-CX"

                    elif string == "X-STAB-BOUND-B-C":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                        stab_to_data[(new_cord1, coords)] = "2-CX"
                        stab_to_data[(new_cord2, coords)] = "1-CX"

            elif merging_type == "AT":
                for coords, string in patch.items():
                    if string == "Z-STAB-BOUND-L-A":
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                        stab_to_data[(coords, new_cord1)] = "1-CX"
                        stab_to_data[(coords, new_cord2)] = "2-CX"

                    elif string == "Z-STAB-BOUND-R-T":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                        stab_to_data[(coords, new_cord1)] = "3-CX"
                        stab_to_data[(coords, new_cord2)] = "4-CX"

                    elif string in {"X-STAB-BOUND-A-A", "X-STAB-BOUND-A-T"}:
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                        stab_to_data[(new_cord1, coords)] = "4-CX"
                        stab_to_data[(new_cord2, coords)] = "3-CX"

                    elif string in {"X-STAB-BOUND-B-A", "X-STAB-BOUND-B-T", "X-STAB-SURGERY-B"}:
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                        stab_to_data[(new_cord1, coords)] = "2-CX"
                        stab_to_data[(new_cord2, coords)] = "1-CX"

    # Running the Functions to populate the Schedule
    attach_interior(merging, merging_type)
    attach_boundary(merging, merging_type)
    return stab_to_data


def _populate_surface(
    patch: dict[Coord, Label],
    is_flipped: bool = False,
    y_basis: bool = False,
    y_switch: bool = False,
    y_memory: bool = False,
    distance: int = 0,
    offset: complex = 0 + 0j,
):
    """
    Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

    * X-Stabs: control ist the **Data** qubit -> (data, stab)
    * Z-Stabs: control is the **stab** qubit -> (stab, data)

    * is_fipped -> Switching the role of X and Z stabilizers
    * y_basis -> Initlizing one z edge and one x edge. Each edge is build up out of two sides of the lattice
    * y_switch -> CX schedule after switch (H & SQRT X DEG gates) -> Implementation of XCY Schedule
    """

    stab_to_data: dict[Pair, str] = {}
    stab_to_data_xcy: dict[Pair, str] = {}

    def attach_interior():
        # ------------------------------
        # Internal helper function
        # ------------------------------

        def _neighbours(c: complex, dx: int, dy: int) -> complex:
            return (c.real + dx) + (c.imag + dy) * 1j

        def _quad(c: complex):
            return (
                _neighbours(c, +1, -1),
                _neighbours(c, -1, -1),
                _neighbours(c, +1, +1),
                _neighbours(c, -1, +1),
            )

        def _assign_orders(table: dict[Pair, str], pairs: list[Pair], orders: list[str]) -> None:
            for (a, b), order in zip(pairs, orders, strict=True):
                table[(a, b)] = order

        orders_x_normal = ["1-CX", "2-CX", "3-CX", "4-CX"]
        orders_x_ybasis = ["4-CX", "3-CX", "2-CX", "1-CX"]
        orders_z_normal = ["1-CX", "3-CX", "2-CX", "4-CX"]
        orders_z_ybasis = ["4-CX", "2-CX", "3-CX", "1-CX"]

        if not is_flipped:
            if not y_switch:
                if not y_basis:
                    for coords, qtype in patch.items():
                        q1, q2, q3, q4 = _quad(coords)
                        if qtype == "X-STAB":
                            pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                            _assign_orders(stab_to_data, pairs, orders_x_normal)
                        elif qtype == "Z-STAB":
                            pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                            _assign_orders(stab_to_data, pairs, orders_z_normal)
                else:
                    for coords, qtype in patch.items():
                        q1, q2, q3, q4 = _quad(coords)
                        if qtype == "X-STAB":
                            pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                            _assign_orders(stab_to_data, pairs, orders_x_ybasis)
                        elif qtype == "Z-STAB":
                            pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                            _assign_orders(stab_to_data, pairs, orders_z_ybasis)
            else:
                # Filtering out the stabs needed for the two XCY Gate TICKS
                filtered_stabs_x = []
                filtered_stabs_z = []

                for i in range(distance - 1 if distance > 0 else 0):
                    stabs_z = (2 + i * 2 + offset.real) + (2 + i * 2 + offset.imag) * 1j
                    stabs_x = (4 + i * 2 + offset.real) + (2 + i * 2 + offset.imag) * 1j
                    filtered_stabs_z.append(stabs_z)
                    filtered_stabs_x.append(stabs_x)

                # Filtering out Stabilizers without H applied -> Normal CX-Schedule
                stabs_norm_dict = {}
                stabs_h_dict = {}

                for cords, qtype in patch.items():
                    # Diagonal Cut
                    if cords.real - offset.real <= cords.imag - offset.imag:
                        stabs_norm_dict[cords] = qtype
                    else:
                        stabs_h_dict[cords] = qtype

                # Implementing Solo XCY Gate (Other one in the CX Schedule)
                for coords, _qtype in patch.items():
                    if coords in filtered_stabs_x:
                        new_cord = _neighbours(coords, -1, +1)
                        stab_to_data_xcy[(coords, new_cord)] = "1TICK"

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
                            _assign_orders(
                                stab_to_data, pairs, ["2TICK", "4TICK", "3TICK", "5TICK"],
                            )

                        else:
                            _assign_orders(stab_to_data, [(q3, coords)], ["3.5TICK"])
                            _assign_orders(
                                stab_to_data, [(coords, q2), (coords, q4)], ["4TICK", "5TICK"],
                            )

                """
                I HAVE NO CLUE WHY ONLY WEIGHT 3 instead of weight 4
                """

                # Implementing regular CX scheduele on H half -> only weight 3
                for coords, qtype in stabs_h_dict.items():
                    # Implementing orientation of CX with sub schedule of XCY Gates
                    if qtype == "X-STAB":
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j

                        # Exclude the x stabs next to the diagonal:
                        if coords not in filtered_stabs_x:
                            stab_to_data[(coords, new_cord1)] = "2TICK"

                        stab_to_data[(coords, new_cord2)] = "3TICK"
                        stab_to_data[(new_cord3, coords)] = "4TICK"
                        stab_to_data[(new_cord1, coords)] = "5TICK"

                    elif qtype == "Z-STAB":
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                        new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j

                        stab_to_data[(new_cord2, coords)] = "3TICK"
                        stab_to_data[(coords, new_cord3)] = "4TICK"
                        stab_to_data[(new_cord1, coords)] = "2TICK"

                        # HOTFIX: DATASET WILL OVERWRITE 2TICK if they are the same
                        stab_to_data[coords, new_cord1, 2] = "5TICK"

        else:
            # flipped roles: swap the X/Z assignment orders
            def _assign_flipped():
                for coords, qtype in patch.items():
                    # Defining Needed Neighbour Qubits
                    q1, q2, q3, q4 = _quad(coords)

                    # Control is stab -> (stab, data)
                    if qtype == "X-STAB":
                        pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                        _assign_orders(stab_to_data, pairs, orders_z_normal)

                    # Control is data -> (data, stab)
                    elif qtype == "Z-STAB":
                        pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                        _assign_orders(stab_to_data, pairs, orders_x_normal)

            _assign_flipped()

    def attach_boundary():
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        # (boundary case reference removed; inlined logic below)

        if not is_flipped:
            if not y_basis:
                for coords, qtype in patch.items():
                    if qtype == "Z-STAB-BOUND-L":
                        new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                        stab_to_data[coords, new_cord1] = "1-CX"
                        stab_to_data[coords, new_cord2] = "2-CX"

                    elif qtype == "Z-STAB-BOUND-R":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                        stab_to_data[coords, new_cord1] = "3-CX"
                        stab_to_data[coords, new_cord2] = "4-CX"

                    elif qtype == "X-STAB-BOUND-U":
                        new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                        stab_to_data[new_cord1, coords] = "4-CX"
                        stab_to_data[new_cord2, coords] = "3-CX"

                    elif qtype == "X-STAB-BOUND-B":
                        new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                        new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                        stab_to_data[new_cord1, coords] = "2-CX"
                        stab_to_data[new_cord2, coords] = "1-CX"

            elif y_basis:
                if not y_switch and not y_memory:
                    for coords, qtype in patch.items():
                        if qtype == "Z-STAB-BOUND-L":
                            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                            stab_to_data[coords, new_cord1] = "4-CX"
                            stab_to_data[coords, new_cord2] = "3-CX"

                        elif qtype == "X-STAB-BOUND-R":
                            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                            new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                            stab_to_data[new_cord1, coords] = "3-CX"
                            stab_to_data[new_cord2, coords] = "1-CX"

                        elif qtype == "Z-STAB-BOUND-U":
                            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                            stab_to_data[coords, new_cord1] = "1-CX"
                            stab_to_data[coords, new_cord2] = "3-CX"

                        elif qtype == "X-STAB-BOUND-B":
                            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                            new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
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
                            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                            stab_to_data[coords, new_cord1] = "2TICK"
                            stab_to_data[coords, new_cord2] = "3TICK"

                        elif qtype == "X-STAB-BOUND-R":
                            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                            new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j

                            # Check for lower boundary condition and exclude the cx which gets replaced by CYX
                            lower_boundary = [i for i in range(4, distance * 2, 4)][-1]
                            lower_coord = (distance * 2 + offset.real) + (
                                lower_boundary + offset.imag
                            ) * 1j

                            if coords != lower_coord:
                                stab_to_data[coords, new_cord2] = "2TICK"

                            stab_to_data[coords, new_cord1] = "3TICK"

                        elif qtype == "X-STAB-BOUND-R-H":
                            # H Boundary has mixed cx direction
                            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                            new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                            stab_to_data[coords, new_cord1] = "4TICK"
                            stab_to_data[new_cord2, coords] = "2TICK"
                            stab_to_data[coords, new_cord2] = "5TICK"

                        elif qtype == "Z-STAB-BOUND-U-H":
                            # H Boundary has mixed direction
                            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j

                            if new_cord1 != 1 + 1j + offset:
                                stab_to_data[coords, new_cord1] = "2TICK"

                            stab_to_data[new_cord2, coords] = "4TICK"
                            stab_to_data[new_cord1, coords] = "5TICK"

                        elif qtype == "Z-STAB-BOUND-U":
                            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                            stab_to_data[new_cord1, coords] = "2TICK"
                            stab_to_data[new_cord2, coords] = "3TICK"

                        elif qtype == "X-STAB-BOUND-B":
                            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                            new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                            stab_to_data[new_cord1, coords] = "3TICK"
                            stab_to_data[new_cord2, coords] = "2TICK"

                elif y_memory:
                    for coords, qtype in patch.items():
                        """
                        In the memory Round we switch to the newly introduced 
                        boundary operators and deactivate the old ones
                        """

                        if qtype == "Z-STAB-BOUND-L":
                            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                            stab_to_data[coords, new_cord1] = "4-CX"
                            stab_to_data[coords, new_cord2] = "3-CX"

                        elif qtype == "X-STAB-BOUND-R-H":
                            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                            new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                            stab_to_data[coords, new_cord1] = "2-CX"
                            stab_to_data[coords, new_cord2] = "1-CX"

                        elif qtype == "Z-STAB-BOUND-U-H":
                            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                            stab_to_data[new_cord1, coords] = "1-CX"
                            stab_to_data[new_cord2, coords] = "2-CX"

                        elif qtype == "X-STAB-BOUND-B":
                            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                            new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                            stab_to_data[new_cord1, coords] = "3-CX"
                            stab_to_data[new_cord2, coords] = "4-CX"

        else:
            for coords, qtype in patch.items():
                # Switching the stabilizers so the Z-Stab have switched
                # CX and therefore act like a X-Stab and vice versa

                if qtype == "Z-STAB-BOUND-L":
                    new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "2-CX"
                    stab_to_data[new_cord2, coords] = "1-CX"

                elif qtype == "Z-STAB-BOUND-R":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real - 1) + (coords.imag + 1) * 1j
                    stab_to_data[new_cord1, coords] = "4-CX"
                    stab_to_data[new_cord2, coords] = "3-CX"

                elif qtype == "X-STAB-BOUND-U":
                    new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
                    stab_to_data[coords, new_cord1] = "3-CX"
                    stab_to_data[coords, new_cord2] = "4-CX"

                elif qtype == "X-STAB-BOUND-B":
                    new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
                    new_cord2 = (coords.real + 1) + (coords.imag - 1) * 1j
                    stab_to_data[coords, new_cord1] = "1-CX"
                    stab_to_data[coords, new_cord2] = "2-CX"

    # Running the Functions to populate the Schedule
    attach_interior()
    attach_boundary()
    if stab_to_data_xcy != {}:
        return stab_to_data, stab_to_data_xcy

    return stab_to_data


def populate_stab_to_data(patch: dict[Coord, Label], *args, **kwargs):
    """
    Unified entry point
    -> Detects which style to run and forwards args.
    """

    # Detect xzzx-style by label names
    labels = set(patch.values())
    if any(
        label in labels
        for label in ("STAB-Ver", "STAB-Hor", "STAB-BOUND-L-Hor", "STAB-BOUND-R-Hor")
    ):
        return _populate_xzzx(patch)

    # If explicit merging argument present, treat as lattice_surgery
    if len(args) >= 1 or "merging" in kwargs or "merging_type" in kwargs:
        # Accept either positional (merging, merging_type) or keyword args
        merging = False
        merging_type = None

        if len(args) >= 1:
            merging = args[0]

        if len(args) >= 2:
            merging_type = args[1]

        merging = kwargs.get("merging", merging)
        merging_type = kwargs.get("merging_type", kwargs.get("merging_type", merging_type))

        return _populate_lattice_surgery(patch, merging=bool(merging), merging_type=merging_type)

    # Otherwise treat as surface/rotated code
    # Map surface kwargs with defaults matching original signature
    is_flipped = kwargs.get("is_flipped", kwargs.get("is_flipped", False))
    y_basis = kwargs.get("y_basis", False)
    y_switch = kwargs.get("y_switch", False)
    y_memory = kwargs.get("y_memory", False)
    distance = kwargs.get("distance", 0)
    offset = kwargs.get("offset", 0 + 0j)

    return _populate_surface(
        patch,
        is_flipped=is_flipped,
        y_basis=y_basis,
        y_switch=y_switch,
        y_memory=y_memory,
        distance=distance,
        offset=offset,
    )
