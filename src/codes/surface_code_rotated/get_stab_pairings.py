from src.core.base_get_stab_pairings import BasePairings

Coord = complex
Label = str
Pair = tuple[Coord, Coord]

__all__ = ["SurfacePairings"]


class StandardPairings(BasePairings):
    def __init__(
        self,
        patch: dict[complex, str],
        distance: int,
        is_flipped: bool = False,
        y_basis: bool = False,
        y_switch: bool = False,
        y_memory: bool = False,
        offset: complex = 0 + 0j,
    ):
        # Initialize Parameters
        self.pairings = patch
        self.is_flipped = is_flipped
        self.y_basis = y_basis
        self.y_switch = y_switch
        self.y_memory = y_memory
        self.distance = distance
        self.offset = offset

    def generate_standard_pairings(self):
        # Initialize Schedule Dictionary
        self.stab_to_data: dict[Pair, str] = {}
        self.stab_to_data_xcy: dict[Pair, str] = {}

        # Running the Functions to populate the Schedule
        self._attach_interior_standard()
        self._attach_boundary_standard()

        if self.y_switch:
            return self.stab_to_data, self.stab_to_data_xcy
        else:
            return self.stab_to_data

    def _attach_interior_standard(self):
        for coords, qtype in self.pairings.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if qtype == "X-STAB":
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX", "3-CX", "4-CX"],
                )

            elif qtype == "Z-STAB":
                pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "3-CX", "2-CX", "4-CX"],
                )

    def _attach_boundary_standard(self):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        for coords, qtype in self.pairings.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if qtype == "Z-STAB-BOUND-L":
                pairs = [(coords, q1), (coords, q3)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX"],
                )

            elif qtype == "Z-STAB-BOUND-R":
                pairs = [(coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["3-CX", "4-CX"],
                )

            elif qtype == "X-STAB-BOUND-U":
                pairs = [(q4, coords), (q3, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "3-CX"],
                )

            elif qtype == "X-STAB-BOUND-B":
                pairs = [(q2, coords), (q1, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["2-CX", "1-CX"],
                )


class YBasisPairings(BasePairings):
    def __init__(
        self,
        patch: dict[complex, str],
        distance: int,
        is_flipped: bool = False,
        y_basis: bool = False,
        y_switch: bool = False,
        y_memory: bool = False,
        offset: complex = 0 + 0j,
    ):
        # Initialize Parameters
        self.pairings = patch
        self.is_flipped = is_flipped
        self.y_basis = y_basis
        self.y_switch = y_switch
        self.y_memory = y_memory
        self.distance = distance
        self.offset = offset

    def generate_ybasis_pairings(self):
        # Initialize Schedule Dictionary
        self.stab_to_data: dict[Pair, str] = {}
        self.stab_to_data_xcy: dict[Pair, str] = {}

        # Running the Functions to populate the Schedule
        self._attach_interior_ybasis()
        self._attach_boundary_ybasis()

        if self.y_switch:
            return self.stab_to_data, self.stab_to_data_xcy
        else:
            return self.stab_to_data

    def _attach_interior_ybasis(self):
        for coords, qtype in self.pairings.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if qtype == "X-STAB":
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "3-CX", "2-CX", "1-CX"],
                )

            elif qtype == "Z-STAB":
                pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "2-CX", "3-CX", "1-CX"],
                )

    def _attach_boundary_ybasis(self):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        if not self.y_switch and not self.y_memory:
            for coords, qtype in self.pairings.items():
                # Get neigbouring data coords
                q1 = self._neighbours(coords, +1, -1)
                q2 = self._neighbours(coords, -1, -1)
                q3 = self._neighbours(coords, +1, +1)
                q4 = self._neighbours(coords, -1, +1)

                if qtype == "Z-STAB-BOUND-L":
                    pairs = [(coords, q1), (coords, q3)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["4-CX", "3-CX"],
                    )

                elif qtype == "X-STAB-BOUND-R":
                    pairs = [(q2, coords), (q4, coords)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["3-CX", "1-CX"],
                    )

                elif qtype == "Z-STAB-BOUND-U":
                    pairs = [(coords, q4), (coords, q3)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["1-CX", "3-CX"],
                    )

                elif qtype == "X-STAB-BOUND-B":
                    pairs = [(q2, coords), (q1, coords)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["3-CX", "4-CX"],
                    )

        elif self.y_switch and not self.y_memory:
            for coords, qtype in self.pairings.items():
                """
                As the H gates was applied we need to flip the corressponding
                stabilizer schedule

                * ATTENTION *
                -> ONLY FLIP THESE WHICH WHERE FLIPPED I.E. on only one diagonal half
                -> Additional Boundary Stabs do CX both ways i.e. detecting z and x 
                    errors!
                """
                # Get neigbouring data coords
                q1 = self._neighbours(coords, +1, -1)
                q2 = self._neighbours(coords, -1, -1)
                q3 = self._neighbours(coords, +1, +1)
                q4 = self._neighbours(coords, -1, +1)

                if qtype == "Z-STAB-BOUND-L":
                    pairs = [(coords, q1), (coords, q3)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["2TICK", "3TICK"],
                    )

                elif qtype == "X-STAB-BOUND-R":
                    # Check for lower boundary condition and exclude
                    # the cx which gets replaced by CYX
                    lower_boundary = [i for i in range(4, self.distance * 2, 4)][-1]
                    lower_coord = (self.distance * 2 + self.offset.real) + (
                        lower_boundary + self.offset.imag
                    ) * 1j

                    if coords != lower_coord:
                        pairs = [(coords, q2), (coords, q4)]
                        self._assign_orders(
                            self.stab_to_data,
                            pairs,
                            ["3TICK", "2TICK"],
                        )
                    else:
                        pairs = [(coords, q2)]
                        self._assign_orders(
                            self.stab_to_data,
                            pairs,
                            ["3TICK"],
                        )

                elif qtype == "X-STAB-BOUND-R-H":
                    # H Boundary has mixed cx direction
                    pairs = [(coords, q2), (q4, coords), (coords, q4)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["4TICK", "2TICK", "5TICK"],
                    )

                elif qtype == "Z-STAB-BOUND-U-H":
                    # H Boundary has mixed direction
                    if q4 != 1 + 1j + self.offset:
                        pairs = [(coords, q4), (q3, coords), (q4, coords)]
                        self._assign_orders(
                            self.stab_to_data,
                            pairs,
                            ["2TICK", "4TICK", "5TICK"],
                        )
                    else:
                        pairs = [(q3, coords), (q4, coords)]
                        self._assign_orders(
                            self.stab_to_data,
                            pairs,
                            ["4TICK", "5TICK"],
                        )

                elif qtype == "Z-STAB-BOUND-U":
                    pairs = [(q4, coords), (q3, coords)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["2TICK", "3TICK"],
                    )

                elif qtype == "X-STAB-BOUND-B":
                    pairs = [(q2, coords), (q1, coords)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["3TICK", "2TICK"],
                    )

        elif self.y_memory:
            for coords, qtype in self.pairings.items():
                """
                In the memory Round we switch to the newly introduced
                boundary operators and deactivate the old ones
                """
                # Get neigbouring data coords
                q1 = self._neighbours(coords, +1, -1)
                q2 = self._neighbours(coords, -1, -1)
                q3 = self._neighbours(coords, +1, +1)
                q4 = self._neighbours(coords, -1, +1)

                if qtype == "Z-STAB-BOUND-L":
                    pairs = [(coords, q1), (coords, q3)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["4-CX", "3-CX"],
                    )

                elif qtype == "X-STAB-BOUND-R-H":
                    pairs = [(coords, q2), (coords, q4)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["2-CX", "1-CX"],
                    )

                elif qtype == "Z-STAB-BOUND-U-H":
                    pairs = [(q4, coords), (q3, coords)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["1-CX", "2-CX"],
                    )

                elif qtype == "X-STAB-BOUND-B":
                    pairs = [(q2, coords), (q1, coords)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["3-CX", "4-CX"],
                    )


class YSwitchPairings(BasePairings):
    def __init__(
        self,
        patch: dict[complex, str],
        distance: int,
        is_flipped: bool = False,
        y_basis: bool = False,
        y_switch: bool = False,
        y_memory: bool = False,
        offset: complex = 0 + 0j,
    ):
        # Initialize Parameters
        self.pairings = patch
        self.is_flipped = is_flipped
        self.y_basis = y_basis
        self.y_switch = y_switch
        self.y_memory = y_memory
        self.distance = distance
        self.offset = offset

    def get_yswitch_schedule(self):
        # Initialize Schedule Dictionary
        self.stab_to_data: dict[Pair, str] = {}
        self.stab_to_data_xcy: dict[Pair, str] = {}

        # Running the Functions to populate the Schedule
        self._attach_interior_yswitch()
        self._attach_boundary_yswitch()

        if self.y_switch:
            return self.stab_to_data, self.stab_to_data_xcy
        else:
            return self.stab_to_data

    def _attach_interior_yswitch(self):
        # Filtering out the stabs needed for the two XCY Gate TICKS
        filtered_stabs_x = []
        filtered_stabs_z = []

        for i in range(self.distance - 1 if self.distance > 0 else 0):
            stabs_z = (2 + i * 2 + self.offset.real) + (2 + i * 2 + self.offset.imag) * 1j
            stabs_x = (4 + i * 2 + self.offset.real) + (2 + i * 2 + self.offset.imag) * 1j
            filtered_stabs_z.append(stabs_z)
            filtered_stabs_x.append(stabs_x)

        # Filtering out Stabilizers without H applied -> Normal CX-Schedule
        stabs_norm_dict = {}
        stabs_h_dict = {}

        for cords, qtype in self.pairings.items():
            # Diagonal Cut
            if cords.real - self.offset.real <= cords.imag - self.offset.imag:
                stabs_norm_dict[cords] = qtype
            else:
                stabs_h_dict[cords] = qtype

        # Implementing Solo XCY Gate (Other one in the CX Schedule)
        for coords, _qtype in self.pairings.items():
            if coords in filtered_stabs_x:
                new_cord = self._neighbours(coords, -1, +1)
                self.stab_to_data_xcy[(coords, new_cord)] = "1TICK"

        for coords, qtype in stabs_norm_dict.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if qtype == "X-STAB":
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["2TICK", "3TICK", "4TICK", "5TICK"],
                )

            elif qtype == "Z-STAB":
                # Checking whether normal CX or the XCY gate
                if coords not in filtered_stabs_z:
                    pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["2TICK", "4TICK", "3TICK", "5TICK"],
                    )

                else:
                    self._assign_orders(self.stab_to_data, [(q3, coords)], ["3.5TICK"])
                    self._assign_orders(
                        self.stab_to_data,
                        [(coords, q2), (coords, q4)],
                        ["4TICK", "5TICK"],
                    )

        """
        I HAVE NO CLUE WHY ONLY WEIGHT 3 instead of weight 4
        """

        # Implementing regular CX scheduele on H half -> only weight 3
        for coords, qtype in stabs_h_dict.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            # Implementing orientation of CX with sub schedule of XCY Gates
            if qtype == "X-STAB":
                # Exclude the x stabs next to the diagonal:
                if coords not in filtered_stabs_x:
                    self.stab_to_data[(coords, q4)] = "2TICK"

                pairs = [(coords, q2), (q3, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["3TICK", "4TICK", "5TICK"],
                )

            elif qtype == "Z-STAB":
                pairs = [(q3, coords), (coords, q2), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["3TICK", "4TICK", "2TICK"],
                )

                # HOTFIX: DATASET WILL OVERWRITE 2TICK if they are the same
                self.stab_to_data[coords, q4, 2] = "5TICK"

    def _attach_boundary_yswitch(self):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        for coords, qtype in self.pairings.items():
            """
            As the H gates was applied we need to flip the corressponding
            stabilizer schedule

            * ATTENTION *
            -> ONLY FLIP THESE WHICH WHERE FLIPPED I.E. on only one diagonal half
            -> Additional Boundary Stabs do CX both ways i.e. detecting z and x 
                errors!
            """
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if qtype == "Z-STAB-BOUND-L":
                pairs = [(coords, q1), (coords, q3)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["2TICK", "3TICK"],
                )

            elif qtype == "X-STAB-BOUND-R":
                # Check for lower boundary condition and exclude
                # the cx which gets replaced by CYX
                lower_boundary = [i for i in range(4, self.distance * 2, 4)][-1]
                lower_coord = (self.distance * 2 + self.offset.real) + (
                    lower_boundary + self.offset.imag
                ) * 1j

                if coords != lower_coord:
                    pairs = [(coords, q2), (coords, q4)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["3TICK", "2TICK"],
                    )
                else:
                    pairs = [(coords, q2)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["3TICK"],
                    )

            elif qtype == "X-STAB-BOUND-R-H":
                # H Boundary has mixed cx direction
                pairs = [(coords, q2), (q4, coords), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4TICK", "2TICK", "5TICK"],
                )

            elif qtype == "Z-STAB-BOUND-U-H":
                # H Boundary has mixed direction
                if q4 != 1 + 1j + self.offset:
                    pairs = [(coords, q4), (q3, coords), (q4, coords)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["2TICK", "4TICK", "5TICK"],
                    )
                else:
                    pairs = [(q3, coords), (q4, coords)]
                    self._assign_orders(
                        self.stab_to_data,
                        pairs,
                        ["4TICK", "5TICK"],
                    )

            elif qtype == "Z-STAB-BOUND-U":
                pairs = [(q4, coords), (q3, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["2TICK", "3TICK"],
                )

            elif qtype == "X-STAB-BOUND-B":
                pairs = [(q2, coords), (q1, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["3TICK", "2TICK"],
                )


class YMemoryPairings(BasePairings):
    def __init__(
        self,
        patch: dict[complex, str],
        distance: int,
        is_flipped: bool = False,
        y_basis: bool = False,
        y_switch: bool = False,
        y_memory: bool = False,
        offset: complex = 0 + 0j,
    ):
        # Initialize Parameters
        self.pairings = patch
        self.is_flipped = is_flipped
        self.y_basis = y_basis
        self.y_switch = y_switch
        self.y_memory = y_memory
        self.distance = distance
        self.offset = offset

    def get_ymemory_schedule(self):
        # Initialize Schedule Dictionary
        self.stab_to_data: dict[Pair, str] = {}
        self.stab_to_data_xcy: dict[Pair, str] = {}

        # Running the Functions to populate the Schedule
        self._attach_interior_ymemory()
        self._attach_boundary_ymemory()

        if self.y_switch:
            return self.stab_to_data, self.stab_to_data_xcy
        else:
            return self.stab_to_data

    def _attach_interior_ymemory(self):
        for coords, qtype in self.pairings.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if qtype == "X-STAB":
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "3-CX", "2-CX", "1-CX"],
                )

            elif qtype == "Z-STAB":
                pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "2-CX", "3-CX", "1-CX"],
                )

    def _attach_boundary_ymemory(self):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        for coords, qtype in self.pairings.items():
            """
            In the memory Round we switch to the newly introduced
            boundary operators and deactivate the old ones
            """
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if qtype == "Z-STAB-BOUND-L":
                pairs = [(coords, q1), (coords, q3)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "3-CX"],
                )

            elif qtype == "X-STAB-BOUND-R-H":
                pairs = [(coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["2-CX", "1-CX"],
                )

            elif qtype == "Z-STAB-BOUND-U-H":
                pairs = [(q4, coords), (q3, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX"],
                )

            elif qtype == "X-STAB-BOUND-B":
                pairs = [(q2, coords), (q1, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["3-CX", "4-CX"],
                )


class FlippedPairings(BasePairings):
    def __init__(
        self,
        patch: dict[complex, str],
        distance: int,
        is_flipped: bool = False,
        y_basis: bool = False,
        y_switch: bool = False,
        y_memory: bool = False,
        offset: complex = 0 + 0j,
    ):
        # Initialize Parameters
        self.pairings = patch
        self.is_flipped = is_flipped
        self.y_basis = y_basis
        self.y_switch = y_switch
        self.y_memory = y_memory
        self.distance = distance
        self.offset = offset

    def get_flipped_schedule(self):
        # Initialize Schedule Dictionary
        self.stab_to_data: dict[Pair, str] = {}
        self.stab_to_data_xcy: dict[Pair, str] = {}

        # Running the Functions to populate the Schedule
        self._attach_interior_flipped()
        self._attach_boundary_flipped()

        if self.y_switch:
            return self.stab_to_data, self.stab_to_data_xcy
        else:
            return self.stab_to_data

    def _attach_interior_flipped(self):
        # flipped roles: swap the X/Z assignment orders

        for coords, qtype in self.pairings.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            # Control is stab -> (stab, data)
            if qtype == "X-STAB":
                pairs = [(coords, q1), (coords, q2), (coords, q3), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "3-CX", "2-CX", "4-CX"],
                )

            # Control is data -> (data, stab)
            elif qtype == "Z-STAB":
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX", "3-CX", "4-CX"],
                )

    def _attach_boundary_flipped(self):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        for coords, qtype in self.pairings.items():
            # Switching the stabilizers so the Z-Stab have switched
            # CX and therefore act like a X-Stab and vice versa
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if qtype == "Z-STAB-BOUND-L":
                pairs = [(q1, coords), (q3, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["2-CX", "1-CX"],
                )

            elif qtype == "Z-STAB-BOUND-R":
                pairs = [(q2, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "3-CX"],
                )

            elif qtype == "X-STAB-BOUND-U":
                pairs = [(coords, q4), (coords, q3)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["3-CX", "4-CX"],
                )

            elif qtype == "X-STAB-BOUND-B":
                pairs = [(coords, q2), (coords, q1)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX"],
                )


class SurfacePairings(
    StandardPairings,
    YBasisPairings,
    YSwitchPairings,
    YMemoryPairings,
    FlippedPairings,
):
    def __init__(
        self,
        patch: dict[complex, str],
        distance: int,
        is_flipped: bool = False,
        y_basis: bool = False,
        y_switch: bool = False,
        y_memory: bool = False,
        offset: complex = 0 + 0j,
    ):
        # Initialize Parameters
        self.pairings = patch
        self.is_flipped = is_flipped
        self.y_basis = y_basis
        self.y_switch = y_switch
        self.y_memory = y_memory
        self.distance = distance
        self.offset = offset

    def get_schedule(self):
        """
        Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

        * X-Stabs: control ist the **Data** qubit -> (data, stab)
        * Z-Stabs: control is the **stab** qubit -> (stab, data)

        * is_flipped -> Switching the role of X and Z stabilizers
        * y_basis -> Initlizing one z edge and one x edge.
                    Each edge is build up out of two sides of the lattice
        * y_switch -> CX schedule after switch (H & SQRT X DEG gates) -> Implementation of XCY
                      Schedule
        """

        # Initialize Schedule Dictionary
        self.stab_to_data: dict[Pair, str] = {}
        self.stab_to_data_xcy: dict[Pair, str] = {}

        # Running the Functions to populate the Schedule
        if not self.is_flipped and not self.y_basis and not self.y_switch and not self.y_memory:
            # Generate Pairings and return dict
            self.generate_standard_pairings()
            return self.stab_to_data

        elif not self.is_flipped and self.y_basis and not self.y_switch and not self.y_memory:
            # Generate Pairings and return dict
            self.generate_ybasis_pairings()
            return self.stab_to_data

        elif not self.is_flipped and self.y_basis and self.y_switch and not self.y_memory:
            # Generate Pairings and return dict
            self.get_yswitch_schedule()
            return self.stab_to_data, self.stab_to_data_xcy

        elif not self.is_flipped and self.y_basis and not self.y_switch and self.y_memory:
            # Generate Pairings and return dict
            self.get_ymemory_schedule()
            return self.stab_to_data

        elif self.is_flipped and not self.y_basis and not self.y_switch and not self.y_memory:
            # Generate Pairings and return dict
            self.get_flipped_schedule()
            return self.stab_to_data
