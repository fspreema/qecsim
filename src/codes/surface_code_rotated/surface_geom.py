from src.core.base_geometry import BaseGeometry

"""
Geometry Class which build all Coordinates and converts them to Indices for the Circuits

-> Can later be added together by the __add__ function if no coordinates overlap
-> If overlap, offset can be adjusted
"""


class SurfaceGeometry(BaseGeometry):
    def __init__(
        self,
        distance: int,
        state_init: str,
        logical_observable: str,
        y_basis: bool = False,
        offset: complex = 0 + 0j,
        q2i: dict[complex, int] = None,
    ):
        """
        Initlizes Geometry Class
        Parameters:
            distance : int
            offset : complex, default 0+0j
                Offset of the block in the overall lattice Layout
            q2i : dict[complex, int], optional
                Case where the indices are already defined (e.g. in Lattice Surgery)
                -> Need to reuse the same indices as in the global lattice
        """
        if state_init not in {"1", "0", "+", "-", "+i", "-i"}:
            raise ValueError("state_init must be one of '1', '0', '+', '-', '+i', '-i'")

        if logical_observable not in {"X", "Y", "Z"}:
            raise ValueError("logical_observable must be one of 'X', 'Y', 'Z'")

        if offset.real != 0 and offset.imag != 0:
            raise ValueError("Offset needs to be either 0 or either real or imaginary")

        # Setting up parameters
        self.distance = distance
        self.state_init = state_init
        self.offset = offset
        self.starting_stabilizer_x = False if y_basis else True
        self.y_basis = y_basis
        self.obs = logical_observable

        # Get Coordinates and Indices
        self.coords = self.get_coords(self.y_basis)

        # If q2i is provided use it, otherwise calculate it
        if q2i is not None:
            self.q2i = q2i
        else:
            self.q2i = self._get_q2i()

        self.i2q = self._get_i2q()

        # Setting up y-Index
        self.y_coord = (1 + self.offset.real) + (1 + self.offset.imag) * 1j
        self.y_index = self.q2i.get(self.y_coord, None)
        self.data_idx = self._get_specific_indices("DATA")

        #########################
        # Stabs for x and z basis
        #########################
        self.stab_x_idx = (
            self._get_specific_indices("X-STAB")
            + self._get_specific_indices("X-STAB-BOUND-U")
            + self._get_specific_indices("X-STAB-BOUND-B")
            + self._get_specific_indices("X-STAB-BOUND-R")
        )
        self.stab_z_idx = (
            self._get_specific_indices("Z-STAB")
            + self._get_specific_indices("Z-STAB-BOUND-L")
            + self._get_specific_indices("Z-STAB-BOUND-R")
            + self._get_specific_indices("Z-STAB-BOUND-U")
        )
        self.stab_idx = self.stab_x_idx + self.stab_z_idx

        #####################
        # Indices for Y Basis
        #####################

        # Stabs for H switch in y-basis
        self.stab_switch_apply_h = (
            self._get_specific_indices("X-STAB")
            + self._get_specific_indices("Z-STAB-BOUND-U")
            + self._get_specific_indices("X-STAB-BOUND-B")
            + self._get_specific_indices("Z-STAB-BOUND-U-H")
        )

        # Stabs for memory rounds in y-basis
        self.stab_x_memory = (
            self._get_specific_indices("X-STAB")
            + self._get_specific_indices("Z-STAB-BOUND-U-H")
            + self._get_specific_indices("X-STAB-BOUND-B")
        )
        self.stab_z_memory = (
            self._get_specific_indices("Z-STAB")
            + self._get_specific_indices("Z-STAB-BOUND-L")
            + self._get_specific_indices("X-STAB-BOUND-R-H")
        )

        # Stabs for the switch to y-basis
        self.upper_h = self._get_specific_indices("Z-STAB-BOUND-U-H")
        self.right_h = self._get_specific_indices("X-STAB-BOUND-R-H")

        # Inidces for reset in y-basis
        self.data_rx_idx, self.data_rz_idx = self._y_basis_initial_reset_qubits()


    def get_coords(self, y_basis: bool = False) -> dict[complex, str]:
        """
        Returns all qubit coordinates with their labels

        Returns:
            dict[complex, str]
                Dictionary with coordinates as keys and labels as values
        """

        qubit_coords = self._get_central_labels()
        bound_coords = self._get_boundary_labels(y_basis=y_basis)
        full_coords = qubit_coords | bound_coords

        return full_coords

    def get_logical_observables(
        self,
        logical_observable: str,
        fixed_coord: int = 1,
    ) -> list[int]:
        """
        Returns the indices of the logical observables

        Parameters:
            logical_observable : str
                'X', 'Y' or 'Z' logical observable
            fixed_coord : int
                Non varying coordinate of the logical string (default 1)
                Can be used to shift the logical string on the lattice
        Returns:
            list[int] with respect to currently applied offset
        """
        if logical_observable == "X":
            # Vertical string at x=1 (odd grid), along imag axis
            return [
                self.q2i[fixed_coord + self.offset.real + (imag + self.offset.imag) * 1j]
                for imag in range(1, 2 * self.distance, 2)
            ]

        elif logical_observable == "Z":
            # Horizontal string at y=1, along real axis
            return [
                self.q2i[real + self.offset.real + (fixed_coord + self.offset.imag) * 1j]
                for real in range(1, 2 * self.distance, 2)
            ]
        elif logical_observable == "Y":
            # Combination of X and Z logical strings + the y qubit in the corner
            z_string = [
                self.q2i[real + self.offset.real + (self.offset.imag + fixed_coord) * 1j]
                for real in range(1, 2 * self.distance, 2)
                if real != fixed_coord
            ]
            x_string = [
                self.q2i[fixed_coord + self.offset.real + (imag + self.offset.imag) * 1j]
                for imag in range(1, 2 * self.distance, 2)
                if imag != fixed_coord
            ]
            y_string = [
                self.q2i[fixed_coord + self.offset.real + (self.offset.imag + fixed_coord) * 1j],
            ]

            return (x_string, y_string, z_string)

        else:
            raise ValueError("logical_observable must be 'X', 'Y' or 'Z'")

    def get_neighbors(self, coords: complex, qtype: str) -> list[int]:
        """
        Returns the list of neighboring qubit coords for a given ancilla qubit.
        """

        offsets = {
            "Z-STAB": [-1 - 1j, +1 - 1j, -1 + 1j, +1 + 1j],
            "Z-STAB-BOUND-L": [+1 - 1j, +1 + 1j],
            "Z-STAB-BOUND-R": [-1 - 1j, -1 + 1j],
            "X-STAB": [-1 - 1j, +1 - 1j, -1 + 1j, +1 + 1j],
            "X-STAB-BOUND-U": [-1 + 1j, +1 + 1j],
            "X-STAB-BOUND-B": [-1 - 1j, +1 - 1j],
        }

        neighbor_coords = [coords + offset for offset in offsets[qtype]]

        return [self.q2i[coord] for coord in neighbor_coords if coord in self.q2i]

    def _y_basis_initial_reset_qubits(self) -> tuple[list[int], list[int]]:
        """
        Resets the data qubits in the needed basis for y-basis initlization

        Look at Crumble circuit for a better understanding
        -> Half Half initlization of x and z basis (Cut diagonal)
        """

        data_rx = []
        data_rz = []

        xs = [self.i2q[i].real - self.offset.real for i in self.data_idx]
        ys = [self.i2q[i].imag - self.offset.imag for i in self.data_idx]

        # Calc threshold for diagonal cut
        s0 = (min(xs) + max(xs)) / 2 + (min(ys) + max(ys)) / 2

        skip_coord = 1 + 1j + self.offset

        for data_index in self.data_idx:
            c = self.i2q[data_index]
            if c == skip_coord:
                continue

            # Diagonal Cut
            if (c.real - self.offset.real + c.imag - self.offset.imag) >= s0:
                data_rz.append(self.q2i[c])
            else:
                data_rx.append(self.q2i[c])

        return data_rx, data_rz

    def _y_basis_get_switch_h_qubits(self) -> list[int]:
        """
        Returns the qubits which need to be applied H during the y-basis switch/ rev switch
        """

        h_qubits = []

        # Diagonal Cut
        for cords, _qtype in self.coords.items():
            if cords != self.y_coord:
                if cords.real - self.offset.real > cords.imag - self.offset.imag:
                    h_qubits.append(self.q2i[cords])

        return h_qubits

    def _y_basis_get_switch_xdag_qubits(self) -> list[int]:
        """
        Returns the qubits which need to be applied X_DAG during the y-basis switch/ rev switch
        """

        xdag_qubits = []

        # Filtering out the X_DAG -> Not on Data
        for cords, qtype in self.coords.items():
            if cords.real - self.offset.real == cords.imag - self.offset.imag:
                if qtype != "DATA":
                    xdag_qubits.append(self.q2i[cords])

        return xdag_qubits

    def _get_central_labels(self) -> dict[complex, str]:
        """
        Returns the qubit coordinates depending on the type of the block

        Parameters:
            starting_stabilizer_x : bool
                Whether the first stabilizer is an X stabilizer (True) or Z stabilizer (False)
        Returns:
            dict[complex, str]
                Dictionary with coordinates as keys and labels as values
        """

        if self.distance <= 2 or self.distance % 2 == 0:
            raise ValueError("distance must be odd and ≥3")

        ox, oy = int(self.offset.real), int(self.offset.imag)

        qubit_coords: dict[tuple[complex, complex], str] = {}
        start_with_x = self.starting_stabilizer_x

        for real in range(self.distance * 2):
            stab_counter = 0

            for imag in range(self.distance * 2):
                coord = complex(real + ox, imag + oy)

                # ---------------------- DATA qubits ---------------------------
                if real % 2 != 0 and imag % 2 != 0:
                    qubit_coords[coord] = "DATA"

                # ----------------- Interior Stabilisers ----------------------
                elif real % 2 == 0 and imag % 2 == 0 and real != 0 and imag != 0:
                    if start_with_x:
                        use_x = stab_counter % 2 == 0
                    else:
                        use_x = stab_counter % 2 == 1

                    qubit_coords[coord] = "X-STAB" if use_x else "Z-STAB"

                    stab_counter += 1

            # flip phase after each even row (except first)
            if real % 2 == 0 and real != 0:
                start_with_x = not start_with_x

        return qubit_coords

    def _get_boundary_labels(self, y_basis: bool = False) -> dict[complex, str]:
        """
        Adds the neseccary Boundary and Surgery Stabilizers needed for the code
        """

        qubit_coords: dict[complex, str] = {}

        if self.offset.real != 0 and self.offset.imag != 0:
            raise ValueError("Offset needs to be either 0 or equal in real and imaginary part")

        max_coord = 2 * self.distance

        if not y_basis:
            # Boundary stabilizers
            for y in range(2, max_coord, 4):
                coord_ancilla = complex(0, y)
                qubit_coords[coord_ancilla] = "Z-STAB-BOUND-L"
                coord_ancilla = complex(y, max_coord)
                qubit_coords[coord_ancilla] = "X-STAB-BOUND-B"

            for y in range(4, max_coord, 4):
                coord_ancilla = complex(max_coord, y)
                qubit_coords[coord_ancilla] = "Z-STAB-BOUND-R"
                coord_ancilla = complex(y, 0)
                qubit_coords[coord_ancilla] = "X-STAB-BOUND-U"

        else:
            # Z-boundary stabilizers
            if self.offset == 0 + 0j:
                # Z-boundary stabilizers
                for y in range(4, max_coord, 4):
                    coord_ancilla = complex(y, max_coord)
                    qubit_coords[coord_ancilla] = "X-STAB-BOUND-B"
                    coord_ancilla = complex(max_coord, y)
                    qubit_coords[coord_ancilla] = "X-STAB-BOUND-R"
                    coord_ancilla = complex(y, 0)
                    qubit_coords[coord_ancilla] = "Z-STAB-BOUND-U"
                    coord_ancilla = complex(0, y)
                    qubit_coords[coord_ancilla] = "Z-STAB-BOUND-L"

                # Additional after H
                for y in range(2, max_coord, 4):
                    coord_ancilla = complex(max_coord, y)
                    qubit_coords[coord_ancilla] = "X-STAB-BOUND-R-H"
                    coord_ancilla = complex(y, 0)
                    qubit_coords[coord_ancilla] = "Z-STAB-BOUND-U-H"

            # Check real offset
            elif self.offset.real != 0:
                # Z-boundary stabilizers
                for y in range(4, max_coord, 4):
                    coord_ancilla = complex(y + self.offset.real, max_coord)
                    qubit_coords[coord_ancilla] = "X-STAB-BOUND-B"
                    coord_ancilla = complex(max_coord + self.offset.real, y)
                    qubit_coords[coord_ancilla] = "X-STAB-BOUND-R"
                    coord_ancilla = complex(y + self.offset.real, 0)
                    qubit_coords[coord_ancilla] = "Z-STAB-BOUND-U"
                    coord_ancilla = complex(self.offset.real, y)
                    qubit_coords[coord_ancilla] = "Z-STAB-BOUND-L"

                # Additional after H
                for y in range(2, max_coord, 4):
                    coord_ancilla = complex(max_coord + self.offset.real, y)
                    qubit_coords[coord_ancilla] = "X-STAB-BOUND-R-H"
                    coord_ancilla = complex(y + self.offset.real, 0)
                    qubit_coords[coord_ancilla] = "Z-STAB-BOUND-U-H"

            # Check imaginary offset
            elif self.offset.imag != 0:
                # Z-boundary stabilizers
                for y in range(4, max_coord, 4):
                    coord_ancilla = complex(y, max_coord + self.offset.imag)
                    qubit_coords[coord_ancilla] = "X-STAB-BOUND-B"
                    coord_ancilla = complex(max_coord, y + self.offset.imag)
                    qubit_coords[coord_ancilla] = "X-STAB-BOUND-R"
                    coord_ancilla = complex(y, 0 + self.offset.imag)
                    qubit_coords[coord_ancilla] = "Z-STAB-BOUND-U"
                    coord_ancilla = complex(0, y + self.offset.imag)
                    qubit_coords[coord_ancilla] = "Z-STAB-BOUND-L"

                # Additional after H
                for y in range(2, max_coord, 4):
                    coord_ancilla = complex(max_coord, y + self.offset.imag)
                    qubit_coords[coord_ancilla] = "X-STAB-BOUND-R-H"
                    coord_ancilla = complex(y, 0 + self.offset.imag)
                    qubit_coords[coord_ancilla] = "Z-STAB-BOUND-U-H"

        return qubit_coords
