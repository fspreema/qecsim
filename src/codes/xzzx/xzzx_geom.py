from src.core.base_geometry import BaseGeometry

"""
Geometry Class which build all Coordinates and converts them to Indices for the Circuits

-> Can later be added together by the __add__ function if no coordinates overlap
-> If overlap, offset can be adjusted
"""


class XZZXGeometry(BaseGeometry):
    def __init__(
        self,
        distance: int,
        state_init: str,
        offset: complex = 0 + 0j,
        starting_stabilizer_x: bool = True,
    ):
        """
        Initlizes Geometry Class
        Parameters:
            distance : int
            type : str
                Type of the block (e.g. XZZX or Surface -> Later additional types)
            offset : complex, default 0+0j
                Offset of the block in the overall lattice Layout
        """

        if state_init not in ["XZZX-VER", "XZZX-HOR"]:
            raise ValueError("type must be either 'XZZX-VER' or 'XZZX-HOR'")

        # Setting up parameters
        self.distance = distance
        self.state_init = state_init
        self.offset = offset
        self.starting_stabilizer_x = starting_stabilizer_x

        # Get Coordinates and Indices
        self.qubit_coords: dict[complex, str] = {}
        self.coords = self.get_coords(self.qubit_coords)
        self.q2i = self._get_q2i()
        self.i2q = self._get_i2q()
        self.data_z_idx = self._get_specific_indices("DATA-Z")
        self.data_x_idx = self._get_specific_indices("DATA-X")
        self.data_idx = self.data_z_idx + self.data_x_idx
        self.stab_ver_idx = (
            self._get_specific_indices("STAB-VER")
            + self._get_specific_indices("STAB-BOUND-A-Ver")
            + self._get_specific_indices("STAB-BOUND-B-Ver")
        )
        self.stab_hor_idx = (
            self._get_specific_indices("STAB-HOR")
            + self._get_specific_indices("STAB-BOUND-L-Hor")
            + self._get_specific_indices("STAB-BOUND-R-Hor")
        )
        self.stab_idx = self.stab_ver_idx + self.stab_hor_idx

    def get_neighbors(self, coords: complex, qtype: str) -> list[int]:
        """
        Returns the list of neighboring qubit coords for a given ancilla qubit.
        """

        offsets = {
            "STAB-VER": [-1 - 1j, 1 - 1j, -1 + 1j, 1 + 1j],
            "STAB-HOR": [-1 - 1j, 1 - 1j, -1 + 1j, 1 + 1j],
            "STAB-BOUND-A-Ver": [1 + 1j, -1 + 1j],
            "STAB-BOUND-B-Ver": [1 - 1j, -1 - 1j],
            "STAB-BOUND-L-Hor": [1 - 1j, 1 + 1j],
            "STAB-BOUND-R-Hor": [-1 - 1j, -1 + 1j],
        }

        neighbor_coords = [coords + offset for offset in offsets[qtype]]

        return [self.q2i[coord] for coord in neighbor_coords if coord in self.q2i]

    def get_coords(self, qubit_coords: dict[complex, str]) -> dict[complex, str]:
        """
        Returns all qubit coordinates with their labels

        Returns:
            dict[complex, str]
                Dictionary with coordinates as keys and labels as values
        """

        qubit_coords = self._get_central_labels(qubit_coords)
        bound_coords = self._get_boundary_labels(qubit_coords)
        full_coords = qubit_coords | bound_coords

        return full_coords

    def get_logical_indexes(self, logical_operator: str) -> dict[str, list[int]]:
        """
        Returns the list of Indices corresponding to the logical operators
        """

        if logical_operator not in {"VER", "HOR"}:
            raise ValueError("logical_operator must be either 'VER' or 'HOR'")

        if logical_operator == "VER":
            log_ver: list[int] = []
            for imag in range(1, (self.distance * 2), 2):
                log_ver.append(self.q2i[1 + imag * 1j])

            return log_ver

        if logical_operator == "HOR":
            log_hor: list[int] = []
            for real in range(1, (self.distance * 2), 2):
                log_hor.append(self.q2i[real + 1j])

            return log_hor

    def _get_central_labels(self, qubit_coords: dict[complex, str]) -> dict[complex, str]:
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
        start_with_x = self.starting_stabilizer_x

        for real in range(self.distance * 2):
            stab_counter = 0
            data_counter = 0

            for imag in range(self.distance * 2):
                coord = complex(real + ox, imag + oy)

                # ---------------------- DATA qubits ---------------------------
                if real % 2 != 0 and imag % 2 != 0:
                    if start_with_x:
                        use_x = data_counter % 2 == 0
                    else:
                        use_x = data_counter % 2 == 1

                    if self.state_init == "XZZX-VER":
                        qubit_coords[coord] = "DATA-X" if use_x else "DATA-Z"
                    else:  # XZZX-HOR
                        qubit_coords[coord] = "DATA-Z" if use_x else "DATA-X"

                    data_counter += 1

                # ----------------- Interior Stabilisers ----------------------
                elif real % 2 == 0 and imag % 2 == 0 and real != 0 and imag != 0:
                    if start_with_x:
                        use_x = stab_counter % 2 == 0
                    else:
                        use_x = stab_counter % 2 == 1

                    qubit_coords[coord] = "STAB-VER" if use_x else "STAB-HOR"

                    stab_counter += 1

            # flip phase after each even row (except first)
            if real % 2 == 0 and real != 0:
                start_with_x = not start_with_x

        return qubit_coords

    def _get_boundary_labels(self, qubit_coords_bound: dict[complex, str]) -> dict[complex, str]:
        """
        Adds the neseccary Boundary and Surgery Stabilizers needed
        """

        max_coord = 2 * self.distance

        # Z-boundary stabilizers
        for y in range(2, max_coord, 4):
            coord_bound = complex(0, y)
            qubit_coords_bound[coord_bound] = "STAB-BOUND-L-Hor"

        for y in range(4, max_coord, 4):
            coord_bound = complex(max_coord, y)
            qubit_coords_bound[coord_bound] = "STAB-BOUND-R-Hor"

        # X-boundary stabilizers
        for y in range(4, max_coord, 4):
            coord_bound = complex(y, 0)
            qubit_coords_bound[coord_bound] = "STAB-BOUND-A-Ver"

        for y in range(2, max_coord, 4):
            coord_bound = complex(y, max_coord)
            qubit_coords_bound[coord_bound] = "STAB-BOUND-B-Ver"

        return qubit_coords_bound
