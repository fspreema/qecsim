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

        self.distance = distance
        self.state_init = state_init
        self.offset = offset
        self.starting_stabilizer_x = starting_stabilizer_x
        self.coords = self._get_qubit_coords() | self._add_boundary_labels()

    def _get_qubit_coords(self) -> dict[complex, str]:
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
                        qubit_coords[coord] = "DATA_X" if use_x else "DATA_Z"
                    else:  # XZZX-HOR
                        qubit_coords[coord] = "DATA_Z" if use_x else "DATA_X"

                    data_counter += 1

                # ----------------- Interior Stabilisers ----------------------
                elif real % 2 == 0 and imag % 2 == 0 and real != 0 and imag != 0:
                    if start_with_x:
                        use_x = stab_counter % 2 == 0
                    else:
                        use_x = stab_counter % 2 == 1

                    qubit_coords[coord] = "STAB-Ver" if use_x else "STAB-Hor"

                    stab_counter += 1

            # flip phase after each even row (except first)
            if real % 2 == 0 and real != 0:
                start_with_x = not start_with_x

        return qubit_coords

    def _add_boundary_labels(self) -> dict[complex, str]:
        """
        Adds the neseccary Boundary and Surgery Stabilizers needed
        """

        max_coord = 2 * self.distance

        qubit_coords_bound: dict[complex, str] = {}

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
