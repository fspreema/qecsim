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
        y_basis: bool = False,
        offset: complex = 0 + 0j,
        starting_stabilizer_x: bool = True,
    ):
        """
        Initlizes Geometry Class
        Parameters:
            distance : int
            offset : complex, default 0+0j
                Offset of the block in the overall lattice Layout
        """

        self.distance = distance
        self.offset = offset
        self.starting_stabilizer_x = starting_stabilizer_x
        self.coords = self.get_coords(y_basis=y_basis)

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
