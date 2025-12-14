from src.core.base_geometry import BaseGeometry

"""
Geometry Class which build all Coordinates and converts them to Indices for the Circuits

-> Can later be added together by the __add__ function if no coordinates overlap
-> If overlap, offset can be adjusted
"""


class SurgeryGeometry(BaseGeometry):
    def __init__(
        self,
        distance: int,
    ):
        """
        Initializes Geometry Class for Lattice Surgery.
        Generates 3 patches (Ancilla, Target, Control) and connecting boundaries.

        Layout:
            Ancilla (Top-Left) | Target (Top-Right)
            -------------------|-------------------
            Control (Bot-Left) | (Empty/Surgery)
        """

        # Initialize Parameters
        self.distance = distance
        self.offset_ancilla = 0 + 0j
        self.offset_target = (self.distance * 2) + 0j
        self.offset_control = 0 + (self.distance * 2) * 1j
        self.start_stab_x_ancilla = True
        self.start_stab_x_target = False
        self.start_stab_x_control = False
        self.type = "Surgery"

        """
        Different Patches are needed, because of different Keywords on 
        identical Coordinates (inside dict.):

        -> X-Stab-Boundary-Above-Control & X-Stab-Boundary-Below-Ancilla f.ex. 
        get Keywords for surgery stabilizers
        """

        # Generate Coordinates for each patch
        self.coords_ancilla = self._get_qubit_coords(self.offset_ancilla, self.start_stab_x_ancilla)
        self.coords_target = self._get_qubit_coords(self.offset_target, self.start_stab_x_target)
        self.coords_control = self._get_qubit_coords(self.offset_control, self.start_stab_x_control)
        self.coords_surgery: dict[complex, str] = {}

        # Add Boundary and Surgery Stabilizers
        self._add_boundary_labels()

        # Full Coordinate Dictionary
        self.coords = (
            self.coords_ancilla | self.coords_target | self.coords_control | self.coords_surgery
        )

    def _get_qubit_coords(self, offset: complex, starting_stabilizer_x: bool) -> dict[complex, str]:
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

        ox, oy = int(offset.real), int(offset.imag)

        qubit_coords: dict[tuple[complex, complex], str] = {}
        start_with_x = starting_stabilizer_x

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

    def _add_boundary_labels(self) -> dict[complex, str]:
        """
        Adds the neseccary Boundary and Surgery Stabilizers needed
        """

        max_coord = 2 * self.distance

        # Z-boundary stabilizers
        for y in range(2, max_coord, 4):
            coord_ancilla = complex(0, y)
            coord_target = complex(self.distance * 2, y + 2)
            coord_control = complex(0, y + 2 + (self.distance * 2))
            coord_surgery = complex(y + 2, max_coord)
            self.coords_ancilla[coord_ancilla] = "Z-STAB-BOUND-L-A"
            self.coords_target[coord_target] = "Z-STAB-BOUND-L-T"
            self.coords_control[coord_control] = "Z-STAB-BOUND-L-C"
            self.coords_surgery[coord_surgery] = "Z-STAB-SURGERY-M"

        for y in range(4, max_coord, 4):
            coord_ancilla = complex(max_coord, y)
            coord_target = complex(max_coord + (self.distance * 2), y - 2)
            coord_control = complex(max_coord, y - 2 + (self.distance * 2))
            self.coords_ancilla[coord_ancilla] = "Z-STAB-BOUND-R-A"
            self.coords_target[coord_target] = "Z-STAB-BOUND-R-T"
            self.coords_control[coord_control] = "Z-STAB-BOUND-R-C"

        # X-boundary stabilizers
        for y in range(4, max_coord, 4):
            coord_ancilla = complex(y, 0)
            coord_target = complex(y - 2 + (self.distance * 2), 0)
            coord_control = complex(y - 2, self.distance * 2)
            coord_surgery = complex(max_coord, y - 2)
            self.coords_ancilla[coord_ancilla] = "X-STAB-BOUND-A-A"
            self.coords_target[coord_target] = "X-STAB-BOUND-A-T"
            self.coords_control[coord_control] = "X-STAB-BOUND-A-C"
            self.coords_surgery[coord_surgery] = "X-STAB-SURGERY-M"

        for y in range(2, max_coord, 4):
            coord_ancilla = complex(y, max_coord)
            coord_target = complex(y + 2 + (self.distance * 2), max_coord)
            coord_control = complex(y + 2, max_coord + (self.distance * 2))
            self.coords_ancilla[coord_ancilla] = "X-STAB-BOUND-B-A"
            self.coords_target[coord_target] = "X-STAB-BOUND-B-T"
            self.coords_control[coord_control] = "X-STAB-BOUND-B-C"

        # Adding X Surgery Stabilizer Between Ancilla & Target
        self.coords_surgery[complex(max_coord, max_coord)] = "X-STAB-SURGERY-B"

        # Adding Z Surgery Stabilizer Between Ancilla & Control
        self.coords_surgery[complex(0, max_coord)] = "Z-STAB-SURGERY-L"

        # Adding aditional side qubits which my be used in order to do the y basis init
        """
        The following dict entries are not really used. 
        -> Instead they are fillers in order for the q2i indexing to work correctly
        -> The real labels are still given by the labeling function from the y circuit
        """

        for y in range(4, max_coord, 4):
            coord_target = complex(y + (self.distance * 2), 0)
            self.coords_target[coord_target] = "Z-STAB-BOUND-U-H"
            coord_target = complex(max_coord + (self.distance * 2), y)
            self.coords_target[coord_target] = "X-STAB-BOUND-R-H"

        for y in range(4, max_coord, 4):
            coord_control = complex((self.distance * 2), y + (self.distance * 2))
            self.coords_control[coord_control] = "X-STAB-BOUND-R-H"
