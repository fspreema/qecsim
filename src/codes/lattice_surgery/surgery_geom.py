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
        control_state_init: str,
        target_state_init: str,
    ):
        """
        Initializes Geometry Class for Lattice Surgery.
        Generates 3 patches (Ancilla, Target, Control) and connecting boundaries.

        Layout:
            Ancilla (Top-Left) | Target (Top-Right)
            -------------------|-------------------
            Control (Bot-Left) |      (Empty)
        """

        # Check valid State Initializations
        valid_states = ["+", "-", "+i", "-i", "0", "1"]
        if control_state_init not in valid_states:
            raise ValueError(
                f"Invalid control_state_init: {control_state_init}. "
                f"Valid options are: {valid_states}",
            )

        if target_state_init not in valid_states:
            raise ValueError(
                f"Invalid target_state_init: {target_state_init}. "
                f"Valid options are: {valid_states}",
            )

        # Initialize Parameters
        self.distance = distance
        self.control_state_init = control_state_init
        self.target_state_init = target_state_init
        self.offset_ancilla = 0 + 0j
        self.offset_target = (self.distance * 2) + 0j
        self.offset_control = 0 + (self.distance * 2) * 1j
        self.start_stab_x_ancilla = True
        self.start_stab_x_target = False
        self.start_stab_x_control = False

        # Get Coordinates and Indices
        self.coords = self.get_coords()
        self.q2i = self._get_q2i()
        self.i2q = self._get_i2q()

        # Get Patch Specific Coordinates
        self.coords_ancilla = self.get_coords(specific_coord="ancilla")
        self.coords_target = self.get_coords(specific_coord="target")
        self.coords_control = self.get_coords(specific_coord="control")

        # Getting Data Qubit Indices
        self.anc_data_idx = self._get_specific_indices("DATA", self.coords_ancilla)
        self.target_data_idx = self._get_specific_indices("DATA", self.coords_target)
        self.control_data_idx = self._get_specific_indices("DATA", self.coords_control)

        # Get Stabilizers for ancilla
        self.anc_x_stb_idx = (
            self._get_specific_indices("X-STAB", self.coords_ancilla)
            + self._get_specific_indices("X-STAB-BOUND-A-A", self.coords_ancilla)
            + self._get_specific_indices("X-STAB-BOUND-B-A", self.coords_ancilla)
        )
        self.anc_z_stb_idx = (
            self._get_specific_indices("Z-STAB", self.coords_ancilla)
            + self._get_specific_indices("Z-STAB-BOUND-L-A", self.coords_ancilla)
            + self._get_specific_indices("Z-STAB-BOUND-R-A", self.coords_ancilla)
        )
        self.anc_x_bdy_b_stb_idx = self._get_specific_indices(
            "X-STAB-BOUND-B-A",
            self.coords_ancilla,
        )
        self.anc_z_bdy_r_stb_idx = self._get_specific_indices(
            "Z-STAB-BOUND-R-A",
            self.coords_ancilla,
        )

        # Get Stabilizers for Control
        self.control_x_stb_idx = (
            self._get_specific_indices("X-STAB", self.coords_control)
            + self._get_specific_indices("X-STAB-BOUND-A-C", self.coords_control)
            + self._get_specific_indices("X-STAB-BOUND-B-C", self.coords_control)
        )
        self.control_z_stb_idx = (
            self._get_specific_indices("Z-STAB", self.coords_control)
            + self._get_specific_indices("Z-STAB-BOUND-L-C", self.coords_control)
            + self._get_specific_indices("Z-STAB-BOUND-R-C", self.coords_control)
        )

        # Get Stabilizers for Target
        self.target_x_stb_idx = (
            self._get_specific_indices("X-STAB", self.coords_target)
            + self._get_specific_indices("X-STAB-BOUND-A-T", self.coords_target)
            + self._get_specific_indices("X-STAB-BOUND-B-T", self.coords_target)
        )
        self.target_z_stb_idx = (
            self._get_specific_indices("Z-STAB", self.coords_target)
            + self._get_specific_indices("Z-STAB-BOUND-L-T", self.coords_target)
            + self._get_specific_indices("Z-STAB-BOUND-R-T", self.coords_target)
        )

        # Get Stabilizers for Surgery
        self.surgery_x_m_stb_idx = self._get_specific_indices(
            "X-STAB-SURGERY-M",
            self.coords_surgery,
        )
        self.surgery_x_b_stb_idx = self._get_specific_indices(
            "X-STAB-SURGERY-B",
            self.coords_surgery,
        )
        self.surgery_z_l_stb_idx = self._get_specific_indices(
            "Z-STAB-SURGERY-L",
            self.coords_surgery,
        )
        self.surgery_z_m_stb_idx = self._get_specific_indices(
            "Z-STAB-SURGERY-M",
            self.coords_surgery,
        )

        # Additional Stabilizer Indices Definitions that are needed
        self.control_target_all_stab_idx = (
            self.control_x_stb_idx
            + self.target_x_stb_idx
            + self.control_z_stb_idx
            + self.target_z_stb_idx
        )

        # Using set to avoid double indices
        self.combined_x_stab_idx_filtered = self._get_filtered_x_stabilizers()
        self.combined_x_stab_idx_filtered_AT = self._get_filtered_x_stabilizers(merging_type="AT")

        # Combining all stabilizers for easier reset/ measurement
        self.all_stab_idx = list(
            set(
                self.anc_x_stb_idx
                + self.anc_z_stb_idx
                + self.control_x_stb_idx
                + self.control_z_stb_idx
                + self.target_x_stb_idx
                + self.target_z_stb_idx,
            ),
        )

    def get_coords(self, specific_coord=None) -> dict[complex, str]:
        """
        Returns all qubit coordinates with their labels

        Returns:
            dict[complex, str]
                Dictionary with coordinates as keys and labels as values
        """

        """
        Different Patches are needed, because of different Keywords on 
        identical Coordinates (inside dict.):

        -> X-Stab-Boundary-Above-Control & X-Stab-Boundary-Below-Ancilla f.ex. 
        get Keywords for surgery stabilizers
        """

        # Generate Coordinates for each patch
        self.coords_ancilla = self._get_central_labels(
            self.offset_ancilla,
            self.start_stab_x_ancilla,
        )
        self.coords_target = self._get_central_labels(
            self.offset_target,
            self.start_stab_x_target,
        )
        self.coords_control = self._get_central_labels(
            self.offset_control,
            self.start_stab_x_control,
        )
        self.coords_surgery: dict[complex, str] = {}

        # Add Boundary and Surgery Stabilizers
        self._get_boundary_labels()

        if specific_coord is None:
            # Full Coordinate Dictionary
            # Careful: Overlapping Coordinates with different labels!
            # -> This means that some Coordinates will be overwritten in the dict.
            return_coords = (
                self.coords_ancilla | self.coords_target | self.coords_control | self.coords_surgery
            )

        elif specific_coord == "ancilla":
            # Get only Ancilla Coordinates
            return_coords = self.coords_ancilla
        elif specific_coord == "target":
            # Get only Target Coordinates
            return_coords = self.coords_target
        elif specific_coord == "control":
            # Get only Control Coordinates
            return_coords = self.coords_control
        else:
            raise ValueError("specific_coord must be one of: None, 'ancilla', 'target', 'control'")

        return return_coords

    def _get_filtered_x_stabilizers(self, merging_type=None) -> list[int]:
        """
        Returns the filtered list of X stabilizer indices depending on the merging type
        -> If no merging type is given, returns filtered list for initialization
        -> If AT merging is selected, additional surgery stabilizers are added
           (AC has only additional Z surgery stabilizers)
        """

        if merging_type == "AT":
            combined_x_stab_idx = (
                self.anc_x_stb_idx
                + self.control_x_stb_idx
                + self.target_x_stb_idx
                + self.surgery_x_b_stb_idx
                + self.surgery_x_m_stb_idx
            )
        elif merging_type in {None, "AC"}:
            combined_x_stab_idx = (
                self.anc_x_stb_idx + self.control_x_stb_idx + self.target_x_stb_idx
            )
        else:
            raise ValueError("merging_type must be one of: None, 'AC', 'AT'")

        # Using set to avoid double indices
        filtered_x_stab_idx = list(set(combined_x_stab_idx))

        return filtered_x_stab_idx

    def get_combined_xz_stabs_merging_lattice(
        self,
        merging_type=None,
        split_type=None,
    ) -> list[int]:
        """
        Returns the combined list of X and Z stabilizer coordinates depending on the merging type
        """

        if merging_type == "AC" or split_type == "AC":
            # Adding h gate for X stabilizers only on merging lattices
            # -> Filtering out double coords in big lattice
            combined_x_stab_merging_lattices = self.anc_x_stb_idx + self.control_x_stb_idx
            combined_z_stab_merging_lattices = (
                self.anc_z_stb_idx
                + self.control_z_stb_idx
                + self.surgery_z_l_stb_idx
                + self.surgery_z_m_stb_idx
            )

        elif merging_type == "AT" or split_type == "AT":
            # Adding h gate for X stabilizers only on merging lattices
            # -> Filtering out double coords in big lattice
            combined_x_stab_merging_lattices = (
                self.anc_x_stb_idx
                + self.target_x_stb_idx
                + self.surgery_x_b_stb_idx
                + self.surgery_x_m_stb_idx
            )
            combined_z_stab_merging_lattices = self.anc_z_stb_idx + self.target_z_stb_idx

        else:
            raise ValueError("merging_type must be one of: 'AC', 'AT'")

        # Using set to avoid double indices
        filtered_x_stab_idx = list(set(combined_x_stab_merging_lattices))
        filtered_z_stab_idx = list(set(combined_z_stab_merging_lattices))

        return filtered_x_stab_idx, filtered_z_stab_idx

    def get_logical_strings(
        self,
        shift_cx_for_y: bool = False,
        shift_tz_for_y: bool = False,
        shift_tx_for_y: bool = False,
        shift_cz_for_y: bool = False,
    ) -> dict[str, list[int]]:
        """
        Args:
            q2i: Mapping from complex coordinates to qubit indices
            distance: Code distance

        Returns:
            Dictionary with keys:
                - 'a_z': Ancilla Z-logical indices
                - 't_x': Target X-logical indices
                - 't_z': Target Z-logical indices
                - 't_y': Target Y-logical indices
                - 'c_x': Control X-logical indices
                - 'c_z': Control Z-logical indices
                - 'c_y': Control Y-logical indices
        """

        # Control Y logical observable components
        # Y is at bottom right corner of control region: (distance*2-1, distance*4-1)
        c_y_z_string = [
            self.q2i[real + (self.distance * 4 - 1) * 1j]
            for real in range(1, self.distance * 2 - 1, 2)
        ]
        c_y_corner = [self.q2i[self.distance * 2 - 1 + (self.distance * 4 - 1) * 1j]]
        c_y_x_string = [
            self.q2i[(self.distance * 2 - 1) + imag * 1j]
            for imag in range(self.distance * 2 + 1, self.distance * 4 - 1, 2)
        ]

        # Target Y logical observable components
        # Y is at bottom right corner of target region: (distance*4-1, distance*2-1)
        t_y_z_string = [
            self.q2i[real + (self.distance * 2 - 1) * 1j]
            for real in range(self.distance * 2 + 1, self.distance * 4 - 1, 2)
        ]
        t_y_corner = [self.q2i[self.distance * 4 - 1 + (self.distance * 2 - 1) * 1j]]
        t_y_x_string = [
            self.q2i[(self.distance * 4 - 1) + imag * 1j]
            for imag in range(1, self.distance * 2 - 1, 2)
        ]

        # Offsets for conditional shifts in Y-including flows (symmetric for both directions)
        cx_shift_real = (self.distance * 2 - 1) if shift_cx_for_y else 1
        tz_shift_imag = (self.distance * 2 - 1) if shift_tz_for_y else 1
        tx_shift_real = (self.distance * 4 - 1) if shift_tx_for_y else (self.distance * 2 + 1)
        cz_shift_imag = (self.distance * 4 - 1) if shift_cz_for_y else (self.distance * 2 + 1)

        return {
            "a_z": [self.q2i[real + 1j] for real in range(1, self.distance * 2, 2)],
            "t_x": [self.q2i[tx_shift_real + imag * 1j] for imag in range(1, self.distance * 2, 2)],
            "t_z": [
                self.q2i[real + tz_shift_imag * 1j]
                for real in range(self.distance * 2 + 1, self.distance * 4, 2)
            ],
            "t_y": {
                "z_string": t_y_z_string,
                "y_corner": t_y_corner,
                "x_string": t_y_x_string,
            },
            "c_x": [
                self.q2i[cx_shift_real + imag * 1j]
                for imag in range(self.distance * 2 + 1, self.distance * 4, 2)
            ],
            "c_z": [self.q2i[real + cz_shift_imag * 1j] for real in range(1, self.distance * 2, 2)],
            "c_y": {
                "z_string": c_y_z_string,
                "y_corner": c_y_corner,
                "x_string": c_y_x_string,
            },
        }

    def _get_central_labels(
        self,
        offset: complex,
        starting_stabilizer_x: bool,
    ) -> dict[complex, str]:
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

    def _get_boundary_labels(self) -> dict[complex, str]:
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
