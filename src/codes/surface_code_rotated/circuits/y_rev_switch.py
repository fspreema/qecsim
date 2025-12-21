import stim

from src.codes.surface_code_rotated.data_geometry import MasterGeometry, MasterPairings

__all__ = ["YRevSwitchCircuit"]


class YRevSwitchCircuit:
    def __init__(
        self,
        master_geometry: MasterGeometry,
        master_pairings: MasterPairings,
    ):
        """
        Initialize the Y-Basis Reversal Switch Circuit
        -> Follows the exact procedure as the Y-Basis Switch but in reverse
            - 1) First cx and xcy gate applications
            - 2) Then half digaonal h gates and SQRT_X_DAG gates

        Parameters:
            geometry : SurfaceGeometry
            pairings : SurfacePairings
        """

        # Set Geometry and Pairings for Y-Basis Switch
        self.geometry = master_geometry.geometry_ybasis
        self.stab_to_data_switch = master_pairings.pairings_yswitch.stab_to_data
        self.stab_to_data_xcy = master_pairings.pairings_yswitch.stab_to_data_xcy

    def build_circuit(self) -> stim.Circuit:
        # Initialize Empty Circuit
        self.switch_circuit = stim.Circuit()

        # Building Circuit
        self.switch_circuit += self._apply_y_rev_switch()

        return self.switch_circuit

    def _apply_y_rev_switch(self) -> stim.Circuit:
        ##########################
        # Applying Y-Basis Reversal Switch
        ##########################
        rev_switch_circuit = stim.Circuit()

        rev_switch_circuit.append(
            "R",
            self.geometry.stab_idx + self.geometry.right_h + self.geometry.upper_h,
        )
        rev_switch_circuit.append("TICK")
        rev_switch_circuit.append("H", self.geometry.stab_switch_apply_h)
        rev_switch_circuit.append("TICK")

        #########################################
        # Adding the XCY gates after the H switch
        #########################################

        # 2) CX Operations
        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "5TICK":
                if len(coord_pairs) == 2:
                    index_pairs = []
                    index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                    index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                    rev_switch_circuit.append("CX", index_pairs)
                else:
                    index_pairs = []
                    index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                    index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                    rev_switch_circuit.append("CX", index_pairs)

        rev_switch_circuit.append("TICK")

        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "4TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                rev_switch_circuit.append("CX", index_pairs)

        rev_switch_circuit.append("TICK")

        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "3.5TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                rev_switch_circuit.append("XCY", index_pairs)

        rev_switch_circuit.append("TICK")

        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "3TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                rev_switch_circuit.append("CX", index_pairs)

        rev_switch_circuit.append("TICK")

        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "2TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                rev_switch_circuit.append("CX", index_pairs)

        rev_switch_circuit.append("TICK")

        for coord_pairs, order in self.stab_to_data_xcy.items():
            # Parallel Implementation of CX
            if order == "1TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                rev_switch_circuit.append("XCY", index_pairs)

        rev_switch_circuit.append("TICK")

        # -------Continue-Circuit------------

        # Adding stabs h
        rev_switch_circuit.append("H", self.geometry._y_basis_get_switch_h_qubits())
        rev_switch_circuit.append("SQRT_X_DAG", self.geometry._y_basis_get_switch_xdag_qubits())
        rev_switch_circuit.append("TICK")

        # Adding half diagonal H
        rev_switch_circuit.append("H", self.geometry.stab_x_idx + self.geometry.right_h)
        rev_switch_circuit.append("TICK")

        # Adding Resets
        rev_switch_circuit.append(
            "MZ",
            self.geometry.stab_idx + self.geometry.right_h + self.geometry.upper_h,
        )
        rev_switch_circuit.append("MY", self.geometry.y_index)
        rev_switch_circuit.append("SHIFT_COORDS", arg=(0, 0, 1))
        rev_switch_circuit.append("TICK")

        return rev_switch_circuit
