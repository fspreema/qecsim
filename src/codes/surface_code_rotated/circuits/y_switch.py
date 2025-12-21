import stim

from src.codes.surface_code_rotated.data_geometry import MasterGeometry, MasterPairings
from src.core.cx_builder import cx_builder

__all__ = ["YSwitchCircuit"]


class YSwitchCircuit:
    def __init__(
        self,
        master_geometry: MasterGeometry,
        master_pairings: MasterPairings,
    ):
        """
        Initialize the Y-Switch Circuit

        The following procedure will be done
        -> Along the mirrored diagonal one half h gate one half not
        -> Boundary on the site of the h gate gets expanded
            -> Every postion is now a boundary! (No 2 coords distance between them)
        -> Along the Y digaonal we do SQRT_X_DAG on all ancilla qubits
        """

        # Set Geometry and Pairings for Y-Basis Switch
        self.geometry = master_geometry.geometry_ybasis
        self.stab_to_data_y_basis = master_pairings.pairings_ybasis.stab_to_data
        self.stab_to_data_switch = master_pairings.pairings_yswitch.stab_to_data
        self.stab_to_data_xcy = master_pairings.pairings_yswitch.stab_to_data_xcy

    def build_circuit(self) -> stim.Circuit:
        # Initialize Empty Circuit
        self.switch_circuit = stim.Circuit()

        # Building Circuit
        self.switch_circuit += self._apply_pre_switch_circ()
        self.switch_circuit += self._apply_y_switch()
        self.switch_circuit += self._add_logical_observables()

        return self.switch_circuit

    def _apply_pre_switch_circ(self) -> stim.Circuit:
        # Initialize Pre-Switch Circuit
        pre_switch_circ = stim.Circuit()

        # 1) Reset/ Basis
        pre_switch_circ.append("R", self.geometry.stab_idx)
        pre_switch_circ.append("TICK")
        pre_switch_circ.append("H", self.geometry.stab_x_idx)
        pre_switch_circ.append("TICK")

        # 2) CX Operations
        # Adding all the CX gates
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_to_data_y_basis,
            circuit=pre_switch_circ,
            excluded_index=self.geometry.y_index,
        )

        # 3) Basis/ Measurement
        pre_switch_circ.append("H", self.geometry.stab_x_idx)
        pre_switch_circ.append("TICK")
        pre_switch_circ.append("MZ", self.geometry.stab_idx)
        pre_switch_circ.append("SHIFT_COORDS", arg=(0, 0, 1))
        pre_switch_circ.append("TICK")

        return pre_switch_circ

    def _apply_y_switch(self) -> stim.Circuit:
        # Initialize Switch Circuit
        switch_circ = stim.Circuit()

        # ---------Adding-RY-Gate----------------
        switch_circ.append("RY", self.geometry.y_index)

        # ---------Adding-Resets-Old&New-Stabs---
        switch_circ.append(
            "RZ",
            self.geometry.stab_idx + self.geometry.right_h + self.geometry.upper_h,
        )

        # 1) Reset/ Basis
        switch_circ.append("TICK")
        switch_circ.append("H", self.geometry.stab_x_idx + self.geometry.right_h)

        # -----------Adding-H-&-X-DAG-Gates------
        switch_circ.append("TICK")
        switch_circ.append("H", self.geometry._y_basis_get_switch_h_qubits())
        switch_circ.append("SQRT_X_DAG", self.geometry._y_basis_get_switch_xdag_qubits())
        switch_circ.append("TICK")

        #########################################
        # Adding the XCY gates after the H switch
        #########################################

        # 2) CX Operations
        for coord_pairs, order in self.stab_to_data_xcy.items():
            # Parallel Implementation of CX
            if order == "1TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                switch_circ.append("XCY", index_pairs)

        switch_circ.append("TICK")

        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "2TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                switch_circ.append("CX", index_pairs)

        switch_circ.append("TICK")

        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "3TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                switch_circ.append("CX", index_pairs)

        switch_circ.append("TICK")

        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "3.5TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                switch_circ.append("XCY", index_pairs)

        switch_circ.append("TICK")

        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "4TICK":
                index_pairs = []
                index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                switch_circ.append("CX", index_pairs)

        switch_circ.append("TICK")

        for coord_pairs, order in self.stab_to_data_switch.items():
            # Parallel Implementation of CX
            if order == "5TICK":
                if len(coord_pairs) == 2:
                    index_pairs = []
                    index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                    index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                    switch_circ.append("CX", index_pairs)
                else:
                    index_pairs = []
                    index_pairs.append(self.geometry.q2i[coord_pairs[1]])
                    index_pairs.append(self.geometry.q2i[coord_pairs[0]])
                    switch_circ.append("CX", index_pairs)

        switch_circ.append("TICK")

        # 3) Basis/ Measurement
        switch_circ.append("H", self.geometry.stab_switch_apply_h)
        switch_circ.append("TICK")
        switch_circ.append(
            "M",
            self.geometry.stab_idx + self.geometry.right_h + self.geometry.upper_h,
        )
        switch_circ.append("SHIFT_COORDS", arg=(0, 0, 1))
        switch_circ.append("TICK")

        return switch_circ

    def _add_logical_observables(self) -> stim.Circuit:
        # Initialize Observable Circuit
        observable_circuit = stim.Circuit()

        ############################################################
        # Adding Observable Includes if logical basis differs from Y
        ############################################################

        """
        We need to remove the added Pauli measurement from the end of the circuit
        -> Else the X/Z paulis tring would anticommute with the RZ/RX reset of the data
        """

        # Getting logical Index
        log_indx = self.geometry.get_logical_observables(
            logical_observable=self.geometry.obs,
            fixed_coord=self.geometry.distance * 2 - 1,
        )

        # XORing the observable away
        if self.geometry.obs == "Z":
            observable_circuit.append("OBSERVABLE_INCLUDE", [f"Z{index}" for index in log_indx], 0)

        elif self.geometry.obs == "X":
            observable_circuit.append("OBSERVABLE_INCLUDE", [f"X{index}" for index in log_indx], 0)

        return observable_circuit
