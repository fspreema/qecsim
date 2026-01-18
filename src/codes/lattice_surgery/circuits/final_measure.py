import stim

from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry

__all__ = ["SurgeryFinalMeasure"]


class SurgeryFinalMeasure:
    def __init__(
        self,
        geometry: SurgeryGeometry,
        flow: str,
    ):
        # Preliminary Setup
        self.geometry = geometry
        self.flow = flow
        self.all_logical_strings = self.geometry.get_logical_strings()
        self.all_logical_strings_shifted = self.geometry.get_logical_strings(
            shift_cx_for_y=True,
            shift_tz_for_y=True,
            shift_cz_for_y=True,
            shift_tx_for_y=True,
        )

    def build_circuit(self) -> stim.Circuit:
        # Init return Circuit
        return_circuit = stim.Circuit()

        # Build Final Measurement Circuit
        is_y_flow = "Y" in self.flow

        if is_y_flow:
            pass
        else:
            return_circuit += self._apply_non_y_measurements()
            return_circuit += self._apply_logical_observables()

        return return_circuit

    def _apply_y_measurements(self):
        # Define central flow dictionary
        y_flow_measurements = {
            "YZ -> XY": [("control", "X", True), ("target", "Y", False)],
            "YI -> YX": [("control", "Y", False), ("target", "X", False)],
            "YX -> YI": [("control", "Y", False)],
            "YY -> XZ": [("control", "X", False), ("target", "Z", False)],
            "IY -> ZY": [("control", "Z", False), ("target", "Y", False)],
            "XY -> YZ": [("control", "Y", False), ("target", "Z", False)],
        }

    def _meassure_patch_operator(self, patch: str, op_type: str, shifted: bool) -> stim.Circuit:
        """
        Returns a circuit that measures the logical operator of a given patch.
        -> Measure exactly along the logical strings of the patch.
        -> Used for Y basis flows.

        Arguments:
            patch (str): Which patch to measure ('control' or 'target')
            op_type (str): Which operator to measure ('X', 'Y', or 'Z')
            shifted (bool): Whether to use the shifted logical strings or not
        """

        # Init Measure Circuit
        measure_circuit = stim.Circuit()
        measure_circuit.append("TICK")

        if patch == "control":
            pass

    def _validate_non_y_flow_selection(self, flow: str, patch: str):
        # Creating Dictionary of valid flows in combination with init states
        valid_flows_control = {
            "+": {"XI -> XX", "XX -> XI", "IX -> IX"},
            "-": {"XI -> XX", "XX -> XI", "IX -> IX"},
            "0": {"IZ -> ZZ", "ZZ -> IZ", "ZI -> ZI", "ZX -> ZX", "IX -> IX"},
            "1": {"IZ -> ZZ", "ZZ -> IZ", "ZI -> ZI", "ZX -> ZX", "IX -> IX"},
        }
        valid_flows_target = {
            "+": {"XI -> XX", "XX -> XI", "IX -> IX", "ZI -> ZI", "ZX -> ZX"},
            "-": {"XI -> XX", "XX -> XI", "IX -> IX", "ZI -> ZI", "ZX -> ZX"},
            "0": {"IZ -> ZZ", "ZZ -> IZ", "ZI -> ZI"},
            "1": {"IZ -> ZZ", "ZZ -> IZ", "ZI -> ZI"},
        }

        if patch == "control" and flow in valid_flows_control.get(
            self.geometry.control_state_init,
            set(),
        ):
            return True
        elif patch == "target" and flow in valid_flows_target.get(
            self.geometry.target_state_init,
            set(),
        ):
            return True
        else:
            return False

    def _get_non_y_measurements_for_flow(self, flow: str):
        # Creating Dictionary of measurement operations per flow
        flow_measurements = {
            "XI -> XX": ["c_x", "t_x"],
            "XX -> XI": ["c_x"],
            "IX -> IX": ["t_x"],
            "IZ -> ZZ": ["c_z", "t_z"],
            "ZZ -> IZ": ["c_z"],
            "ZI -> ZI": ["t_z"],
            "ZX -> ZX": ["c_z", "t_x"],
        }

        return flow_measurements.get(flow, [])

    def _apply_non_y_measurements(self):
        # Init measure Circuit
        measure_circuit = stim.Circuit()

        if self.geometry.control_state_init in {"+", "-"}:
            measure_circuit.append("MX", self.geometry.control_data_idx)

        elif self.geometry.control_state_init in {"0", "1"}:
            measure_circuit.append("MZ", self.geometry.control_data_idx)

        if self.geometry.target_state_init in {"+", "-"}:
            measure_circuit.append("MX", self.geometry.target_data_idx)

        elif self.geometry.target_state_init in {"0", "1"}:
            measure_circuit.append("MZ", self.geometry.target_data_idx)

        return measure_circuit

    def _apply_logical_observables(self):
        # Init Observable Circuit
        observable_circuit = stim.Circuit()

        # Checking validity of flow selection
        if not self._validate_non_y_flow_selection(self.flow, "control"):
            raise ValueError("Invalid flow selected for control patch initial state!")
        if not self._validate_non_y_flow_selection(self.flow, "target"):
            raise ValueError("Invalid flow selected for target patch initial state!")

        # Get indices of logical operator (can be on both patches) depending on flow
        flow_measurements = self._get_non_y_measurements_for_flow(self.flow)

        logical_string = [self.all_logical_strings[i] for i in flow_measurements]

        tar_rec = []

        for rec_pos, index in enumerate(
            self.geometry.control_data_idx + self.geometry.target_data_idx,
        ):
            if index in logical_string:
                tar_rec.append(rec_pos)

        observable_circuit.append(
            "OBSERVABLE_INCLUDE",
            [
                stim.target_rec(
                    -len(self.geometry.control_data_idx + self.geometry.target_data_idx) + k,
                )
                for k in tar_rec
            ],
            0,
        )

        return observable_circuit
