import stim

from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry

__all__ = ["SurgeryPauliObservables"]


class SurgeryPauliObservables:
    def __init__(
        self,
        geometry: SurgeryGeometry,
        control_measure_basis: str,
        target_measure_basis: str,
        control_state_init: str,
        target_state_init: str,
        flow_type: str,
    ):
        # Check validity of input arguments
        if flow_type not in {"incoming_flow", "outgoing_flow"}:
            raise ValueError("Invalid type. Must be 'incoming_flow' or 'outgoing_flow'.")

        # Preliminary Setup
        self.geometry = geometry
        self.all_logical_strings = self.geometry.get_logical_strings()
        self.all_logical_strings_shifted = self.geometry.get_logical_strings(
            shift_cx_for_y=True,
            shift_tz_for_y=True,
            shift_cz_for_y=True,
            shift_tx_for_y=True,
        )

        # Select Basis depending on incoming or outgoing flow
        if flow_type == "outgoing_flow":
            self.control_measure_basis = control_measure_basis
            self.target_measure_basis = target_measure_basis
        elif flow_type == "incoming_flow":
            self.control_measure_basis = control_state_init[0]
            self.target_measure_basis = target_state_init[0]

    def build_circuit(self) -> stim.Circuit:
        # Init return Circuit
        return_circuit = stim.Circuit()

        # Adding Pauli Observables for both control and target patches
        return_circuit += self._add_pauli_observables(
            patch="control",
        )
        return_circuit += self._add_pauli_observables(
            patch="target",
        )

        return return_circuit

    def _add_pauli_observables(self, patch: str) -> stim.Circuit:
        # Check validity of input arguments
        if patch not in {"control", "target"}:
            raise ValueError("Invalid patch type. Must be 'control' or 'target'.")

        # Initilize measurement Circuit
        measurement_circuit = stim.Circuit()

        # Select Patch attributes
        if patch == "control":
            measure_basis = self.control_measure_basis
            x_string = self.all_logical_strings_shifted["c_x"]
            z_string = self.all_logical_strings_shifted["c_z"]
            y_logical_string = [
                self.all_logical_strings["c_y"]["x_string"],
                self.all_logical_strings["c_y"]["y_corner"],
                self.all_logical_strings["c_y"]["z_string"],
            ]

        elif patch == "target":
            measure_basis = self.target_measure_basis
            x_string = self.all_logical_strings_shifted["t_x"]
            z_string = self.all_logical_strings_shifted["t_z"]
            y_logical_string = [
                self.all_logical_strings["t_y"]["x_string"],
                self.all_logical_strings["t_y"]["y_corner"],
                self.all_logical_strings["t_y"]["z_string"],
            ]

        # Adding Measurements and Observables based on measurement basis and observable type
        if measure_basis == "X":
            measurement_circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_pauli(i, "X") for i in x_string],
                0,
            )
        elif measure_basis == "Z":
            measurement_circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_pauli(i, "Z") for i in z_string],
                0,
            )
        elif measure_basis == "Y":
            measurement_circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_pauli(i, "X") for i in y_logical_string[0]]
                + [stim.target_pauli(i, "Y") for i in y_logical_string[1]]
                + [stim.target_pauli(i, "Z") for i in y_logical_string[2]],
                0,
            )

        # Return Circuit
        return measurement_circuit
