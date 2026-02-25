import stim

from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry

__all__ = ["SurgeryFinalMeasure"]


class SurgeryFinalMeasure:
    def __init__(
        self,
        geometry: SurgeryGeometry,
        curr_flow: str,
        control_measure_basis: str,
        target_measure_basis: str,
    ):
        # Preliminary Setup
        self.geometry = geometry
        self.curr_flow = curr_flow
        self.control_measure_basis = control_measure_basis
        self.target_measure_basis = target_measure_basis
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

        # Apply Measurements for both control and target patches
        return_circuit += self._add_measurements_and_observables(
            patch="control",
        )
        return_circuit += self._add_measurements_and_observables(
            patch="target",
        )

        return return_circuit

    def _add_measurements_and_observables(self, patch: str) -> stim.Circuit:
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
        measured_qubits = []
        measurement_circuit.append("TICK")

        if measure_basis == "X":
            measurement_circuit.append("MX", x_string)
            measured_qubits = x_string

        elif measure_basis in {"Z", "I"}:
            measurement_circuit.append("MZ", z_string)
            measured_qubits = z_string

        elif measure_basis == "Y":
            # Apply Y measurement on patch
            measurement_circuit.append("MX", y_logical_string[0])
            measurement_circuit.append("MY", y_logical_string[1])
            measurement_circuit.append("MZ", y_logical_string[2])
            measured_qubits = y_logical_string[0] + y_logical_string[1] + y_logical_string[2]

        # Only add Observable if no Identity is selected
        if measure_basis != "I":
            measurement_circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_rec(-len(measured_qubits) + k) for k in range(len(measured_qubits))],
                0,
            )

        # Return Circuit
        return measurement_circuit
