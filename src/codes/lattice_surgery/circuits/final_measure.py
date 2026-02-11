import stim

from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry

__all__ = ["SurgeryFinalMeasure"]


class SurgeryFinalMeasure:
    def __init__(
        self,
        geometry: SurgeryGeometry,
        curr_flow: str,
        valid_flow: bool,
        control_measure_basis: str,
        target_measure_basis: str,
    ):
        # Preliminary Setup
        self.geometry = geometry
        self.valid_flow = valid_flow
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

        # Check if flow is valid or non deterministic
        if self.valid_flow:
            # If flow is valid, we can apply measurements and observables based on the flow
            return_circuit += self._add_measurements_and_observables(
                patch="control",
                observable_type="record",
            )
            return_circuit += self._add_measurements_and_observables(
                patch="target",
                observable_type="record",
            )
        else:
            return_circuit += self._add_measurements_and_observables(
                patch="control",
                observable_type="pauli",
            )
            return_circuit += self._add_measurements_and_observables(
                patch="target",
                observable_type="pauli",
            )

        return return_circuit

    def _add_measurements_and_observables(self, patch: str, observable_type: str) -> stim.Circuit:
        # Check validity of input arguments
        if observable_type not in {"pauli", "record"}:
            raise ValueError("Invalid observable type. Must be 'pauli' or 'record'.")
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
            measurement_circuit.append("TICK")
            measurement_circuit.append("MX", x_string)
            if observable_type == "record":
                measurement_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-len(x_string) + k) for k in range(len(x_string))],
                    0,
                )
            else:
                measurement_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_pauli(i, "X") for i in x_string],
                    0,
                )
        elif measure_basis == "Z":
            measurement_circuit.append("TICK")
            measurement_circuit.append("MZ", z_string)
            if observable_type == "record":
                measurement_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-len(z_string) + k) for k in range(len(z_string))],
                    0,
                )
            else:
                measurement_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_pauli(i, "Z") for i in z_string],
                    0,
                )
        elif measure_basis == "Y":
            # Apply Y measurement on patch
            measurement_circuit.append("TICK")
            measurement_circuit.append("MX", y_logical_string[0])
            measurement_circuit.append("MY", y_logical_string[1])
            measurement_circuit.append("MZ", y_logical_string[2])

            if observable_type == "record":
                total_measurements = (
                    len(y_logical_string[0]) + len(y_logical_string[1]) + len(y_logical_string[2])
                )
                measurement_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-total_measurements + k) for k in range(total_measurements)],
                    0,
                )
            else:
                measurement_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_pauli(i, "X") for i in y_logical_string[0]]
                    + [stim.target_pauli(i, "Y") for i in y_logical_string[1]]
                    + [stim.target_pauli(i, "Z") for i in y_logical_string[2]],
                    0,
                )

        # Return Circuit
        return measurement_circuit
