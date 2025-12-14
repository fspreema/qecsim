import stim

"""
Class which is used to build the stim sub circuits

-> Used to mark which detectors are used and what measurements are included
-> Allows mapping of detectors which later get used for filtration
"""


class GenCircuit:
    def __init__(self, circuit: stim.Circuit, auto_detectors: bool = True):
        self.circuit = circuit
        self.auto_detectors = auto_detectors

    def __add__(self, other: "GenCircuit") -> "GenCircuit":
        return self.add_circuit(other, add_overlapping_det=True)

    def add_circuit(self, other: "GenCircuit", add_overlapping_det: bool) -> "GenCircuit":
        """
        If two GenCircuits are added together, return a new GenCircuit

        Also handle automatic Detector Placement between the two circuits if selected
        -> If Measurement on same Index happens on both circuits
            1. Find Both measurement record indices
            2. Add a detector after both circuits with these record targets
        """
        new_circuit = stim.Circuit()
        new_circuit += self.circuit
        new_circuit += other.circuit

        # Checking for overlapping measurement indices
        if add_overlapping_det and self.auto_detectors and other.auto_detectors:
            overlapping_meas = self.find_overlapping_measurements(other)

            for self_rec_index, other_rec_index in overlapping_meas:
                # Other target is the newest edition to the circuit
                other_target = stim.target_rec(
                    -(other.circuit.num_measurements - other_rec_index + 1),
                )

                # Self target is the older measurement and therefore needs to be shifted back
                self_target = stim.target_rec(
                    -(
                        other.circuit.num_measurements
                        + self.circuit.num_measurements
                        - self_rec_index
                        + 1
                    ),
                )

                new_circuit.append(
                    "DETECTOR",
                    [self_target, other_target],
                )

        return GenCircuit(new_circuit)

    def add_operator(self, op: str, qubits: list[int]) -> None:
        # Add an operator to the circuit

        self.circuit.append(op, qubits)

    def add_coords(self, qubit: int, coord: complex) -> None:
        # Add coordinates to the circuit
        self.circuit.append("QUBIT_COORDS", [qubit], [coord.real, coord.imag])

    def add_measurement(self, qubits: list[int]) -> None:
        # Add a measurement to the circuit
        self.circuit.append("M", qubits)

    def add_detector(self, args: tuple[float, float, float], qubit_records: list[int]) -> None:
        """
        This method is only used at the beginning when the automatic detectors can not be used
        -> First detectors are only assoicated with one measurement instead of a pair
        """
        # Add Detector to circuit
        self.circuit.append("DETECTOR", qubit_records, args)

    def find_overlapping_measurements(self, other: "GenCircuit") -> list[tuple[int, int]]:
        """
        Find overlapping measurement record indices between two circuits

        -> Returns list of tuples with (self_rec_index, other_rec_index)
        """
        self_meas_map: dict[int, int] = {}
        rec_idx = 0
        for instruction in self.circuit:
            if instruction.name in {"M", "MX", "MY", "MZ", "MR"}:
                for target in instruction.targets_copy():
                    rec_idx += 1
                    # stim.target_value is for qubit index
                    self_meas_map[target.value] = rec_idx

        other_meas_map: dict[int, int] = {}
        rec_idx = 0
        for instruction in other.circuit:
            if instruction.name in {"M", "MX", "MY", "MZ", "MR"}:
                for target in instruction.targets_copy():
                    rec_idx += 1
                    other_meas_map[target.value] = rec_idx

        matched_indices: list[tuple[int, int]] = []
        for qubit, self_rec in self_meas_map.items():
            if qubit in other_meas_map:
                other_rec = other_meas_map[qubit]
                matched_indices.append((self_rec, other_rec))

        return matched_indices

    def get_circuit(self) -> stim.Circuit:
        """
        Returns the built stim Circuit currently stored in the class

        -> Before the circuit is returned, TICKS are added to ensure no overlapping
        """

        finished_circ = stim.Circuit()

        # Filter out any accidental TICKS or SHIFT_COORDS already in the circuit
        double_tick = False
        double_shift = False

        for instruction in self.circuit:
            if instruction.name == "TICK" and not double_tick:
                finished_circ.append(instruction)
                double_tick = True
            elif instruction.name == "SHIFT_COORDS" and not double_shift:
                finished_circ.append(instruction)
                double_shift = True
            else:
                double_tick = False
                double_shift = False

        return finished_circ

    def get_flagged_detectors(self) -> dict:
        """
        Returns the dictionary of flagged detectors needed for post selection
        """

        pass
