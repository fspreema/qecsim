import stim

"""
Class which is used to build the stim sub circuits

-> Used to mark which detectors are used and what measurements are included
-> Allows mapping of detectors which later get used for filtration
"""


class GenCircuit:
    def __init__(self, circuit: stim.Circuit):
        self.circuit = circuit

    def __add__(self, other: "GenCircuit", add_overlapping_det: bool) -> "GenCircuit":
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
        if add_overlapping_det:
            overlapping_meas = self.find_overlapping_measurements(other)

            for self_rec_index, other_rec_index in overlapping_meas:
                new_circuit.append(
                    "DETECTOR",
                    [self_rec_index, other_rec_index - self.circuit.num_measurements],
                )

        return GenCircuit(new_circuit)

    def add_operator(self, op: str, qubits: list[int]) -> None:
        # Add an operator to the circuit
        self.circuit.append(op, qubits)

    def add_measurement(self, qubits: list[int]) -> None:
        # Add a measurement to the circuit
        self.circuit.append("M", qubits)

    def add_detector(
        self,
        args: tuple[float, float, float],
        qubit_records: list[int],
        flag: bool,
    ) -> None:
        """
        Flag is used to mark if detector result is needed for filtration
        """

        pass

    def find_overlapping_measurements(self, other: "GenCircuit") -> list[tuple[int, int]]:
        """
        Find overlapping measurement record indices between two circuits

        -> Returns list of tuples with (self_rec_index, other_rec_index)
        """
        pass

    def get_circuit(self) -> stim.Circuit:
        """
        Returns the built stim Circuit currently stored in the class

        -> Before the cricuit is returned, TICKS are added to ensure no overlappint
        """

        finished_circ = stim.Circuit()

        for instruction in self.circuit:
            if instruction.name in {
                "QUBIT_COORDS",
                "DETECTOR",
                "OBSERVABLE_INCLUDE",
                "SHIFT_COORDS",
            }:
                continue
            else:
                finished_circ.append(instruction)
                finished_circ.append("TICK")

        return finished_circ

    def get_flagged_detectors(self) -> dict:
        """
        Returns the dictionary of flagged detectors needed for post selection
        """

        pass
