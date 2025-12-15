import stim

from src.codes.xzzx.get_stab_pairings import XZZXPairings
from src.codes.xzzx.xzzx_geom import XZZXGeometry
from src.core.cx_builder import cx_builder

__all__ = ["RepetitionCircuit"]


class RepetitionCircuit:
    def __init__(
        self,
        geometry: XZZXGeometry,
        stab_pairings: XZZXPairings,
    ):
        # Setting up parameters
        self.geometry = geometry
        self.stab_pairings = stab_pairings

    def build_repetition_circuit(self) -> stim.Circuit:
        self.repetition_circuit = stim.Circuit()

        # Building Circuit
        self.repetition_circuit += self._build_repetition_circuit()
        self.repetition_circuit += self._get_detectors()

        # Return repeated Circuit
        return self.repetition_circuit * (self.geometry.distance - 1)

    def _build_repetition_circuit(self):
        ###########################################
        # Adding Repeat Block
        ###########################################

        repeat_circuit = stim.Circuit()

        repeat_circuit.append("SHIFT_COORDS", arg=(0, 0, 1))
        repeat_circuit.append("R", self.geometry.stab_idx)

        repeat_circuit.append("TICK")
        repeat_circuit.append("H", self.geometry.stab_idx)

        repeat_circuit.append("TICK")

        ####################################################
        # CX Operations
        ####################################################

        cx_builder(
            circuit=repeat_circuit,
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_pairings.stab_to_data,
            orders=("1-CX", "2-CZ", "3-CZ", "4-CX"),
        )

        # Retreive Boundary + Normal Stabilizers Ancilla:
        repeat_circuit.append("H", self.geometry.stab_idx)
        repeat_circuit.append("TICK")

        repeat_circuit.append("M", self.geometry.stab_idx)
        repeat_circuit.append("TICK")

        return repeat_circuit

    def _get_detectors(self):
        # Initlizing Detector Circuit
        det_circuit = stim.Circuit()

        ######################################################
        # Implementing Detectors (All Basis are deterministic)
        ######################################################

        # Determining Position in the measurement Run of only the Ancilla
        pos_to_index_ver: list = []
        pos_to_index_hor: list = []

        for pos, index in enumerate(self.geometry.stab_idx):
            if index in self.geometry.stab_ver_idx:
                pos_to_index_ver.append([pos, index])
            elif index in self.geometry.stab_hor_idx:
                pos_to_index_hor.append([pos, index])

        """
        As intial circuit run is completed, now we define all stabilizers in every basis
        -> Detectors on both basis are now deterministic
        """

        # Adding the needed Detectors (Z-Basis)
        for index_pos in pos_to_index_ver:
            current_tar = index_pos[0] - len(self.geometry.stab_idx)
            previous_tar = index_pos[0] - 2 * len(self.geometry.stab_idx)
            q_index = index_pos[1]
            det_circuit.append(
                "DETECTOR",
                [stim.target_rec(current_tar), stim.target_rec(previous_tar)],
                (self.geometry.i2q[q_index].real, self.geometry.i2q[q_index].imag, 0),
            )

        # Adding the needed Detectors (X-Basis)
        for index_pos in pos_to_index_hor:
            current_tar = index_pos[0] - len(self.geometry.stab_idx)
            previous_tar = index_pos[0] - 2 * len(self.geometry.stab_idx)
            q_index = index_pos[1]
            det_circuit.append(
                "DETECTOR",
                [stim.target_rec(current_tar), stim.target_rec(previous_tar)],
                (self.geometry.i2q[q_index].real, self.geometry.i2q[q_index].imag, 0),
            )

        return det_circuit
