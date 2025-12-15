import stim

from src.codes.xzzx.get_stab_pairings import XZZXPairings
from src.codes.xzzx.xzzx_geom import XZZXGeometry
from src.core.cx_builder import cx_builder

__all__ = ["InitialCircuit"]


class InitialCircuit:
    def __init__(
        self,
        geometry: XZZXGeometry,
        stab_pairings: XZZXPairings,
    ):
        self.geometry = geometry
        self.stab_pairings = stab_pairings

    def build_initial_circuit(self) -> stim.Circuit:
        self.initial_circuit = stim.Circuit()

        self.initial_circuit += self._build_initial_circuit()
        self.initial_circuit += self._get_detectors()

        return self.initial_circuit

    def _build_initial_circuit(self):
        ##########################
        # Defining Initial Circuit
        ##########################

        initial_circuit = stim.Circuit()

        initial_circuit.append("TICK")
        initial_circuit.append("R", self.geometry.stab_idx)

        initial_circuit.append("TICK")
        initial_circuit.append("H", self.geometry.stab_idx)

        initial_circuit.append("TICK")

        ###############
        # CX Operations
        ###############

        cx_builder(
            circuit=initial_circuit,
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_pairings.stab_to_data,
            orders=("1-CX", "2-CZ", "3-CZ", "4-CX"),
        )

        # Retreive Boundary + Normal Stabilizers Ancilla:
        initial_circuit.append("H", self.geometry.stab_idx)
        initial_circuit.append("TICK")

        initial_circuit.append("M", self.geometry.stab_idx)
        initial_circuit.append("TICK")

        return initial_circuit

    def _get_detectors(self):
        # Initlize Detector Circuit
        det_circuit = stim.Circuit()

        ##########################################################################
        # Implementing Detectors for Ancilla (inital Basis is deterministic)
        ##########################################################################

        # Determining Position in the measurement Run of only the Ancilla
        pos_to_index_ver: list = []
        pos_to_index_hor: list = []

        for pos, index in enumerate(self.geometry.stab_idx):
            if index in self.geometry.stab_ver_idx:
                pos_to_index_ver.append([pos, index])
            elif index in self.geometry.stab_hor_idx:
                pos_to_index_hor.append([pos, index])

        if self.geometry.state_init in {"XZZX-VER"}:
            # Adding the needed Detectors
            for index_pos in pos_to_index_ver:
                current_tar = index_pos[0] - len(self.geometry.stab_idx)
                q_index = index_pos[1]
                det_circuit.append(
                    "DETECTOR",
                    [stim.target_rec(current_tar)],
                    (self.geometry.i2q[q_index].real, self.geometry.i2q[q_index].imag, 0),
                )

        elif self.geometry.state_init in {"XZZX-HOR"}:
            # Adding the needed Detectors
            for index_pos in pos_to_index_hor:
                current_tar = index_pos[0] - len(self.geometry.stab_idx)
                q_index = index_pos[1]
                det_circuit.append(
                    "DETECTOR",
                    [stim.target_rec(current_tar)],
                    (self.geometry.i2q[q_index].real, self.geometry.i2q[q_index].imag, 0),
                )

        else:
            raise ValueError("Not a valid Basis for initlization in the Control Lattice")

        return det_circuit
