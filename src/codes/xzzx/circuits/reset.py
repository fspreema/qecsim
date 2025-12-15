import numpy as np
import stim

from src.codes.xzzx.get_stab_pairings import XZZXPairings
from src.codes.xzzx.xzzx_geom import XZZXGeometry
from src.core.data_models import XZZXNoise

__all__ = ["ResetCircuit"]


class ResetCircuit:
    def __init__(
        self,
        geometry: XZZXGeometry,
        stab_pairings: XZZXPairings,
        noise: XZZXNoise,
    ):
        self.geometry = geometry
        self.stab_pairings = stab_pairings
        self.noise = noise

    def build_reset_circuit(self) -> stim.Circuit:
        # Initialize Empty Circuit
        self.reset_circ = stim.Circuit()

        # Apply Reset Operations
        self.reset_circ = self._apply_reset()

        # Apply Pre-Round Noise
        self.reset_circ = self._apply_pre_round_noise()

        return self.reset_circ

    def _apply_reset(self):
        ########################
        # Define Initial Circuit
        ########################

        """
        Looking at every state preperation seperatly seems to be inefficient
        ->  If not all Operators only once used one gets an incorrect formatting in the
            timeslice view because of the Operations being in different timeslices in each TICK!
        """

        ##################
        # Appending Coords
        ##################

        for q, i in self.geometry.q2i.items():
            self.reset_circ.append("QUBIT_COORDS", [i], [q.real, q.imag])

        ########################################################################
        # Inilizing Ancilla in Plus (Reset) and Control/ Target in desired State
        ########################################################################

        init_patterns = {
            ("XZZX-VER"): [
                ("RZ", self.geometry.data_z_idx),
                ("RX", self.geometry.data_x_idx),
            ],
            ("XZZX-HOR"): [
                ("RX", self.geometry.data_x_idx),
                ("RZ", self.geometry.data_z_idx),
            ],
        }

        # Apply the initialization pattern
        if self.geometry.state_init not in init_patterns:
            raise ValueError(f"Invalid basis combination: {self.geometry.state_init}")

        for gate, qubits in init_patterns[self.geometry.state_init]:
            self.reset_circ.append(gate, qubits)

        return self.reset_circ

    def _apply_pre_round_noise(self):
        # -------Adding Before Round Data Depol.------------
        if self.noise.before_round_depol > 0:
            self.reset_circ.append(
                "DEPOLARIZE1",
                self.geometry.data_idx,
                self.noise.before_round_depol,
            )
        # --------------------------------------------------

        # -------Adding Before Round Data Depol.------------
        if np.any(self.noise.before_round_p_xyz):
            self.reset_circ.append(
                "PAULI_CHANNEL_1",
                self.geometry.data_idx,
                self.noise.before_round_p_xyz,
            )
        # --------------------------------------------------

        return self.reset_circ
