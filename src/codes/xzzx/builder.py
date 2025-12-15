import stim

from src.codes.xzzx.circuits.final_measure import FinalMeasureCircuit
from src.codes.xzzx.circuits.initial import InitialCircuit
from src.codes.xzzx.circuits.repetition import RepetitionCircuit
from src.codes.xzzx.circuits.reset import ResetCircuit
from src.codes.xzzx.get_stab_pairings import XZZXPairings
from src.codes.xzzx.xzzx_geom import XZZXGeometry
from src.core.data_models import XZZXNoise
from src.core.noise_models import BiasNoise, CircuitNoise

Coord = complex
Label = str

__all__ = ["XZZXBuilder"]

# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------


class XZZXBuilder:
    def __init__(self, distance: int, state_init: str, noise: XZZXNoise = None):
        self.geometry = XZZXGeometry(
            distance=distance,
            state_init=state_init,
        )
        self.pairings = XZZXPairings(patch=self.geometry.coords)

        # Initialize Noise Model
        if noise is None:
            noise = XZZXNoise(
                after_c_depol_prob=0.0,
                before_round_depol=0.0,
                before_m_flip_prob=0.0,
                after_r_flip=0.0,
                after_c_pauli_channel_prob=0.0,
                noise_bias=None,
            )
        self.noise = noise

    def build_circuit(self) -> stim.Circuit:
        # Initialize Empty Circuit
        self.full_circuit = stim.Circuit()

        # 1) Adding State Reset Circuit
        reset_circ = ResetCircuit(
            self.geometry,
            self.pairings,
            self.noise,
        )
        self.full_circuit += reset_circ.build_reset_circuit()

        # 2) Adding Initial Circuit
        init_circ = InitialCircuit(
            self.geometry,
            self.pairings,
        )
        self.full_circuit += init_circ.build_initial_circuit()

        # 3) Adding Repetition Circuit
        rep_circ = RepetitionCircuit(
            self.geometry,
            self.pairings,
        )
        self.full_circuit += rep_circ.build_repetition_circuit()

        # 4) Adding Final Measurement Circuit
        final_circ = FinalMeasureCircuit(
            self.geometry,
        )
        self.full_circuit += final_circ.build_final_measurement_circuit()

        # 5) Applying Noise Model
        return self._apply_noise(input_circuit=self.full_circuit)

    def _apply_noise(self, input_circuit: stim.Circuit) -> stim.Circuit:
        # 1) Circuit Noise Model
        if (
            self.noise.before_m_flip_prob > 0.0
            or self.noise.after_r_flip > 0.0
            or self.noise.after_c_depol_prob > 0.0
            or self.noise.before_round_depol > 0.0
        ):
            noise_dict = {
                "before_round_depol": self.noise.before_round_depol,
                "before_m_flip_prob": self.noise.before_m_flip_prob,
                "after_r_flip": self.noise.after_r_flip,
                "after_c_depol_prob": self.noise.after_c_depol_prob,
            }

            # Apply Noise Model
            circuit_noise_builder = CircuitNoise(circuit=input_circuit, noise=noise_dict)
            input_circuit = circuit_noise_builder.apply()

        # 2) Biased Noise Model
        if self.noise.after_c_pauli_channel_prob not in (
            0.0,
            None,
        ) or self.noise.noise_bias not in (None, []):
            noise_dict = {
                "after_c_custom_noise": self.noise.after_c_pauli_channel_prob,
                "bias": self.noise.noise_bias,
            }

            # Apply Noise Model
            bias_noise_builder = BiasNoise(circuit=input_circuit, noise=noise_dict)
            input_circuit = bias_noise_builder.apply()

        return input_circuit
