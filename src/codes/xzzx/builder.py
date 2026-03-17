import stim

from src.codes.xzzx.circuits.final_measure import FinalMeasureCircuit
from src.codes.xzzx.circuits.initial import InitialCircuit
from src.codes.xzzx.circuits.repetition import RepetitionCircuit
from src.codes.xzzx.circuits.reset import ResetCircuit
from src.codes.xzzx.get_stab_pairings import XZZXPairings
from src.codes.xzzx.xzzx_geom import XZZXGeometry
from src.core.base_class_builder import BaseClassBuilder
from src.core.data_models import NoiseParameters

Coord = complex
Label = str

__all__ = ["XZZXBuilder"]

# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------


class XZZXBuilder(BaseClassBuilder):
    def __init__(self, distance: int, state_init: str, noise: NoiseParameters = None):
        self.geometry = XZZXGeometry(
            distance=distance,
            state_init=state_init,
        )
        self.pairings = XZZXPairings(patch=self.geometry.coords)
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
        return self.apply_noise(input_circuit=self.full_circuit, 
                                noise=self.noise, 
                                geometry=self.geometry, 
                                distance=self.geometry.distance)
