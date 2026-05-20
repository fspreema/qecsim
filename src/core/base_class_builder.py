from abc import ABC, abstractmethod

import stim

from src.core.base_geometry import BaseGeometry
from src.core.data_models import NoiseParameters
from src.core.get_measurement_recs import get_measurement_recs
from src.core.noise_models import CircuitNoise, PenomenologicalNoise


class BaseClassBuilder(ABC):
    @abstractmethod
    def __init__(self):
        self.noise = None
        self.return_circuit = None
        pass

    @abstractmethod
    def build_circuit(self) -> stim.Circuit:
        """
        Abstract Method to build the circuit.

        Returns:
            stim.Circuit
        """
        pass

    def get_logical_meas_rec(self, observable_index: int) -> list[int]:
        """
        Returns the measurement records that build up the logical operator
        """

        if self.return_circuit is None:
            raise ValueError("Circuit has not been built yet. Please build the circuit first.")

        measurement_records = get_measurement_recs(
            circuit=self.return_circuit,
            observable_index=observable_index,
        )

        return measurement_records

    def apply_noise(
        self,
        input_circuit: stim.Circuit,
        geometry: BaseGeometry,
        num_tick_first_noise: int,
        num_tick_last_noise: int,
        noise: NoiseParameters = None,
        ft_init: bool = False,
        ft_meas: bool = False,
    ) -> stim.Circuit:
        """
        Method to apply noise models to the circuit.
        """

        # Initlize Noise Model and Check what noise to apply
        if noise is None:
            noise = NoiseParameters()
        self.noise = noise

        use_circuit_noise = self.noise.circuit_noise_prob > 0.0
        use_pheno_noise   = self.noise.phenomenological_noise_prob > 0.0
        use_bias          = self.noise.noise_bias is not None

        shared_kwargs = dict(
            circuit=input_circuit,
            geometry=geometry,
            ft_init=ft_init,
            ft_measurements=ft_meas,
            num_tick_first_noise=num_tick_first_noise,
            num_tick_last_noise=num_tick_last_noise,
        )

        if use_bias and use_circuit_noise:
            # Init Circuit Noise but with Biased
            noise_dict = {
                "CircuitNoiseProbability": self.noise.circuit_noise_prob,
                "bias": self.noise.noise_bias,
            }
            input_circuit = CircuitNoise(noise=noise_dict, **shared_kwargs).apply()

        elif use_bias and use_pheno_noise:
            # Init Pheno Noise but with Biased
            noise_dict = {
                "PhemoNoiseProbability": self.noise.phenomenological_noise_prob,
                "bias": self.noise.noise_bias,
            }
            input_circuit = PenomenologicalNoise(noise=noise_dict, **shared_kwargs).apply()

        elif use_circuit_noise:
            noise_dict = {"CircuitNoiseProbability": self.noise.circuit_noise_prob}
            input_circuit = CircuitNoise(noise=noise_dict, **shared_kwargs).apply()

        elif use_pheno_noise:
            noise_dict = {
                "PhemoNoiseProbability": self.noise.phenomenological_noise_prob,
            }
            input_circuit = PenomenologicalNoise(noise=noise_dict, **shared_kwargs).apply()

        return input_circuit
