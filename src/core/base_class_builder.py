from abc import ABC, abstractmethod

import stim

from src.core.data_models import NoiseParameters
from src.core.noise_models import BiasNoise, CircuitNoise


class BaseClassBuilder(ABC):
    @abstractmethod
    def __init__(self):
        self.noise = None
        pass

    @abstractmethod
    def build_circuit(self) -> stim.Circuit:
        """
        Abstract Method to build the circuit.

        Returns:
            stim.Circuit
        """
        pass

    @abstractmethod
    def get_logical_meas_rec(self, observable_index: int) -> list[int]:
        """
        Abstract Method to get the list of measurement record positions that need to be
        xored together to get the final logical measurement.

        Returns:
            list[int]: List of measurement record positions
        """
        pass

    def apply_noise(
        self,
        input_circuit: stim.Circuit,
        distance:int ,
        noise: NoiseParameters = None,
        ft_init: bool = False,
        ft_meas: bool = False,
    ) -> stim.Circuit:
        """
        Method to apply noise models to the circuit.
        """

        # Initialize Noise Model
        if noise is None:
            noise = NoiseParameters(
                after_c_depol_prob=0.0,
                before_round_depol=0.0,
                before_m_flip_prob=0.0,
                after_r_flip=0.0,
                after_c_pauli_channel_prob=0.0,
                noise_bias=None,
            )
        self.noise = noise

        # 1) Circuit Noise Model
        if self.noise.noise_bias is None and self.noise.after_c_pauli_channel_prob == 0.0:
            noise_dict = {
                "before_round_depol": self.noise.before_round_depol,
                "before_m_flip_prob": self.noise.before_m_flip_prob,
                "after_r_flip": self.noise.after_r_flip,
                "after_c_depol_prob": self.noise.after_c_depol_prob,
            }

            # Apply Noise Model
            circuit_noise_builder = CircuitNoise(circuit=input_circuit,
                                                 noise=noise_dict,
                                                 distance = distance,
                                                 ft_init = ft_init,
                                                 ft_measurements = ft_meas)
            input_circuit = circuit_noise_builder.apply()

        # 2) Biased Noise Model
        if self.noise.after_c_pauli_channel_prob > 0.0 or self.noise.noise_bias is not None:
            
            noise_dict = {
                "after_c_custom_noise": self.noise.after_c_pauli_channel_prob,
                "bias": self.noise.noise_bias,
                "before_m_flip_prob": self.noise.before_m_flip_prob,
                "after_r_flip": self.noise.after_r_flip,
            }

            # Apply Noise Model
            bias_noise_builder = BiasNoise(circuit=input_circuit, 
                                            noise=noise_dict,
                                            distance = distance,
                                            ft_init = ft_init,
                                            ft_measurements = ft_meas)
            input_circuit = bias_noise_builder.apply()

        return input_circuit
