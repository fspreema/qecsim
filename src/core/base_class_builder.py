from abc import ABC, abstractmethod


class BaseClassBuilder(ABC):
    @abstractmethod
    def build_circuit(self) -> None:
        """
        Abstract Method to build the circuit.

        Returns:
            None
        """
        pass

    @abstractmethod
    def _get_detectors(self) -> None:
        """
        Abstract Method to get the detectors.

        Returns:
            None
        """
        pass
