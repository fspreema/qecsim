import stim


class GrowCircuit:
    """
    Class to gorw a given lattice by a multiple of distance 2
    """

    def __init__(self, circuit: stim.Circuit, distance: int):
        self.circuit = circuit
        self.distance = distance

    def grow(self) -> stim.Circuit:
        pass
