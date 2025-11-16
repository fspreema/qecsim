import stim

"""
Generates the Coordinates needed for a Color Code
-> Labels the Data Qubits
-> Labels the Stabilizers

Also allows offset in real and imaginary axis in order to allow for multiple patches
"""

__all__ = ["ColorCode"]

# For now the Coordiantes are Hardcoded as we will fixaed on the distance 3 Color Code
COLOR_0: complex = complex(5, 15)
COLOR_1: complex = complex(6, 12)
COLOR_2: complex = complex(3, 15)
COLOR_3: complex = complex(3, 13)
COLOR_4: complex = complex(2, 10)
COLOR_5: complex = complex(2, 12)
COLOR_6: complex = complex(1, 15)

# Also Define the X and Z Ancillas

X_COLOR_0: complex = complex(4, 14)
Z_COLOR_0: complex = complex(6, 14)
X_COLOR_1: complex = complex(5, 11)
Z_COLOR_1: complex = complex(7, 11)
X_COLOR_2: complex = complex(2, 14)
Z_COLOR_2: complex = complex(4, 16)


# Create fixed Dict for Color Code
COLOR_CODE_DICT = {
    COLOR_0: "DATA",
    COLOR_1: "DATA",
    COLOR_2: "DATA",
    COLOR_3: "DATA",
    COLOR_4: "DATA",
    COLOR_5: "DATA",
    COLOR_6: "DATA",
    X_COLOR_0: "X",
    Z_COLOR_0: "Z",
    X_COLOR_1: "X",
    Z_COLOR_1: "Z",
    X_COLOR_2: "X",
    Z_COLOR_2: "Z",
}

# Define Mapping Needed for Color Code


class ColorCode:
    def __init__(self, distance: int = 3, cords_dict: dict = COLOR_CODE_DICT):
        self.distance = distance

        # Get Index Mappings
        self.q2i: dict[complex, int] = {
            q: i
            for i, q in enumerate(
                sorted(cords_dict.keys(), key=lambda v: (v.real, v.imag)),
            )
        }

    def initilize(self, state: str) -> stim.Circuit:
        if state != "+":
            raise NotImplementedError("Only + state initialization is implemented.")

        init_circ = stim.Circuit()

        # First add Qubit Coords
        for q, i in self.q2i.items():
            init_circ.append("QUBIT_COORDS", [i], [q.real, q.imag])

        return init_circ

    def injection(self) -> stim.Circuit:
        inj_circ = stim.Circuit()

        return inj_circ

    def verify(self) -> stim.Circuit:
        ver_circ = stim.Circuit()

        return ver_circ
