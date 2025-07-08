from typing import Dict

Label = str
Coord = complex

__all__ = ["build_lattice"]

# -----------------------
# Public Function
# -----------------------

def build_lattice(distance: int, state_init : str, *, starting_stabilizer_A: bool = True) -> Dict[Coord, Label]:
    """
    Adding the geometry of the XZZX-Code as a helper function:
        ->Returns the {coords: label} dict used for building the qubit coords

    Parameters:
    -----------
    distance : int
        Code distance (must be odd and ≥3)
    offset : complex, default 0+0j
        Upper left corner displacement -> Reusement for the ancilla/target/control blocks
    """

    """
    In theroy one needs to change the application of the boundary stab. 
    for even and odd distances but odd distances are the only practical code 
    length so we check for odd distance right at the beginning(makes no sene to 
    include them here for the boundary stab)
    """

    if distance <= 2 or distance % 2 == 0:
        raise ValueError("distance must be odd and ≥3")

    qubit_coords: Dict[Coord, Label] = {}
    start_with_x = starting_stabilizer_A

    for real in range(distance * 2):
        stab_counter = 0
        data_counter = 0

        for imag in range(distance * 2):
            # ---------------------- DATA qubits ---------------------------
            if real % 2 != 0 and imag % 2 != 0:
                coord = complex(real, imag)
                qubit_coords[coord] = "DATA"

                if start_with_x:
                    use_x = data_counter % 2 == 0
                else:
                    use_x = data_counter % 2 == 1

                if state_init in {"Z0", "Z1"}:
                    qubit_coords[coord] = "DATA_X" if use_x else "DATA_Z"
                    data_counter += 1

                elif state_init in {"X+", "X-"}:
                    qubit_coords[coord] = "DATA_Z" if use_x else "DATA_X"
                    data_counter += 1

                else:
                    ValueError("No valid inital state")

            # ----------------- Interior Stabilisers ----------------------
            elif real % 2 == 0 and imag % 2 == 0 and real != 0 and imag != 0:
                coord = complex(real, imag)

                if start_with_x:
                    use_x = stab_counter % 2 == 0
                else:
                    use_x = stab_counter % 2 == 1

                '''
                Stab A -> First Round Stabs are first one top left, then third top...
                Stab B -> First Round Stabs are second one top left, then fifth top...
                '''

                qubit_coords[coord] = "STAB-A" if use_x else "STAB-B"
                stab_counter += 1

        # flip phase after each even row (except first)
        if real % 2 == 0 and real != 0:
            start_with_x = not start_with_x

    return qubit_coords