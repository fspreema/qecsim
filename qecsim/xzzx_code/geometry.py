from typing import Dict

Label = str
Coord = complex

__all__ = ["build_lattice"]

# -----------------------
# Public Function
# -----------------------

def build_lattice(distance: int, state_init : str) -> Dict[Coord, Label]:
    """
    Adding the geometry of the XZZX-Code as a helper function:
        ->Returns the {coords: label} dict used for building the qubit coords

    *Parameters:*
    ------------
    distance : int
        Code distance (must be odd and ≥3)
    state_init : str
        What basis (Vertical/ Horizontal) is the code initialized in?

    *Comments:*
    ----------
    *Basis:*
    As there are no trivial loops build up by only X/Z Paulis whe specify 
    the given basis by the orientation of the logical Operator i.e. vertical/ horizontal

    *Distance:*
    In theroy one needs to change the application of the boundary stab. 
    for even and odd distances but odd distances are the only practical code 
    length so we check for odd distance right at the beginning(makes no sene to 
    include them here for the boundary stab)
    """

    if distance <= 2 or distance % 2 == 0:
        raise ValueError("distance must be odd and ≥3")

    qubit_coords: Dict[Coord, Label] = {}

    start_with_x = True

    for real in range(distance * 2):
        stab_counter = 0
        data_counter = 0

        for imag in range(distance * 2):

            #######################
            # Labeling DATA QUBITS
            #######################

            if real % 2 != 0 and imag % 2 != 0:
                coord = complex(real, imag)
                qubit_coords[coord] = "DATA"

                if start_with_x:
                    use_x = data_counter % 2 == 0
                else:
                    use_x = data_counter % 2 == 1

                """
                Labeling of the Data qubits
                -> Depending on the initlization in the Horizontal/ Vertical Basis one has to initlize the data qubits differently
                -> Labeling ensures which Reset to apply on which qubit
                """

                if state_init in {"Ver"}:
                    qubit_coords[coord] = "DATA_X" if use_x else "DATA_Z"
                    data_counter += 1

                elif state_init in {"Hor"}:
                    qubit_coords[coord] = "DATA_Z" if use_x else "DATA_X"
                    data_counter += 1

                else:
                    ValueError("No valid inital state")

            ######################################
            # Labeling Stabilizer (Aniclla) Qubits
            ######################################

            elif real % 2 == 0 and imag % 2 == 0 and real != 0 and imag != 0:
                coord = complex(real, imag)

                if start_with_x:
                    use_x = stab_counter % 2 == 0
                else:
                    use_x = stab_counter % 2 == 1

                """
                Stab-Ver -> Deterministic Stabilizers (Detectors) with respect to the Vertical Operator/ Basis
                Stab-Hor -> Deterministic Stabilizers (Detectors) with respect to the Horizontal Operator/ Basis
                """

                qubit_coords[coord] = "STAB-Ver" if use_x else "STAB-Hor"
                stab_counter += 1

        # flip phase after each even row (except first)
        if real % 2 == 0 and real != 0:
            start_with_x = not start_with_x

    return qubit_coords