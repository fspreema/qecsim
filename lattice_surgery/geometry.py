from typing import Dict

Label = str
Coord = complex

__all__ = ["build_lattice"]

# -----------------------
# Public Function
# -----------------------

def build_lattice(distance: int, *, offset: Coord = 0 + 0j, starting_stabilizer_x: bool = True) -> Dict[Coord, Label]:
    """
    Adding the geometry of the Lattice Surgery as a helper function:
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

    ox, oy = int(offset.real), int(offset.imag)

    qubit_coords: Dict[Coord, Label] = {}
    start_with_x = starting_stabilizer_x

    for real in range(distance * 2):
        stab_counter = 0

        for imag in range(distance * 2):
            # ---------------------- DATA qubits ---------------------------
            if real % 2 != 0 and imag % 2 != 0:
                coord = complex(real + ox, imag + oy)
                qubit_coords[coord] = "DATA"

            # ----------------- Interior Stabilisers ----------------------
            elif real % 2 == 0 and imag % 2 == 0 and real != 0 and imag != 0:
                coord = complex(real + ox, imag + oy)

                if start_with_x:
                    use_x = stab_counter % 2 == 0
                else:
                    use_x = stab_counter % 2 == 1

                qubit_coords[coord] = "X-STAB" if use_x else "Z-STAB"
                stab_counter += 1

        # flip phase after each even row (except first)
        if real % 2 == 0 and real != 0:
            start_with_x = not start_with_x

    return qubit_coords

# Demo
if __name__ == "__main__":
    patch = build_lattice(3)
    print("Coords:", len(patch))
    # print first 5 entries to show format
    for k, v in patch.items():
        print(k, v)