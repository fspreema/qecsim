import stim

from qecsim.core.data_models import (
    CircuitResult,
    Context,
    NoiseModel,
    Patch,
)

Coord = complex

__all__ = ["y_rev_switch_circ"]


def y_rev_switch_circ(
    *, lct: Context, patches: dict[str, Patch], offset: complex = 0 + 0j, noise: NoiseModel,
) -> CircuitResult:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # -Loading in Patches
    patch = patches["patch"]

    # -Retrieving Global Infomration
    q2i = lct.q2i
    stab_to_data_switch = lct.stab_to_data_modified
    stab_to_data_switch_xcy = lct.stab_to_data_modified2

    # Getting infromation about the additional newly formed boundary operators:
    r_h_stabs = patch.right_h
    u_h_stabs = patch.upper_h

    # -Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab
    switch_stab_apply_h = patch.stab_switch_apply_h

    # Finding Upper right qubit index -> need to look in 2-CX
    y_coords = 1 + 1j + offset
    y_index = q2i[y_coords]

    #################################
    # Creating final MPP measurements
    #################################

    final_measurement = stim.Circuit()
    final_measurement.append("TICK")

    #####################################
    # Define First Stab X MPP measurement
    #####################################

    reversed_switch_circ = stim.Circuit()
    reversed_switch_circ.append("TICK")

    """
    The following procedure will be done
    -> Along the mirrored diagonal one half h gate one half not
    -> Boundary on the site of the h gate gets expanded
        -> Every postion is now a boundary! (No 2 coords distance between them)
    -> Along the Y digaonal we do SQRT_X_DAG on all ancilla qubits
    """

    # Getting infromation about the additional newly formed boundary operators:
    r_h_stabs = patch.right_h
    u_h_stabs = patch.upper_h

    # Defining new diagonal and the corresponding stabilizers after the switch:

    index_h = []
    index_x_deg = []
    index_nh = []

    for cords, qtype in patch.coords.items():
        if cords != y_coords:
            # Diagonal Cut
            if cords.real - offset.real > cords.imag - offset.imag:
                index_h.append(q2i[cords])

            # Filtering out the X_DAG -> Not on Data
            elif cords.real - offset.real == cords.imag - offset.imag:
                if qtype != "DATA":
                    index_x_deg.append(q2i[cords])

            else:
                index_nh.append(q2i[cords])

    h_gates_rep: list = []

    for qubit, q_type in patch.coords.items():
        if q_type in {"X-STAB", "X-STAB-BOUND-B"}:
            if q2i[qubit] in index_h:
                continue
            else:
                h_gates_rep.append(q2i[qubit])

        elif q_type == "Z-STAB":
            if q2i[qubit] in index_h:
                h_gates_rep.append(q2i[qubit])
            else:
                continue

        if q_type in {"Z-STAB-BOUND-U-H", "Z-STAB-BOUND-U"}:
            h_gates_rep.append(q2i[qubit])

    reversed_switch_circ.append("R", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)
    reversed_switch_circ.append("TICK")
    reversed_switch_circ.append("H", switch_stab_apply_h)
    reversed_switch_circ.append("TICK")

    #########################################
    # Adding the XCY gates after the H switch
    #########################################

    # 2) CX Operations

    for coord_pairs, order in stab_to_data_switch.items():
        # Parallel Implementation of CX
        if order == "5TICK":
            if len(coord_pairs) == 2:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                reversed_switch_circ.append("CX", index_pairs)
            else:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                reversed_switch_circ.append("CX", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
        # Parallel Implementation of CX
        if order == "4TICK":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("CX", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
        # Parallel Implementation of CX
        if order == "3.5TICK":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("XCY", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
        # Parallel Implementation of CX
        if order == "3TICK":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("CX", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
        # Parallel Implementation of CX
        if order == "2TICK":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("CX", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch_xcy.items():
        # Parallel Implementation of CX
        if order == "1TICK":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("XCY", index_pairs)

    reversed_switch_circ.append("TICK")

    # -------Continue-Circuit------------

    # Adding stabs h
    reversed_switch_circ.append("H", index_h)
    reversed_switch_circ.append("SQRT_X_DAG", index_x_deg)
    reversed_switch_circ.append("TICK")

    # Adding half diagonal H
    reversed_switch_circ.append("H", x_stab_index + r_h_stabs)
    reversed_switch_circ.append("TICK")

    # Adding locial readout error -> Y logical meassured thorugh the stabilizers

    if noise.before_m_flip_prob > 0:
        reversed_switch_circ.append(
            "X_ERROR", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs + [y_index], noise.before_m_flip_prob,
        )

    # Adding Resets
    reversed_switch_circ.append("MZ", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)
    reversed_switch_circ.append("MY", y_index)

    # -> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    reversed_switch_circ.append("SHIFT_COORDS", arg=(0, 0, 1))

    return CircuitResult(circuit=reversed_switch_circ)
