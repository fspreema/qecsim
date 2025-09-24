import stim
from .dataclasses import Config, Patch, Context, NoiseModel, CircuitResult

Coord = complex

__all__ = ["y_memory_circ"]

def y_memory_circ(*, 
                  lct: Context, 
                  patches: dict[str, Patch], 
                  cfg: Config, 
                  noise: NoiseModel) -> CircuitResult:
    
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    patch = patches["patch"]

    #-Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    stab_to_data_switch = lct.stab_to_data_modified
    stab_to_data_switch_xcy = lct.stab_to_data_modified2

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab
    switch_stab_apply_h = patch.stab_switch_apply_h

    # Finding Upper right qubit index -> need to look in 2-CX
    y_coords = 1 + 1j
    y_index = q2i[y_coords]

    # Getting infromation about the additional newly formed boundary operators:
    r_h_stabs = patch.right_h
    u_h_stabs = patch.upper_h

    # Define Index which h where applied to
    index_h = []
    index_x_deg = []
    index_nh = []

    for cords, qtype in patch.coords.items():

        if cords != y_coords:

            # Diagonal Cut
            if cords.real > cords.imag:
                index_h.append(q2i[cords])

            # Filtering out the X_DAG -> Not on Data
            elif cords.real == cords.imag:
                if qtype != "DATA":
                    index_x_deg.append(q2i[cords])

            else:
                index_nh.append(q2i[cords])

    ###########################
    # Define Repetition Circuit
    ###########################

    memory_y = stim.Circuit()

    #########################################
    # Adding the XCY gates after the H switch
    #########################################

    """
    As now the h gates where applied in the switch round we need to relable HALF of the stabilizers!
    """

    h_gates_rep : list = []

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


    #3) Appl H gate after reset
    memory_y.append("H", h_gates_rep) 
    memory_y.append("TICK")

    #2) CX Operations

    for coord_pairs, order in stab_to_data_switch_xcy.items():
   
        #Parallel Implementation of CX
        if order == "1TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            memory_y.append("XCY", index_pairs)

    memory_y.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "2TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            memory_y.append("CX", index_pairs)

    memory_y.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            memory_y.append("CX", index_pairs)

    memory_y.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3.5TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            memory_y.append("XCY", index_pairs)

    memory_y.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "4TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            memory_y.append("CX", index_pairs)
   
    memory_y.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "5TICK":

            if len(coord_pairs) == 2:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                memory_y.append("CX", index_pairs)
            else:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                memory_y.append("CX", index_pairs)

    #-------Continue-Circuit------------
    
    memory_y.append("TICK")

    #3) Basis/ Measurement
    memory_y.append("H", switch_stab_apply_h)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        memory_y.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    memory_y.append("TICK")

    memory_y.append("MR", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)

    memory_y.append("TICK")

    # Adding missing detectors
    for index, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
        prev_tar = -2 * memory_y.num_measurements + index
        current_tar = -1 * memory_y.num_measurements + index
        # memory_y.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
        #                      (i2q[q_index].real, i2q[q_index].imag, 0))


    for index, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
        role = patch.coords[i2q[q_index]]
        coord = i2q[q_index]

        if role == "X-STAB":
            # skip stabs on the Y-cut
            if q_index in index_nh:
                prev_tar = -2 * memory_y.num_measurements + index
                current_tar = -1 * memory_y.num_measurements + index
                memory_y.append(
                    "DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0)
                )

        elif role in {"X-STAB-BOUND-R"}:

            if q2i[(coord.real - 1 + 1j * coord.imag + 1j)] in index_h:
                prev_tar = -2 * memory_y.num_measurements + index
                current_tar = -1 * memory_y.num_measurements + index
                memory_y.append(
                    "DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0)
                )

        elif role in {"X-STAB-BOUND-B", "Z-STAB-BOUND-L", "Z-STAB-BOUND-U"}:
            prev_tar = -2 * memory_y.num_measurements + index
            current_tar = -1 * memory_y.num_measurements + index
            memory_y.append(
                "DETECTOR",
                [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                (coord.real, coord.imag, 0)
            )

        elif role == "Z-STAB":
            # skip stabs on the Y-cut
            if q_index in index_nh:
                prev_tar = -2 * memory_y.num_measurements + index
                current_tar = -1 * memory_y.num_measurements + index
                memory_y.append(
                    "DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0)
                )

    # memory_y += memory_y.missing_detectors()

    return CircuitResult(circuit=memory_y)
