from itertools import chain
import stim
from .data_models import Config, Patch, Context, NoiseModel, CircuitResult

Coord = complex

__all__ = ["y_rev_switch_circ"]

def y_rev_switch_circ(*, 
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
    rounds = cfg.rounds
    distance = cfg.distance
    stab_to_data = lct.stab_to_data
    stab_to_data_switch = lct.stab_to_data_modified
    stab_to_data_switch_xcy = lct.stab_to_data_modified2

    #-Retrieving Data Coords
    data = patch.data

    # Getting infromation about the additional newly formed boundary operators:
    r_h_stabs = patch.right_h
    u_h_stabs = patch.upper_h

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab
    x_stab_index_memory = patch.x_stab_memory
    z_stab_index_memory = patch.z_stab_memory
    switch_stab_apply_h = patch.stab_switch_apply_h

    # Finding Upper right qubit index -> need to look in 2-CX
    y_coords = 1 + 1j
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
            if cords.real > cords.imag:
                index_h.append(q2i[cords])

            # Filtering out the X_DAG -> Not on Data
            elif cords.real == cords.imag:
                if qtype != "DATA":
                    index_x_deg.append(q2i[cords])

            else:
                index_nh.append(q2i[cords])

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

    reversed_switch_circ.append("R", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)
    reversed_switch_circ.append("TICK")
    reversed_switch_circ.append("H", switch_stab_apply_h)
    reversed_switch_circ.append("TICK")

    #########################################
    # Adding the XCY gates after the H switch
    #########################################

    #2) CX Operations

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
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
   
        #Parallel Implementation of CX
        if order == "4TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("CX", index_pairs)
   
    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3.5TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("XCY", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("CX", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "2TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("CX", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch_xcy.items():
   
        #Parallel Implementation of CX
        if order == "1TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("XCY", index_pairs)

    reversed_switch_circ.append("TICK")

    #-------Continue-Circuit------------

    # Adding stabs h
    reversed_switch_circ.append("H", index_h)
    reversed_switch_circ.append("SQRT_X_DAG", index_x_deg)
    reversed_switch_circ.append("TICK")

    #Adding half diagonal H
    reversed_switch_circ.append("H", x_stab_index + r_h_stabs)
    reversed_switch_circ.append("TICK")

    #Adding Resets
    reversed_switch_circ.append("MZ", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)
    reversed_switch_circ.append("MY", y_index)

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    reversed_switch_circ.append("SHIFT_COORDS", arg = (0,0,1))

    ###########################################################
    # Adding Detectors (Newly Gen Boundary and Old Stabilizers)
    ###########################################################

    #1) Deterministic Detectors which can be build up by the old stabilizers
    """
    Reuse the same procedure as for the switch but reversed
    """

    previous_round =  len(x_stab_index + z_stab_index)
    current_round = reversed_switch_circ.num_measurements

    for index, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs + [y_index]):

        # Calc important info
        coord = i2q[q_index]
        role = patch.coords[coord]
        
        if role == "X-STAB":

            # skip stabs on the Y-cut
            if q_index in index_nh:

                # find this stabilizers position in the previous block
                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                reversed_switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                
                
            # Everything but the off diagonal
            elif i2q[q_index] not in [2 + i + 1j * i for i in range(2,distance * 2 + 1, 2)]:

                left_neighbour = coord.real - 2 + coord.imag*1j
                left_index = q2i[left_neighbour]

                # find this stabilizers position in the previous block (Z stab ancilla nex to it)
                for curr_idx, q_index in enumerate(patch.x_stab_memory + patch.z_stab_memory):
                    if q_index == left_index:
                        prev_tar = - current_round - previous_round + curr_idx
                
                current_tar = - current_round + index
                reversed_switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                
            # Off Diagonal
            else:
                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                #reversed_switch_circ.append("DETECTOR",
                #    [stim.target_rec(current_tar)],
                #    (coord.real, coord.imag, 0))

        elif role == "Z-STAB":

            diagonal_indices = [2 + i + 2j + 1j * i for i in range(0,distance * 2 - 2, 2)]

            #skip stabs on the Y-cut
            if q_index in index_nh:
                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                reversed_switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                
                
            # Skip Diagonal
            elif i2q[q_index] not in diagonal_indices:

                bottom_neighbour = coord.real + coord.imag*1j + 2j
                bottom_index = q2i[bottom_neighbour]

                # find this stabilizers position in the previous block (Z stab ancilla nex to it)
                for curr_idx, q_index in enumerate(patch.x_stab_memory + patch.z_stab_memory):
                    if q_index == bottom_index:
                        prev_tar = - current_round - previous_round + curr_idx
                
                current_tar = - current_round + index
                reversed_switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))

            # Diagonal 
            else:
                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                #reversed_switch_circ.append("DETECTOR",
                #    [stim.target_rec(current_tar)],
                #    (coord.real, coord.imag, 0))

        elif role in {"X-STAB-BOUND-B", "Z-STAB-BOUND-L"}:
            prev_tar = - current_round - previous_round + index
            current_tar = - current_round + index
            reversed_switch_circ.append("DETECTOR",
                [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                (coord.real, coord.imag, 0))
            
        elif role in {"X-STAB-BOUND-R-H"}:

            # find this stabilizers position
            for curr_idx, curr_q_index in enumerate(patch.x_stab_memory + patch.z_stab_memory):
                if q_index == curr_q_index:
                    prev_tar = - current_round - previous_round + curr_idx
                    
            current_tar = - current_round + index
            reversed_switch_circ.append("DETECTOR",
                [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                (coord.real, coord.imag, 0))


        elif role in {"Z-STAB-BOUND-U-H"}:

            # find this stabilizers position
            for curr_idx, curr_q_index in enumerate(patch.x_stab_memory + patch.z_stab_memory):
                if q_index == curr_q_index:
                    prev_tar = - current_round - previous_round + curr_idx
                    
            current_tar = - current_round + index
            reversed_switch_circ.append("DETECTOR",
                [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                (coord.real, coord.imag, 0))


    return CircuitResult(circuit=reversed_switch_circ)