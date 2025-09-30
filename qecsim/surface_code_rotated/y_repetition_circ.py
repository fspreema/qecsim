import stim
from .data_models import Config, Patch, Context, NoiseModel, CircuitResult

Coord = complex

__all__ = ["y_repetition_circ"]

def y_repetition_circ(*, 
                      lct: Context, 
                      patches: dict[str, Patch], 
                      cfg: Config, 
                      noise: NoiseModel,
                      memory_round: bool = False,
                      ft_round:bool = False) -> CircuitResult:
    
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    patch = patches["patch"]

    #-Retrieving Global Information
    if not memory_round:
        stab_to_data = lct.stab_to_data
    else:
        stab_to_data = lct.stab_to_data_modified3

    q2i = lct.q2i
    i2q = lct.i2q
    rounds = cfg.rounds
    distance = cfg.distance

    # Getting infromation about the additional newly formed boundary operators:
    r_h_stabs = patch.right_h
    u_h_stabs = patch.upper_h

    # Finding Upper right qubit index -> need to look in 2-CX
    y_coords = 1 + 1j
    y_index = q2i[y_coords]

    #What half is an h applied?
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

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    if not memory_round:
        x_stab_index = patch.x_stab
        z_stab_index = patch.z_stab
    else:
        x_stab_index = patch.x_stab_memory
        z_stab_index = patch.z_stab_memory

    ###########################
    # Define Repetition Circuit
    ###########################

    #-----BUILDING-REPETITION-CIRC------

    round_circuit = stim.Circuit()

    #-------Adding-Before-Round-Depol.-Data------------

    if noise.before_round_depol > 0:
        round_circuit.append("DEPOLARIZE1", data, noise.before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    round_circuit.append("TICK")
    round_circuit.append("R", x_stab_index + z_stab_index)
    round_circuit.append("TICK")
    round_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        round_circuit.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    round_circuit.append("TICK")

    #2) CX Operations

    def _pairs_for(order: str) -> list[list[int]]:
        # Return the Pairs needed at the current order
        return [[q2i[cp[1]], q2i[cp[0]]] for cp, o in stab_to_data.items() if o == order]

    def _append_by_order(op: str, order: str, noise: float = 0.0) -> None:
        # Getting pair info
        for pair in _pairs_for(order):
            if not memory_round:
                # Checking for upper corner CX and leave it out!
                if y_index not in pair:
                    #Adding Pair on Operation
                    if op == "CX":
                        round_circuit.append(op, pair)
                    elif op == "DEPOLARIZE2":
                        round_circuit.append(op, pair, noise)
            else:
                #Adding Pair on Operation
                if op == "CX":
                    round_circuit.append(op, pair)
                elif op == "DEPOLARIZE2":
                    round_circuit.append(op, pair, noise)


    # Adding all the CX gates
    for order in ("1-CX", "2-CX", "3-CX", "4-CX"):
        _append_by_order("CX", order)
        if noise.after_c_depol_prob > 0:
            _append_by_order("DEPOLARIZE2", order, noise.after_c_depol_prob)
        round_circuit.append("TICK")

    #-------Continue-Circuit------------

    #3) Basis/ Measurement
    round_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        round_circuit.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    round_circuit.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.-------

    if noise.before_m_flip_prob > 0:
        round_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.before_m_flip_prob)

    #-------Continue-Circuit----------

    round_circuit.append("M", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if noise.after_r_flip > 0:
        round_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.after_r_flip)


    #-------Continue-Circuit------------

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    round_circuit.append("SHIFT_COORDS", arg= (0,0,1))

    #4) Detectors
    num_measurements_repeat = len(x_stab_index + z_stab_index)

    if not memory_round and not ft_round:

        ###########################################
        # Adding normal detectors in repetition run
        ###########################################

        for index, q_index in enumerate(x_stab_index + z_stab_index):
            prev_tar = -2 * num_measurements_repeat + index
            current_tar = -1 * num_measurements_repeat + index
            round_circuit.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                                (i2q[q_index].real, i2q[q_index].imag, 0))
            
        rep_circ = round_circuit * int((rounds - 2) / 2)
        
        return CircuitResult(circuit=rep_circ)
    
    elif memory_round and not ft_round:

        #####################################################
        # Adding Detector from switch to repeating transition
        #####################################################

        pre_det = stim.Circuit()

        previous_round = len(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)
        current_round = round_circuit.num_measurements

        for index, q_index in enumerate(x_stab_index + z_stab_index):

            # Calc important info
            coord = i2q[q_index]
            role = patch.coords[coord]
            
            if role == "X-STAB":

                # skip where x-stab & boundary meet
                if coord.imag != 2:

                    # find this stabilizers position in the previous block
                    prev_tar = - current_round - previous_round + index
                    current_tar = - current_round + index
                    pre_det.append("DETECTOR",
                        [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                        (coord.real, coord.imag, 0))
                    
                # combine current record with old z stab boundary h
                else:

                    #Finding postion of boundary:
                    right_neighbour = coord.real + 2 + coord.imag * 1j
                    right_index = q2i[right_neighbour]

                    #Determingin old Index
                    for curr_idx, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
                        if q_index == right_index:
                            prev_tar1 = - current_round - previous_round + curr_idx

                    # Determining current index
                    current_tar = - current_round + index

                    #pre_det.append("DETECTOR",
                    #    [stim.target_rec(current_tar), stim.target_rec(prev_tar1)],
                    #    (coord.real, coord.imag, 0))



            elif role == "Z-STAB":

                # skip where z-stab & boundary meet
                if coord.real != distance * 2 - 2:

                    # find this stabilizers position in the previous block
                    prev_tar = - current_round - previous_round + index
                    current_tar = - current_round + index
                    pre_det.append("DETECTOR",
                        [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                        (coord.real, coord.imag, 0))
                    
                # Meeting points becomes weight 3 -> include boundary h 
                else:

                    #Finding postion of boundary:
                    right_upper_neighbour = coord.real + 2 + coord.imag * 1j - 2j
                    right_upper_index = q2i[right_upper_neighbour]

                    for curr_idx, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
                        if q_index == right_upper_index:
                            prev_tar1 = - current_round - previous_round + curr_idx

                    prev_tar2 = - current_round - previous_round + index
                    current_tar = - current_round + index

                    #pre_det.append("DETECTOR",
                    #    [stim.target_rec(current_tar), stim.target_rec(prev_tar2)],
                    #    (coord.real, coord.imag, 0))
                    
            elif role in {"X-STAB-BOUND-B", "Z-STAB-BOUND-L"}:
                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                pre_det.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                
            elif role in {"Z-STAB-BOUND-U-H"}:
                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                #pre_det.append("DETECTOR",
                #    [stim.target_rec(current_tar)],
                #    (coord.real, coord.imag, 0))

        #######################################  
        # Adding detectors for repeating rounds
        #######################################

        pre_round = round_circuit
        det_round = stim.Circuit()
        det_round += round_circuit

        for index, q_index in enumerate(x_stab_index + z_stab_index):
            prev_tar = -2 * num_measurements_repeat + index
            current_tar = -1 * num_measurements_repeat + index
            det_round.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                                (i2q[q_index].real, i2q[q_index].imag, 0))

        full_run = pre_round
        full_run += pre_det 
        full_run += det_round * (rounds - 1)

        return CircuitResult(circuit=full_run)
    
    elif ft_round and not memory_round:

        """
        We need to add the parity measurement needed to determine the y measurement (MY is done inside the switch circ as the last measurement)
        """

        ###########################################  
        # Adding detectors for fault tolerant round
        ###########################################

        # Adding one circuit with non trivial detectors and then adding the trivial repeating sam index detectors
        pre_round = round_circuit
        det_round = stim.Circuit()
        det_round += round_circuit

        # Non trivial round
        diagonal_indices = [2 + i + 2j + 1j * i for i in range(0,distance * 2 - 2, 2)]
        previous_round = len(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs + [y_index])
        current_round = round_circuit.num_measurements

        for index, q_index in enumerate(x_stab_index + z_stab_index):

            # Calc important info
            coord = i2q[q_index]
            role = patch.coords[coord]
            
            if role == "X-STAB":

                # skip stabs on the Y-cut
                if q_index in index_nh:

                    # find this stabilizers position in the previous block
                    prev_tar = - current_round - previous_round + index
                    current_tar = - current_round + index
                    pre_round.append("DETECTOR",
                        [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                        (coord.real, coord.imag, 0))
                    

            elif role == "Z-STAB":

                #skip stabs on the Y-cut
                if q_index in index_nh:

                    prev_tar = - current_round - previous_round + index
                    current_tar = - current_round + index
                    pre_round.append("DETECTOR",
                        [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                        (coord.real, coord.imag, 0))
                    
                # Everything but the off diagonal
                elif i2q[q_index] not in diagonal_indices and coord.real == distance * 2 - 2:

                    top_neighbour = coord.real + coord.imag*1j - 2j
                    top_index = q2i[top_neighbour]

                    # find this stabilizers position in the previous block (Z stab ancilla nex to it)
                    for curr_idx, q_index in enumerate(x_stab_index + z_stab_index):
                        if q_index == top_index:
                            prev_tar = - current_round - previous_round + curr_idx
                    
                    current_tar = - current_round + index

                    #pre_round.append("DETECTOR",
                    #    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    #    (coord.real, coord.imag, 0))    
                    

        # Trivial Det round
        for index, q_index in enumerate(x_stab_index + z_stab_index):
            prev_tar = -2 * num_measurements_repeat + index
            current_tar = -1 * num_measurements_repeat + index
            det_round.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                                (i2q[q_index].real, i2q[q_index].imag, 0))

        full_run = pre_round 
        full_run += det_round * (rounds - 1)

        return CircuitResult(circuit=full_run)
    
    else:
        return ValueError
