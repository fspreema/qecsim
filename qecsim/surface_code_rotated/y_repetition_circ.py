import stim
import numpy as np
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
    distance = cfg.distance
    rounds = cfg.rounds
    log_obs = cfg.obs
    state_init = cfg.state_init

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

    ############################################
    # Define Pre-Rep Circuit -> for memory round
    ############################################

    #-----BUILDING-REPETITION-CIRC------
    before_round_circuit = stim.Circuit()
    round_circuit = stim.Circuit()

    #1) Reset/ Basis
    before_round_circuit.append("TICK")
    before_round_circuit.append("R", x_stab_index + z_stab_index)

    #-------Adding-Before-Round-Depol.-Data------------

    if noise.before_round_depol > 0 and memory_round:
        before_round_circuit.append("TICK")
        before_round_circuit.append("DEPOLARIZE1", data, noise.before_round_depol)

    #-------Continue-Circuit------------

    ############################################
    # Choosing logical state of the memory round
    ############################################
    """
    Adding logical Observable flip if instead of +i, -i was selected
    """

    if state_init == "-i" and memory_round:

        def _logical_y_indices(offset:int = 1) -> list[int]:
            # Both
            z_string =  [q2i[real + offset * 1j] for real in range(1, 2 * distance - 1, 2)]
            x_string = [q2i[offset + imag * 1j] for imag in range(1, 2 * distance - 1, 2)]
            y_string = [q2i[offset + offset * 1j]]

            return(x_string, y_string, z_string)
        
        x_string, y_string, z_string = _logical_y_indices(distance * 2 - 1)

        before_round_circuit.append("TICK")
        before_round_circuit.append("X", x_string)
        before_round_circuit.append("Y", y_string)
        before_round_circuit.append("Z", z_string)

    #-------Adding-After-Reset-Flip-Prob.------------

    if noise.after_r_flip > 0 and memory_round:
        before_round_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.after_r_flip)

    before_round_circuit.append("TICK")
    before_round_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0 and memory_round:
        before_round_circuit.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    before_round_circuit.append("TICK")

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
                        before_round_circuit.append(op, pair)
                    elif op == "DEPOLARIZE2":
                        before_round_circuit.append(op, pair, noise)
            else:
                #Adding Pair on Operation
                if op == "CX":
                    before_round_circuit.append(op, pair)
                elif op == "DEPOLARIZE2":
                    before_round_circuit.append(op, pair, noise)


    # Adding all the CX gates
    for order in ("1-CX", "2-CX", "3-CX", "4-CX"):
        _append_by_order("CX", order)
        if noise.after_c_depol_prob > 0 and memory_round:
            _append_by_order("DEPOLARIZE2", order, noise.after_c_depol_prob)
        before_round_circuit.append("TICK")

    #-------Continue-Circuit------------

    #3) Basis/ Measurement
    before_round_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0 and memory_round:
        before_round_circuit.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    before_round_circuit.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.-------

    if noise.before_m_flip_prob > 0 and memory_round:
        before_round_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.before_m_flip_prob)

    #-------Continue-Circuit----------

    before_round_circuit.append("M", x_stab_index + z_stab_index)

    #-------Continue-Circuit------------

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    before_round_circuit.append("SHIFT_COORDS", arg= (0,0,1))

    ###########################
    # Define Repetition Circuit
    ###########################

    #1) Reset/ Basis
    round_circuit.append("TICK")
    round_circuit.append("R", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if noise.after_r_flip > 0 and memory_round:
        round_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.after_r_flip)

    round_circuit.append("TICK")
    round_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0 and memory_round:
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
        if noise.after_c_depol_prob > 0 and memory_round:
            _append_by_order("DEPOLARIZE2", order, noise.after_c_depol_prob)
        round_circuit.append("TICK")

    #-------Continue-Circuit------------

    #3) Basis/ Measurement
    round_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0 and memory_round:
        round_circuit.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    round_circuit.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.-------

    if noise.before_m_flip_prob > 0 and memory_round:
        round_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.before_m_flip_prob)

    #-------Continue-Circuit----------

    round_circuit.append("M", x_stab_index + z_stab_index)

    #-------Continue-Circuit------------

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    round_circuit.append("SHIFT_COORDS", arg= (0,0,1))

    # If FT round or Beggining Rounds
    if not memory_round:

        rep_circ = round_circuit
        rep_circ += round_circuit * int((rounds - 2) / 2)
        
        return CircuitResult(circuit=rep_circ)
    
    elif memory_round:

        rep_circ = before_round_circuit
        rep_circ += round_circuit * (rounds - 1)

        ###################################
        # Adding X/Z oBservables if chosen:
        ###################################

        def _logical_x_indices(fixed_coord) -> list[int]:
            # Vertical string at x=1 (odd grid), along imag axis
            return [q2i[fixed_coord + imag * 1j] for imag in range(1, 2 * distance, 2)]

        def _logical_z_indices(fixed_coord) -> list[int]:
            # Horizontal string at y=1, along real axis
            return [q2i[real + fixed_coord * 1j] for real in range(1, 2 * distance, 2)]

        if log_obs == "X":

            # Getting corresponding logical string and rec
            log_x = _logical_x_indices(distance * 2 - 1)

            rep_circ.append("MX", log_x)
            rep_circ.append("OBSERVABLE_INCLUDE", [f"X{index}" for index in log_x], 0)

            # For later decoding we need the measurement record postiions of the logical operator
            rec_list = [- i - 1 for i in range(len(log_x))]

            return CircuitResult(circuit= rep_circ, obs_indices= rec_list)
        
        elif log_obs == "Z":

            # Getting corresponding logical string and rec
            log_z = _logical_z_indices(distance * 2 - 1)

            rep_circ.append("MZ", log_z)
            rep_circ.append("OBSERVABLE_INCLUDE", [f"Z{index}" for index in log_z], 0)

            # For later decoding we need the measurement record postiions of the logical operator
            rec_list = [- i - 1 for i in range(len(log_z))]

            return CircuitResult(circuit= rep_circ, obs_indices= rec_list)
                
        return CircuitResult(circuit=rep_circ)        
