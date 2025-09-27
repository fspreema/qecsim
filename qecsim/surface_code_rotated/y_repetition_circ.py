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
    

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    if not memory_round:
        x_stab_index = patch.x_stab
        z_stab_index = patch.z_stab
    else:
        x_stab_index = patch.x_stab_memory
        z_stab_index = patch.z_stab_memory

    # Finding Upper right qubit index -> need to look in 2-CX
    y_index = q2i[1 + 1j]

    #------------------------------------------------------
    # Creating list of Logical X/Z string and their indices
    #------------------------------------------------------
    """
    -> Used for swithcing of the state in a given basis
    """

    # Ancilla
    a_log_obs_z_index : list[complex] = []

    for real in range(1, (distance * 2), 2):
        a_log_obs_z_index.append(q2i[real + 1j])

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
        for index, q_index in enumerate(x_stab_index + z_stab_index):
            prev_tar = -2 * num_measurements_repeat + index
            current_tar = -1 * num_measurements_repeat + index
            round_circuit.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                                (i2q[q_index].real, i2q[q_index].imag, 0))
            
        rep_circ = round_circuit * int((rounds - 2) / 2)
        
        return CircuitResult(circuit=rep_circ)
            
    else:

        # Adding one circuit without detector and then adding the detectors in the second round
        pre_round = round_circuit
        det_round = stim.Circuit()
        det_round += round_circuit

        for index, q_index in enumerate(x_stab_index + z_stab_index):
            prev_tar = -2 * num_measurements_repeat + index
            current_tar = -1 * num_measurements_repeat + index
            det_round.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                                (i2q[q_index].real, i2q[q_index].imag, 0))

        full_run = pre_round 
        full_run += det_round * (rounds - 1)

        return CircuitResult(circuit=full_run)
