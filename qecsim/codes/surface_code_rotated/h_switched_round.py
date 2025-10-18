import stim
from qecsim.core.data_models import ConfigSurface as Config, Patch, Context, NoiseModel, CircuitResult

__all__ = ["h_switched_circ"]

def h_switched_circ(*, 
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
    stab_to_data = lct.stab_to_data_modified

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.z_stab
    z_stab_index = patch.x_stab

    ###########################
    # Define Repetition Circuit
    ###########################

    #-----BUILDING-REPETITION-CIRC------

    switched_round_circ = stim.Circuit()

    #-------Adding-Before-Round-Depol.-Data------------

    if noise.before_round_depol > 0:
        switched_round_circ.append("DEPOLARIZE1", data, noise.before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    switched_round_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        switched_round_circ.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    switched_round_circ.append("TICK")

    #2) CX Operations

    def _pairs_for(order: str) -> list[list[int]]:
        # Return the Pairs needed at the current order
        return [[q2i[cp[1]], q2i[cp[0]]] for cp, o in stab_to_data.items() if o == order]

    def _append_by_order(op: str, order: str, noise: float = 0.0) -> None:
        # Getting pair info
        for pair in _pairs_for(order):
            #Adding Pair on Operation
            if op == "CX":
                switched_round_circ.append(op, pair)
            elif op == "DEPOLARIZE2":
                switched_round_circ.append(op, pair, noise)

    # Adding all the CX gates
    for order in ("1-CX", "2-CX", "3-CX", "4-CX"):
        _append_by_order("CX", order)
        if noise.after_c_depol_prob > 0:
            _append_by_order("DEPOLARIZE2", order, noise.after_c_depol_prob)
        switched_round_circ.append("TICK")

    #-------Continue-Circuit------------

    #3) Basis/ Measurement
    switched_round_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        switched_round_circ.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    switched_round_circ.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.-------

    if noise.before_m_flip_prob > 0:
        switched_round_circ.append("X_ERROR", x_stab_index + z_stab_index, noise.before_m_flip_prob)

    #-------Continue-Circuit----------

    switched_round_circ.append("MR", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if noise.after_r_flip > 0:
        switched_round_circ.append("X_ERROR", x_stab_index + z_stab_index, noise.after_r_flip)

    #-------Continue-Circuit------------

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    switched_round_circ.append("SHIFT_COORDS", arg= (0, 0, 1))

    #4) Detectors
    num_measurements_repeat = len(x_stab_index + z_stab_index)

    for index, q_index in enumerate(x_stab_index + z_stab_index):
        prev_tar = -2 * num_measurements_repeat + index
        current_tar = -1 * num_measurements_repeat + index
        switched_round_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                             (i2q[q_index].real, i2q[q_index].imag, 0))
        
    switched_round_circ.append("TICK")

    rep_circ = switched_round_circ * (rounds - 1)

    return CircuitResult(circuit=rep_circ)
