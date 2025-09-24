import stim
from .dataclasses import Config, Patch, Context, NoiseModel, CircuitResult

Coord = complex

__all__ = ["initial"]

def initial(*, 
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
    distance = cfg.distance
    init_state = cfg.state_init
    stab_to_data = lct.stab_to_data

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab

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

    ########################
    # Define Initial Circuit
    ########################

    initial_circuit = stim.Circuit()

    initial_circuit.append("TICK")

    #-------Adding-Before-Round-Depol.-Data------------

    if noise.before_round_depol > 0:
        initial_circuit.append("DEPOLARIZE1", data, noise.before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    initial_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)
        
    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #2) CX Operations

    def _pairs_for(order: str) -> list[list[int]]:
        # Return the Pairs needed at the current order
        return [[q2i[cp[1]], q2i[cp[0]]] for cp, o in stab_to_data.items() if o == order]

    def _append_by_order(op: str, order: str, noise: float = 0.0) -> None:
        # Getting pair info
        for pair in _pairs_for(order):
            #Adding Pair on Operation
            if op == "CX":
                initial_circuit.append(op, pair)
            elif op == "DEPOLARIZE2":
                initial_circuit.append(op, pair, noise)

    # Adding all the CX gates
    for order in ("1-CX", "2-CX", "3-CX", "4-CX"):
        _append_by_order("CX", order)
        if noise.after_c_depol_prob > 0:
            _append_by_order("DEPOLARIZE2", order, noise.after_c_depol_prob)
        initial_circuit.append("TICK")

    #-------Continue-Circuit------------
    
    initial_circuit.append("TICK")

    #3) Basis/ Measurement
    initial_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.------------

    if noise.before_m_flip_prob > 0:
        initial_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.before_m_flip_prob)

    #-------Continue-Circuit------------
    
    initial_circuit.append("MR", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if noise.after_r_flip > 0:
        initial_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.after_r_flip)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #4) DETECTORS -> Measure only deterministic-Stabilizers!

    if init_state in {"0", "1"}:
        num_measurements_initial = len(z_stab_index)
        
        for index, q_index in enumerate(z_stab_index):
            current_tar = -1 * num_measurements_initial + index
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    elif init_state in {"+", "-"}:
        num_measurements_initial = len(x_stab_index + z_stab_index)
        
        for index, q_index in enumerate(x_stab_index):
            current_tar = -1 * num_measurements_initial + index
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    return CircuitResult(circuit=initial_circuit)
