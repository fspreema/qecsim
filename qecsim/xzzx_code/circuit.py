import stim
import numpy as np
from typing import Dict, Tuple
from .geometry import build_lattice
from .stabilizers import populate_stab_to_data

Coord = complex
Label = str
Index = int
Pair = Tuple[Coord, Coord]

__all__ = ["XZZX_code"]

# -------------------------
# Helper Functions
# -------------------------

def _add_boundary_labels(distance: int, ancilla: Dict[Coord, Label]) -> None:
    
    """
    Adds the neseccary Boundary and Surgery Stabilizers needed
    """

    max_coord = 2 * distance

    # Z-boundary stabilizers
    for y in range(2, max_coord, 4):
        coord_ancilla = complex(0, y)
        ancilla[coord_ancilla] = "STAB-BOUND-L-Hor"

    for y in range(4, max_coord, 4):
        coord_ancilla = complex(max_coord, y)
        ancilla[coord_ancilla] = "STAB-BOUND-R-Hor"

    # X-boundary stabilizers
    for y in range(4, max_coord, 4):
        coord_ancilla = complex(y, 0)
        ancilla[coord_ancilla] = "STAB-BOUND-A-Ver"

    for y in range(2, max_coord, 4):
        coord_ancilla = complex(y, max_coord)
        ancilla[coord_ancilla] = "STAB-BOUND-B-Ver"

def _noise_mode_creator(bias : list, after_c_custom_noise : float) -> list[list]:

    """
    Creates the noise needed fot the Curstom Pauli Channels
    -> Given a bias list and full chance of a physical_error
    -> type of bias is [b_x, b_y, b_z] with b_x + b_y + b_z = 1
    -> Returns:
        - List of 3 probabilities for PAULI_CHANNEL_1
        - List of 15 probabilities for PAULI_CHANNEL_2  
    """

    # Check if Prob under 3/4 else BLoch sphere turned inside out
    if after_c_custom_noise > 3/4:
        return 0

    if np.any(bias):
        if np.isclose(sum(bias), 1.0):
            #######################################
            # Adding Noise for single Pauli Channel
            #######################################

            after_c_p_xyz = [after_c_custom_noise * bias[i] for i in range(3)]

            ######################################
            # Adding Noise for multi Pauli Channel
            ######################################

            bx, by, bz = bias 
            single_probs = np.array([1, bx, by, bz])

            after_c_p_xyz_multi_unnorm : list = []

            # Probabilities for I{I,X,Y,Z}
            after_c_p_xyz_multi_unnorm += list(single_probs)

            # Remove II prob.
            after_c_p_xyz_multi_unnorm.pop(0)

            # Probabilities for X{I,X,Y,Z}
            after_c_p_xyz_multi_unnorm += list(single_probs * bx)

            # Probabilities for Y{I,X,Y,Z}
            after_c_p_xyz_multi_unnorm += list(single_probs * by)

            # Probabilities for Z{I,X,Y,Z}
            after_c_p_xyz_multi_unnorm += list(single_probs * bz)

            #Normalize Weights
            total = sum(after_c_p_xyz_multi_unnorm)

            after_c_p_xyz_multi = [weights * (after_c_custom_noise / total) for weights in after_c_p_xyz_multi_unnorm]

        else:
            return 0

    else:
        after_c_p_xyz = [0] * 3
        after_c_p_xyz_multi = [0] * 15

    return after_c_p_xyz, after_c_p_xyz_multi

# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------

def XZZX_code(distance: int, rounds : int, *, state_init: str, before_round_depol : float = 0.0, 
              before_round_p_xyz : list = [0,0,0], before_m_flip_prob : float = 0.0,
                    after_r_flip : float = 0.0, after_c_depol_prob : float = 0.0,
                    noise_bias : list = [0,0,0], after_c_pauli_channel_prob : float = 0.0) -> stim.Circuit:
    """
    Returns XZZX-Code circuit

    Arguments:
                -> state_init: In which basis should the lattice be initlized?

    Returns:
                -> Fully implemented XZZX-Code in stim.Circuit format
    """

    ######################################
    # Creating Noise Model with given Bias
    ######################################

    if after_c_pauli_channel_prob > 3/4:
        return ValueError("Prob too high for custom channel")
    
    if after_c_pauli_channel_prob != 0 and noise_bias == [0] * 3:
        return ValueError("Please input a prob and a bias")
    
    if not np.isclose(sum(noise_bias), 1.0):
        return ValueError("Bias does not add up to 1!")

    after_c_p_xyz, after_c_p_xyz_multi = _noise_mode_creator(bias = noise_bias, after_c_custom_noise = after_c_pauli_channel_prob)

    ###############################################################
    # 1. Build independent square patches (using geometry function)
    ###############################################################
    qubit_coords: Dict[Coord, Label] = build_lattice(distance, state_init = state_init)

    #######################################################################################
    # 2. Insert boundary & surgery labels (Only get activated in splitting/merging process)
    #######################################################################################
    _add_boundary_labels(distance, qubit_coords)

    ################################################################################
    # 3. Adding the Mapping from Stabilizer to Data for later CX gate implementation
    ################################################################################
    stab_to_data: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords)

    ###############################################
    # 4. Indexing All Qubits From given Coordinates
    ###############################################

    #Indexing Qubits
    q2i: dict[complex, int] = {q: i for i, q in enumerate(
    sorted(qubit_coords, key=lambda v: (v.real, v.imag))
    )}

    #Reverse Indexing
    i2q: dict[int, complex] = {i: q for q, i in q2i.items()}

    #------------------
    # Building Circuit
    #------------------

    #Indexing Z and X Stabilizers
    stab_index_ver = [q2i[q] for q, qtype in qubit_coords.items() if qtype in {"STAB-Ver", "STAB-BOUND-A-Ver", "STAB-BOUND-B-Ver"}]
    stab_index_hor = [q2i[q] for q, qtype in qubit_coords.items() if qtype in {"STAB-Hor", "STAB-BOUND-L-Hor", "STAB-BOUND-R-Hor"}]
    stab_index = [q2i[q] for q, qtype in qubit_coords.items() if qtype in {"STAB-Hor", "STAB-Ver", "STAB-BOUND-L-Hor", 
                                                                           "STAB-BOUND-R-Hor", "STAB-BOUND-A-Ver", "STAB-BOUND-B-Ver"}]

    #Indexing Data-Qubits
    data_Z = [q2i[q] for q, qtype in qubit_coords.items() if qtype == "DATA_Z"]
    data_X = [q2i[q] for q, qtype in qubit_coords.items() if qtype == "DATA_X"]
    data = [q2i[q] for q, qtype in qubit_coords.items() if qtype in {"DATA_X", "DATA_Z"}]

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

    """
    Looking at every state preperation seperatly seems to be inefficient
    ->  If not all Operators only once used one gets an incorrect formatting in the 
        timeslice view because of the Operations being in different timeslices in each TICK!
    """

    ##################
    # Appending Coords
    ##################

    for q, i in q2i.items():
        initial_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

    ########################################################################
    # Inilizing Ancilla in Plus (Reset) and Control/ Target in desired State
    ########################################################################
    
    init_patterns = {
    ("Ver"): [("RZ", data_Z + stab_index), ("RX", data_X)],
    ("Hor"): [("RX", data_X), ("RZ", data_Z + stab_index)],
    }
    
    # Apply the initialization pattern
    key = (state_init)

    if key not in init_patterns:
        raise ValueError(f"Invalid basis combination: {key}")

    for gate, qubits in init_patterns[key]:
        initial_circuit.append(gate, qubits)

    #-------Adding Before Round Data Depol.------------
    if before_round_depol > 0:
        initial_circuit.append("DEPOLARIZE1", data, before_round_depol)
    #--------------------------------------------------

    #-------Adding Before Round Data Depol.------------
    if np.any(before_round_p_xyz):
        initial_circuit.append("PAULI_CHANNEL_1", data, before_round_p_xyz)
    #--------------------------------------------------
    
    initial_circuit.append("TICK")
    initial_circuit.append("H", stab_index)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", stab_index, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.----
    if np.any(after_c_p_xyz):
        initial_circuit.append("PAULI_CHANNEL_1", stab_index, after_c_p_xyz)
    #-----------------------------------------------

    initial_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "1-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.-----
    if np.any(after_c_p_xyz_multi):
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "1-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("PAULI_CHANNEL_2", index_pairs, after_c_p_xyz_multi)
    #-----------------------------------------------

    initial_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CZ":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CZ", index_pairs)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "2-CZ":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.-----
    if np.any(after_c_p_xyz_multi):
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "2-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("PAULI_CHANNEL_2", index_pairs, after_c_p_xyz_multi)
    #-----------------------------------------------

    initial_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CZ":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CZ", index_pairs)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "3-CZ":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.-----
    if np.any(after_c_p_xyz_multi):
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "3-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("PAULI_CHANNEL_2", index_pairs, after_c_p_xyz_multi)
    #-----------------------------------------------
    
    initial_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "4-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.-----
    if np.any(after_c_p_xyz_multi):
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "4-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("PAULI_CHANNEL_2", index_pairs, after_c_p_xyz_multi)
    #-----------------------------------------------

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!):
    initial_circuit.append("TICK")
    initial_circuit.append("H", stab_index)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", stab_index, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.----
    if np.any(after_c_p_xyz):
        initial_circuit.append("PAULI_CHANNEL_1", stab_index, after_c_p_xyz)
    #-----------------------------------------------

    initial_circuit.append("TICK")

    #-------Adding Measurement Flip--------------------
    if before_m_flip_prob > 0:
        initial_circuit.append("X_ERROR", stab_index, before_m_flip_prob)
    #--------------------------------------------------

    initial_circuit.append("MR", stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------
    if after_r_flip > 0:
        initial_circuit.append("X_ERROR", stab_index, after_r_flip)
    #------------------------------------------------

    initial_circuit.append("TICK")

    ##########################################################################
    # Implementing Detectors for Ancilla (inital Basis is deterministic)
    ##########################################################################

    #Determining Position in the measurement Run of only the Ancilla
    pos_to_index_ver : list = []
    pos_to_index_hor : list = []

    for pos, index in enumerate(stab_index):
        if index in stab_index_ver:
            pos_to_index_ver.append([pos, index])
        elif index in stab_index_hor:
            pos_to_index_hor.append([pos, index])

    if state_init in {"Ver"}:

        #Adding the needed Detectors
        for index_pos in pos_to_index_ver:
            current_tar = index_pos[0] - len(stab_index)
            q_index = index_pos[1]
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    elif state_init in {"Hor"}:

        #Adding the needed Detectors
        for index_pos in pos_to_index_hor:
            current_tar = index_pos[0] - len(stab_index)
            q_index = index_pos[1]
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")

    ###########################################
    # Adding Repeat Block
    ###########################################

    repeat_circuit = stim.Circuit()

    repeat_circuit.append("SHIFT_COORDS", arg=(0,0,1))
    repeat_circuit.append("H", stab_index)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
        repeat_circuit.append("DEPOLARIZE1", stab_index, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.----
    if np.any(after_c_p_xyz):
        repeat_circuit.append("PAULI_CHANNEL_1", stab_index, after_c_p_xyz)
    #-----------------------------------------------

    repeat_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            repeat_circuit.append("CX", index_pairs)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "1-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                repeat_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.-----
    if np.any(after_c_p_xyz_multi):
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "1-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                repeat_circuit.append("PAULI_CHANNEL_2", index_pairs, after_c_p_xyz_multi)
    #-----------------------------------------------

    repeat_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CZ":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            repeat_circuit.append("CZ", index_pairs)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "2-CZ":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                repeat_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.-----
    if np.any(after_c_p_xyz_multi):
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "2-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                repeat_circuit.append("PAULI_CHANNEL_2", index_pairs, after_c_p_xyz_multi)
    #-----------------------------------------------

    repeat_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CZ":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            repeat_circuit.append("CZ", index_pairs)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "3-CZ":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                repeat_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.-----
    if np.any(after_c_p_xyz_multi):
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "3-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                repeat_circuit.append("PAULI_CHANNEL_2", index_pairs, after_c_p_xyz_multi)
    #-----------------------------------------------
    
    repeat_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            repeat_circuit.append("CX", index_pairs)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "4-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                repeat_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.-----
    if np.any(after_c_p_xyz_multi):
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "4-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                repeat_circuit.append("PAULI_CHANNEL_2", index_pairs, after_c_p_xyz_multi)
    #-----------------------------------------------

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!):
    repeat_circuit.append("TICK")
    repeat_circuit.append("H", stab_index)

    #-------Adding-After-Clifford-Depol.------------
    if after_c_depol_prob > 0:
        repeat_circuit.append("DEPOLARIZE1", stab_index, after_c_depol_prob)
    #-----------------------------------------------

    #-------Adding-After-Clifford-Pauli-Channel.----
    if np.any(after_c_p_xyz):
        repeat_circuit.append("PAULI_CHANNEL_1", stab_index, after_c_p_xyz)
    #-----------------------------------------------

    repeat_circuit.append("TICK")

    #-------Adding Measurement Flip--------------------
    if before_m_flip_prob > 0:
        repeat_circuit.append("X_ERROR", stab_index, before_m_flip_prob)
    #--------------------------------------------------

    repeat_circuit.append("MR", stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------
    if after_r_flip > 0:
        repeat_circuit.append("X_ERROR", stab_index, after_r_flip)
    #------------------------------------------------

    repeat_circuit.append("TICK")

    ######################################################
    # Implementing Detectors (All Basis are deterministic)
    ######################################################

    #Determining Position in the measurement Run of only the Ancilla
    pos_to_index_ver : list = []
    pos_to_index_hor : list = []

    for pos, index in enumerate(stab_index):
        if index in stab_index_ver:
            pos_to_index_ver.append([pos, index])
        elif index in stab_index_hor:
            pos_to_index_hor.append([pos, index])

    """
    As intial circuit run is completed, now we define all stabilizers in every basis
    -> Detectors on both basis are now deterministic
    """

    #Adding the needed Detectors (Z-Basis)
    for index_pos in pos_to_index_ver:
        current_tar = index_pos[0] - len(stab_index)
        previous_tar = index_pos[0] - 2 * len(stab_index)
        q_index = index_pos[1]
        repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    
    #Adding the needed Detectors (X-Basis)
    for index_pos in pos_to_index_hor:
        current_tar = index_pos[0] - len(stab_index)
        previous_tar = index_pos[0] - 2 * len(stab_index)
        q_index = index_pos[1]
        repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    
    ##########################
    # Adding final Measurement
    ##########################

    final_circuit = stim.Circuit()

    final_circuit.append("SHIFT_COORDS", arg = (0,0,1))

    #-------Adding-Before-Measurement-Flip-Prob.----------
    if before_m_flip_prob > 0:
        final_circuit.append("X_ERROR", data,before_m_flip_prob)
    #-------Continue-Circuit----------

    ###########################################
    #6) Implementing Final Measurement-Round -> Detectors compromised by last round stab. measurements & Data parity checks from final measurement
    ###########################################

    final_circuit.append("MZ", data_Z)
    final_circuit.append("MX", data_X)

    # -> Defining Data to measurement indexing
    Index_to_rec_data : dict[int, int] = {q: i for i, q in enumerate(reversed(data_Z + data_X))}
    Index_to_rec_ancilla: dict[int, int] = {q: i for i, q in enumerate(reversed(stab_index))}

    """
    Why do we only check the Z stabilizers in the last measurement round and not also the x stabilizers as we did before
    -> We need to meassure the X stabilizers in the X basis but we already meassured in the z basis (XZ do not commute)
    """

    if state_init in {"Ver"}:
    
        for q, qtype in qubit_coords.items():
            
            if qtype == "STAB-Ver":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_upper_right] - 1,
                                -Index_to_rec_data[index_lower_left] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_Z + data_X)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "STAB-BOUND-A-Ver":
                #Needed Data Qubits
                lower_right = q.real + 1 + (q.imag + 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_lower_right = q2i[lower_right]
                index_lower_left = q2i[lower_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_lower_right] - 1, -Index_to_rec_data[index_lower_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_Z + data_X)]

                #Combining the record targets
                final_record = current_record + last_record
                    
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "STAB-BOUND-B-Ver":
                #Needed Data Qubits
                upper_right = (q.real + 1) + (q.imag - 1) * 1j
                upper_left = q.real - 1 + (q.imag - 1) * 1j

                #Finding the Correct Data Index
                index_upper_right = q2i[upper_right]
                index_upper_left = q2i[upper_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_right] - 1, -Index_to_rec_data[index_upper_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_Z + data_X)]

                #Combining the record targets
                final_record = current_record + last_record
                    
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))
            

    elif state_init in {"Hor"}:
            
        for q, qtype in qubit_coords.items():

            if qtype == "STAB-Hor":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_upper_right] - 1,
                                    -Index_to_rec_data[index_lower_left] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_Z + data_X)]

                #Combining the record targets
                final_record = current_record + last_record
                    
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "STAB-BOUND-L-Hor":
                #Needed Data Qubits
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_right = q2i[upper_right]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_right] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_Z + data_X)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "STAB-BOUND-R-Hor":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_lower_left = q2i[lower_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_lower_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_Z + data_X)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

    ########################################################
    #7) Defining logical Operators (Final Measurement-Round)
    ########################################################
    
    if state_init in {"Ver"}:

        # Control stabilized by x logical
        log_ver = []

        for imag in range(1, (distance * 2), 2):
            log_ver.append(q2i[1 + imag*1j])

        tar_rec = []

        for rec_pos, index in enumerate(data_Z + data_X):
            if index in log_ver:
                tar_rec.append(rec_pos)

        final_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data_Z + data_X) + k) for k in tar_rec], 0)

    elif state_init in {"Hor"}:

        # Control stabilized by x logical
        log_hor = []

        for real in range(1, (distance * 2), 2):
            log_hor.append(q2i[real + 1j])

        tar_rec = []

        for rec_pos, index in enumerate(data_Z + data_X):
            if index in log_hor:
                tar_rec.append(rec_pos)

        final_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data_Z + data_X) + k) for k in tar_rec], 0)

    ################################
    #8) Adding all circuits together
    ################################

    initial_circuit += repeat_circuit * (rounds - 1)
    initial_circuit += final_circuit

    return initial_circuit