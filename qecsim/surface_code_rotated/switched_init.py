from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch, Context

Coord = complex

__all__ = ["switched_circ_init"]

def switched_circ_init(*, lct : Context, patches: dict[str, Patch], 
            cfg : Config, before_round_depol : float, before_m_flip_prob : float, 
            after_r_flip : float, after_c_depol_prob : float) -> stim.Circuit:
    
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
    init_state = cfg.state_init
    stab_to_data = lct.stab_to_data_flipped

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.z_stab
    z_stab_index = patch.x_stab

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

    switched_init_circ = stim.Circuit()

    # Adding H gate on data

    switched_init_circ.append("H", data)
    switched_init_circ.append("TICK")

    #-------Adding-Before-Round-Depol.-Data------------

    if before_round_depol > 0:
        switched_init_circ.append("DEPOLARIZE1", data, before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    switched_init_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
        switched_init_circ.append("DEPOLARIZE1", x_stab_index, after_c_depol_prob)

    #-------Continue-Circuit------------

    switched_init_circ.append("TICK")

    #2) CX Operations

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switched_init_circ.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "1-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                switched_init_circ.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------

    switched_init_circ.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switched_init_circ.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "2-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                switched_init_circ.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------

    switched_init_circ.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switched_init_circ.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "3-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                switched_init_circ.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------

    switched_init_circ.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switched_init_circ.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "4-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                switched_init_circ.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------
    
    switched_init_circ.append("TICK")

    #3) Basis/ Measurement
    switched_init_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
        switched_init_circ.append("DEPOLARIZE1", x_stab_index, after_c_depol_prob)

    #-------Continue-Circuit------------

    switched_init_circ.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.-------

    if before_m_flip_prob > 0:
        switched_init_circ.append("X_ERROR", x_stab_index + z_stab_index, before_m_flip_prob)

    #-------Continue-Circuit----------

    """
    Due to the fact that x and z indicies are flipped, measurement order is changed!
    """
    
    switched_init_circ.append("MR", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if after_r_flip > 0:
        switched_init_circ.append("X_ERROR", x_stab_index + z_stab_index, after_r_flip)

    #-------Continue-Circuit------------

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    switched_init_circ.append("SHIFT_COORDS", arg= (0,0,1))

    """
    For the detectors we now compare the newly formed stabilizers with the old ones on the same position!
    -> Z stabilizers (Here after switch) lays on the old x stabilizers
    -> We compare this Z stabilizer to the old X stabilizer
    -> We compare current qubit index with the same qubit index of old emasurement round
    """

    #4) Detectors
    num_measurements_repeat = len(x_stab_index + z_stab_index)

    # X -> Z Stabilizers & Z -> X Stabilizers
    for index, q_index in enumerate(x_stab_index + z_stab_index):

        # Currently Z Stab
        if q_index in z_stab_index:
            prev_tar = int(-2.5 * num_measurements_repeat + index)
            current_tar = -1 * num_measurements_repeat + index
            switched_init_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                                (i2q[q_index].real, i2q[q_index].imag, 0))
            
        # Currently Z Stab
        if q_index in x_stab_index:
            prev_tar = int(-1.5 * num_measurements_repeat + index)
            current_tar = -1 * num_measurements_repeat + index
            switched_init_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                                (i2q[q_index].real, i2q[q_index].imag, 0))


    switched_init_circ.append("TICK")

    return(switched_init_circ)
