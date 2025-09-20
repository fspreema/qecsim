from typing import Dict, Tuple, Mapping
from itertools import chain
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch, Context

Coord = complex

__all__ = ["y_switch_circ"]

def y_switch_circ(*, lct : Context, patches: dict[str, Patch], 
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
    stab_to_data = lct.stab_to_data
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

    pre_switch_circ = stim.Circuit()

    #-------Adding-Before-Round-Depol.-Data------------

    if before_round_depol > 0:
        pre_switch_circ.append("DEPOLARIZE1", data, before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    pre_switch_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
        pre_switch_circ.append("DEPOLARIZE1", x_stab_index, after_c_depol_prob)

    #-------Continue-Circuit------------

    pre_switch_circ.append("TICK")

    #2) CX Operations

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            pre_switch_circ.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "1-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                pre_switch_circ.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------

    pre_switch_circ.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():

        #Parallel Implementation of CX
        if order == "2-CX":

            # Skip if either endpoint equals y_index
            if any(p == y_coords for p in coord_pairs):
                continue
                
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            pre_switch_circ.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "2-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                pre_switch_circ.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------

    pre_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            pre_switch_circ.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "3-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                pre_switch_circ.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------

    pre_switch_circ.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            pre_switch_circ.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "4-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                pre_switch_circ.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------
    
    pre_switch_circ.append("TICK")

    #3) Basis/ Measurement
    pre_switch_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
        pre_switch_circ.append("DEPOLARIZE1", x_stab_index, after_c_depol_prob)

    #-------Continue-Circuit------------

    pre_switch_circ.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.-------

    if before_m_flip_prob > 0:
        pre_switch_circ.append("X_ERROR", x_stab_index + z_stab_index, before_m_flip_prob)

    #-------Continue-Circuit----------

    pre_switch_circ.append("MZ", x_stab_index + z_stab_index)

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    pre_switch_circ.append("SHIFT_COORDS", arg = (0,0,1))

    #4) Detectors
    num_measurements_repeat = len(x_stab_index + z_stab_index)

    for index, q_index in enumerate(x_stab_index + z_stab_index):
        prev_tar = -2 * num_measurements_repeat + index
        current_tar = -1 * num_measurements_repeat + index
        pre_switch_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                             (i2q[q_index].real, i2q[q_index].imag, 0))
        
    pre_switch_circ.append("TICK")

    ############################################
    # Adding the Switch after the init of the RY
    ############################################

    switch_circ = stim.Circuit()

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

    #---------Adding-RY-Gate----------------

    switch_circ.append("RY", y_index)

    #---------Adding-Resets-Old&New-Stabs---

    switch_circ.append("RZ", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)

    #-------Adding-Before-Round-Depol.-Data------------

    if before_round_depol > 0:
        switch_circ.append("DEPOLARIZE1", data, before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    switch_circ.append("TICK")
    switch_circ.append("H", x_stab_index + r_h_stabs)

    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
        switch_circ.append("DEPOLARIZE1", x_stab_index, after_c_depol_prob)

    #-------Continue-Circuit------------

    #-----------Adding-H-&-X-DAG-Gates------

    switch_circ.append("TICK")
    switch_circ.append("H", index_h)
    switch_circ.append("SQRT_X_DAG", index_x_deg)
    switch_circ.append("TICK")

    #########################################
    # Adding the XCY gates after the H switch
    #########################################

    #2) CX Operations

    for coord_pairs, order in stab_to_data_switch_xcy.items():
   
        #Parallel Implementation of CX
        if order == "1TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("XCY", index_pairs)

    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "2TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("CX", index_pairs)

    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("CX", index_pairs)

    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3.5TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("XCY", index_pairs)

    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "4TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("CX", index_pairs)
   
    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "5TICK":

            if len(coord_pairs) == 2:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                switch_circ.append("CX", index_pairs)
            else:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                switch_circ.append("CX", index_pairs)

    #-------Continue-Circuit------------
    
    switch_circ.append("TICK")

    #3) Basis/ Measurement
    switch_circ.append("H", switch_stab_apply_h)

    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
        switch_circ.append("DEPOLARIZE1", x_stab_index, after_c_depol_prob)

    #-------Continue-Circuit------------

    switch_circ.append("TICK")

    switch_circ.append("MR", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    switch_circ.append("SHIFT_COORDS", arg = (0,0,1))

    ###########################################################
    # Adding Detectors (Newly Gen Boundary and Old Stabilizers)
    ###########################################################

    #1) Deterministic Detectors which can be build up by the old stabilizers
    """
    Here we switch analogue to logical H; this means the following: 
    -> x and z stabs from the upper half of the lattice are flipped (i.e. where the h gates where applied)
    -> X and Z are the same for the lower half

    --> As we meassure x with old z on one half (So the cx have switched roles) and x with x, z with z on the other half
        there is no need to do anything i.e. old ancilla index is compared with new ancilla index

    """


    num_measurements_repeat = len(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)
    current_mes_offset = len(r_h_stabs + u_h_stabs)

    for index, q_index in enumerate(chain(x_stab_index, z_stab_index, r_h_stabs, u_h_stabs)):

        if q_index in (x_stab_index + z_stab_index):
            prev_tar = -2 * num_measurements_repeat + current_mes_offset + index
            current_tar = -1 * num_measurements_repeat + index
            
            #switch_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
    #2) Newly generated Stabilizer -> Detectors with only one measurement record

    for index, q_index in enumerate(chain(x_stab_index, z_stab_index, r_h_stabs, u_h_stabs)):

        if q_index in (r_h_stabs + u_h_stabs):
            current_tar = -1 * num_measurements_repeat + index
            #switch_circ.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))


    full_switch = pre_switch_circ + switch_circ

    return(full_switch)
