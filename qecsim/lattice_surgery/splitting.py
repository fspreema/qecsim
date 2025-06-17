from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext
from .stabilizers import populate_stab_to_data

Coord = complex

def split(*, lct : LatticeContext, patches: dict[str, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery,], cfg : Config, split_type : str) -> stim.Circuit:

    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    ancilla_patch = patches["ancilla"]
    target_patch = patches["target"]
    control_patch = patches["control"]
    surgery_patch = patches["surgery"]

    #-Retrieving Global Infomration
    distance = cfg.distance
    q2i = lct.q2i
    i2q = lct.i2q
    control_state_init = cfg.control_state_init
    target_state_init = cfg.target_state_init

    #-Retrieving CX Gate Orders
    stab_to_data = lct.stab_to_data
    stab_to_data_surgery_ac = lct.stab_to_data_surgery_ac
    stab_to_data_surgery_at = lct.stab_to_data_surgery_at

    #-Retrieving Lattice Coords
    qubit_coords_ancilla = ancilla_patch.coords
    qubit_coords_control = control_patch.coords
    qubit_coords_target = target_patch.coords
    qubit_coords_surgery = surgery_patch.coords

    ancilla_set = set(qubit_coords_ancilla)
    target_set  = set(qubit_coords_target)
    control_set = set(qubit_coords_control)
    
    # helper to filter the big dict
    get_view = lambda region: {
        pair: order
        for pair, order in stab_to_data.items()
        if pair[0] in region or pair[1] in region
    }

    stab_to_data_ancilla = get_view(ancilla_set)
    stab_to_data_target  = get_view(target_set)
    stab_to_data_control = get_view(control_set)

    #-Retrieving Data Coords
    data_ancilla = ancilla_patch.data
    data_control = control_patch.data
    data_target = target_patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index_ancilla = ancilla_patch.x_stab
    z_stab_index_ancilla = ancilla_patch.z_stab
    x_stab_boundary_b_index_ancilla = ancilla_patch.x_bdyB
    z_stab_boundary_r_index_ancilla = ancilla_patch.z_bdyR
    x_stab_index_control = control_patch.x_stab
    z_stab_index_control = control_patch.z_stab
    x_stab_index_target = target_patch.x_stab
    z_stab_index_target = target_patch.z_stab

    ###################################################
    # Set general global settings for each merging type
    ###################################################

    if split_type == "AC":
        stab_to_data_curr_merg = stab_to_data_surgery_ac
        stab_to_data_untouched_circ = stab_to_data_target
        x_stab_index_untouched_circ = x_stab_index_target
        z_stab_index_untouched_circ = z_stab_index_target

        #Shortcut for initilized lattice states
        untouched_state_init = target_state_init
        merged_state_init = control_state_init

    elif split_type == "AT":
        stab_to_data_curr_merg = stab_to_data_surgery_at
        stab_to_data_untouched_circ = stab_to_data_control
        x_stab_index_untouched_circ = x_stab_index_control
        z_stab_index_untouched_circ = z_stab_index_control

        #Shortcut for initilized lattice states
        untouched_state_init = control_state_init
        merged_state_init = target_state_init

    else:
        return ValueError("No valid merging Type in Function selected!")

    ################################
    # Define initial split Circuit
    ################################

    split_init_circuit = stim.Circuit()

    #Indexing of the additional Stabilizers included in the merging/splitting process
    z_stab_index_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "Z-STAB-SURGERY-M"]
    z_stab_boundary_l_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "Z-STAB-SURGERY-L"]
    x_stab_index_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-M"]
    x_stab_boundary_b_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-B"]

    #Adding h gate for X stabilizers -> Filtering out double coords
    combined_x_stab : list = []
    for coords in (x_stab_index_ancilla + x_stab_index_control + x_stab_index_target):
        if coords not in combined_x_stab:
            combined_x_stab.append(coords)

    split_init_circuit.append("TICK")
    split_init_circuit.append("TICK")
    split_init_circuit.append("TICK")

    split_init_circuit.append("H", combined_x_stab)

    ####################################################
    # CX Operations
    ####################################################

    split_init_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_init_circuit.append("CX", index_pairs)

    split_init_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_init_circuit.append("CX", index_pairs)

    split_init_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_init_circuit.append("CX", index_pairs)
    
    split_init_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_init_circuit.append("CX", index_pairs)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!):
    split_init_circuit.append("TICK")
    split_init_circuit.append("H", x_stab_index_ancilla)
    split_init_circuit.append("TICK")
    split_init_circuit.append("MR", x_stab_index_ancilla + z_stab_index_ancilla)
    split_init_circuit.append("TICK")
    split_init_circuit.append("H", x_stab_boundary_b_index_ancilla)
    split_init_circuit.append("TICK")

    ##########################################################################
    # Implementing Detectors for Ancilla (+ State -> X Basis is deterministic)
    ##########################################################################

    #Shifting time coords for correct dimension
    split_init_circuit.append("SHIFT_COORDS", arg=(0,0,1))

    #Determining Position in the measurement Run of only the Ancilla
    pos_to_index_ancilla_x : list = []
    pos_to_index_ancilla_z : list = []
    combined_x_stab_merging_lattices : list = []
    combined_z_stab_merging_lattices : list = []

    """
    Only Comparing the inner stabilizers and the boundarys who are not involed in the merging with one another
    """

    if split_type == "AC":

        #Coords of the combined merged lattice
        for coords in (x_stab_index_ancilla + x_stab_index_control):
            if coords not in combined_x_stab_merging_lattices:
                combined_x_stab_merging_lattices.append(coords)

        for coords in (z_stab_index_ancilla + z_stab_index_control + z_stab_boundary_l_surgery + z_stab_index_surgery):
            if coords not in combined_z_stab_merging_lattices:
                combined_z_stab_merging_lattices.append(coords)

        #Determing index of the measurement in the disconnected lattice configuration and the corresponding position in the joint merged lattice configuration
        for pos_single, index_single in enumerate(x_stab_index_ancilla + z_stab_index_ancilla):
            for pos_lattice, index_lattice in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
                if index_single in x_stab_index_ancilla:
                    if index_lattice in x_stab_index_ancilla:
                        if index_single not in x_stab_boundary_b_index_ancilla:
                            if index_lattice not in x_stab_boundary_b_index_ancilla:
                                if index_lattice == index_single:
                                    pos_to_index_ancilla_x.append([[pos_single, index_single], [pos_lattice, index_lattice]])

                elif index_single in z_stab_index_ancilla:
                    if index_lattice in z_stab_index_ancilla:
                        if index_lattice == index_single:
                            pos_to_index_ancilla_z.append([[pos_single, index_single], [pos_lattice, index_lattice]])

        #X-Stabs
        for joint_index_pos in pos_to_index_ancilla_x:
            current_tar = joint_index_pos[0][0] - len(x_stab_index_ancilla + z_stab_index_ancilla)
            previous_tar = joint_index_pos[1][0] - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

        #Z-Stabs
        for joint_index_pos in pos_to_index_ancilla_z:
            current_tar = joint_index_pos[0][0] - len(x_stab_index_ancilla + z_stab_index_ancilla)
            previous_tar = joint_index_pos[1][0] - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    elif split_type == "AT":

        #Adding h gate for X stabilizers only on merging lattices -> Filtering out double coords in big lattice
        for coords in (x_stab_index_ancilla + x_stab_index_target + x_stab_boundary_b_surgery + x_stab_index_surgery):
            if coords not in combined_x_stab_merging_lattices:
                combined_x_stab_merging_lattices.append(coords)

        for coords in (z_stab_index_ancilla + z_stab_index_target):
            if coords not in combined_z_stab_merging_lattices:
                combined_z_stab_merging_lattices.append(coords)

        #Determing index of the measurement in the disconnected lattice configuration and the corresponding position in the joint merged lattice configuration
        for pos_single, index_single in enumerate(x_stab_index_ancilla + z_stab_index_ancilla):
            for pos_lattice, index_lattice in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
                if index_single in x_stab_index_ancilla:
                    if index_lattice in x_stab_index_ancilla:
                        if index_lattice == index_single:
                            pos_to_index_ancilla_x.append([[pos_single, index_single], [pos_lattice, index_lattice]])

                elif index_single in z_stab_index_ancilla:
                    if index_lattice in z_stab_index_ancilla:
                        if index_single not in z_stab_boundary_r_index_ancilla:
                            if index_lattice not in z_stab_boundary_r_index_ancilla:
                                if index_lattice == index_single:
                                    pos_to_index_ancilla_z.append([[pos_single, index_single], [pos_lattice, index_lattice]])

        #X-Stabs
        for index_pos in pos_to_index_ancilla_x:
            current_tar = index_pos[0][0] - len(x_stab_index_ancilla + z_stab_index_ancilla)
            previous_tar = index_pos[1][0] - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

        #Z-Stabs
        for index_pos in pos_to_index_ancilla_z:
            current_tar = index_pos[0][0] - len(x_stab_index_ancilla + z_stab_index_ancilla)
            previous_tar = index_pos[1][0] - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))


    #Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "5-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_init_circuit.append("CX", index_pairs)

    split_init_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "6-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_init_circuit.append("CX", index_pairs)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    split_init_circuit.append("TICK")
    split_init_circuit.append("H", x_stab_index_control +  x_stab_index_target)
    split_init_circuit.append("TICK")
    split_init_circuit.append("MR", control_target_stabs)

    ################################################################
    # Determining Postion in the measurement Run of Target & Control
    ################################################################

    # Determeing Position with respect to current unmerged and previous merged run
    pos_to_index_control_x : list = []
    pos_to_index_control_z : list = []
    pos_to_index_target_x : list = []
    pos_to_index_target_z : list = []

    """
    Again only comparing stabs to the last round (from the merge) if they do not share the stabilizer
    -> These shared stabilizers are handled at the end seperately
    """

    if split_type == "AC":

        for pos_single, index_single in enumerate(control_target_stabs):
            for pos_lattice, index_lattice in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
                
                #Determining joint postiion in the current splitted measurement round and the prvious merged run
                #X-Stabs
                if index_single in x_stab_index_control:
                    if index_lattice in x_stab_index_control:
                        if index_single not in x_stab_boundary_b_index_ancilla:
                            if index_lattice not in x_stab_boundary_b_index_ancilla:
                                if index_single == index_lattice:
                                    pos_to_index_control_x.append([[pos_single, index_single], [pos_lattice, index_lattice]])

                #Z-Stabs
                elif index_single in z_stab_index_control:
                    if index_lattice in z_stab_index_control:
                        if index_single == index_lattice:
                            pos_to_index_control_z.append([[pos_single, index_single], [pos_lattice, index_lattice]])
        
        for pos_single, index_single in enumerate(control_target_stabs):
            for pos_lattice, index_lattice in enumerate(x_stab_index_untouched_circ + z_stab_index_untouched_circ):
                
                #Target always left alone -> joint pos of current meassurement and measurement of untouched lattice measurement
                if index_single in x_stab_index_target:
                    if index_lattice in x_stab_index_target:
                        if index_single == index_lattice:
                            pos_to_index_target_x.append([[pos_single, index_single], [pos_lattice, index_lattice]])

                elif index_single in z_stab_index_target:
                    if index_lattice in z_stab_index_target:
                        if index_single == index_lattice:
                            pos_to_index_target_z.append([[pos_single, index_single], [pos_lattice, index_lattice]])

        ####################################
        # Implementing Detectors for Control
        ####################################

        #Z-Stabs
        for joint_index_pos in pos_to_index_control_z:
            current_tar = joint_index_pos[0][0] - len(control_target_stabs)
            previous_tar = joint_index_pos[1][0] - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
        #X-Stabs
        for joint_index_pos in pos_to_index_control_x:
            current_tar = joint_index_pos[0][0] - len(control_target_stabs)
            previous_tar = joint_index_pos[1][0] - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
            
        ########################################
        # Implementing Detectors for Target
        ########################################

        #Z-Stabs
        for joint_index_pos in pos_to_index_target_z:
            current_tar = joint_index_pos[0][0] - len(control_target_stabs)
            previous_tar = joint_index_pos[1][0] - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
        #X-Stabs
        for joint_index_pos in pos_to_index_target_x:
            current_tar = joint_index_pos[0][0] - len(control_target_stabs)
            previous_tar = joint_index_pos[1][0] - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

        #################################################################
        # Implementing Detectors from the stabs wieght 4 ->  2 x weight 2
        #################################################################

    elif split_type == "AT":

        for pos_single, index_single in enumerate(control_target_stabs):
            for pos_lattice, index_lattice in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
                
                #Determining joint postiion in the current splitted measurement round and the prvious merged run
                #Z-Stabs
                if index_single in z_stab_index_target:
                    if index_lattice in z_stab_index_target:
                        if index_single not in z_stab_boundary_r_index_ancilla:
                            if index_lattice not in z_stab_boundary_r_index_ancilla:
                                if index_single == index_lattice:
                                    pos_to_index_control_z.append([[pos_single, index_single], [pos_lattice, index_lattice]])

                #X-Stabs
                elif index_single in x_stab_index_control:
                    if index_lattice in x_stab_index_control:
                        if index_single == index_lattice:
                            pos_to_index_control_x.append([[pos_single, index_single], [pos_lattice, index_lattice]])
        
        for pos_single, index_single in enumerate(control_target_stabs):
            for pos_lattice, index_lattice in enumerate(x_stab_index_untouched_circ + z_stab_index_untouched_circ):
                
                #Control always left alone -> joint pos of current meassurement and measurement of untouched lattice measurement
                if index_single in x_stab_index_control:
                    if index_lattice in x_stab_index_control:
                        if index_single == index_lattice:
                            pos_to_index_control_x.append([[pos_single, index_single], [pos_lattice, index_lattice]])

                elif index_single in z_stab_index_control:
                    if index_lattice in z_stab_index_control:
                        if index_single == index_lattice:
                            pos_to_index_control_z.append([[pos_single, index_single], [pos_lattice, index_lattice]])

        ####################################
        # Implementing Detectors for Target
        ####################################

        #Z-Stabs
        for joint_index_pos in pos_to_index_target_z:
            current_tar = joint_index_pos[0][0] - len(control_target_stabs)
            previous_tar = joint_index_pos[1][0] - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
        #X-Stabs
        for joint_index_pos in pos_to_index_target_x:
            current_tar = joint_index_pos[0][0] - len(control_target_stabs)
            previous_tar = joint_index_pos[1][0] - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
            
        ########################################
        # Implementing Detectors for Control
        ########################################

        #Z-Stabs
        for joint_index_pos in pos_to_index_control_z:
            current_tar = joint_index_pos[0][0] - len(control_target_stabs)
            previous_tar = joint_index_pos[1][0] - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
        #X-Stabs
        for joint_index_pos in pos_to_index_control_x:
            current_tar = joint_index_pos[0][0] - len(control_target_stabs)
            previous_tar = joint_index_pos[1][0] - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            q_index = joint_index_pos[0][1]
            split_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
        #################################################################
        # Implementing Detectors from the stabs wieght 4 ->  2 x weight 2
        #################################################################
    

    ###########################
    # Implementing Repeat Block
    ###########################

    split_repeat_circuit = stim.Circuit()

    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("SHIFT_COORDS", arg=(0,0,1))
    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("H", combined_x_stab)
    split_repeat_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_repeat_circuit.append("CX", index_pairs)

    split_repeat_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_repeat_circuit.append("CX", index_pairs)

    split_repeat_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_repeat_circuit.append("CX", index_pairs)
    
    split_repeat_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_repeat_circuit.append("CX", index_pairs)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!):
    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("H", x_stab_index_ancilla)
    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("MR", x_stab_index_ancilla + z_stab_index_ancilla)
    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("H", x_stab_boundary_b_index_ancilla)
    split_repeat_circuit.append("TICK")

    ##########################################################################
    # Implementing Detectors for Ancilla (+ State -> X Basis is deterministic)
    ##########################################################################

    #Determining Position in the measurement Run of only the Ancilla
    pos_to_index_ancilla_x : list = []
    pos_to_index_ancilla_z : list = []

    for pos, index in enumerate(x_stab_index_ancilla + z_stab_index_ancilla):
        if index in x_stab_index_ancilla:
            pos_to_index_ancilla_x.append([pos, index])

        elif index in z_stab_index_ancilla:
            pos_to_index_ancilla_z.append([pos, index])

    #Adding the needed Detectors
    for index_pos in pos_to_index_ancilla_x:
        current_tar = index_pos[0] - len(x_stab_index_ancilla + z_stab_index_ancilla)
        previous_tar = index_pos[0] - 2 * len(x_stab_index_ancilla + z_stab_index_ancilla) - len(control_target_stabs)
        q_index = index_pos[1]
        split_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))


    #Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "5-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_repeat_circuit.append("CX", index_pairs)

    split_repeat_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "6-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            split_repeat_circuit.append("CX", index_pairs)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("H", x_stab_index_control +  x_stab_index_target)
    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("MR", control_target_stabs)

    ################################################################
    # Determining Postion in the measurement Run of Target & Control
    ################################################################

    pos_to_index_control_x : list = []
    pos_to_index_control_z : list = []
    pos_to_index_target_x : list = []
    pos_to_index_target_z : list = []

    for pos, index in enumerate(control_target_stabs):
        if index in x_stab_index_control:
            pos_to_index_control_x.append([pos, index])

        elif index in z_stab_index_control:
            pos_to_index_control_z.append([pos, index])

        elif index in x_stab_index_target:
            pos_to_index_target_x.append([pos, index])

        elif index in z_stab_index_target:
            pos_to_index_target_z.append([pos, index])
    
    ####################################
    # Implementing Detectors for Control
    ####################################

    #Z-Basis (0/1 - state)
    if control_state_init in {"Z0", "Z1"}:

        for index_pos in pos_to_index_control_z:
            current_tar = index_pos[0] - len(control_target_stabs)
            previous_tar = index_pos[0] - 2 * len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
            q_index = index_pos[1]
            split_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    
    #X-Basis (+/- - state)
    elif control_state_init in {"X-", "X+"}:

        for index_pos in pos_to_index_control_x:
            current_tar = index_pos[0] - len(control_target_stabs)
            previous_tar = index_pos[0] - 2 * len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
            q_index = index_pos[1]
            split_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")

    ########################################
    # Implementing Detectors for Target
    ########################################

    #Z-Basis (0/1 - state)
    if target_state_init in {"Z0", "Z1"}:
  
        for index_pos in pos_to_index_target_z:
            current_tar = index_pos[0] - len(control_target_stabs)
            previous_tar = index_pos[0] - 2 * len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
            q_index = index_pos[1]
            split_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    
    #X-Basis (+/- - state)
    elif target_state_init in {"X-", "X+"}:

        for index_pos in pos_to_index_target_x:
            current_tar = index_pos[0] - len(control_target_stabs)
            previous_tar = index_pos[0] - 2 * len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
            q_index = index_pos[1]
            split_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Target Lattice")

    ##########################################
    # Adding Repeat Circ and returning circuit
    ##########################################

    split_init_circuit += split_repeat_circuit * (distance - 1)

    return split_init_circuit