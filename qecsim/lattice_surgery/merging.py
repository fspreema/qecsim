from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext
from .stabilizers import populate_stab_to_data

Coord = complex

def merge(*, lct : LatticeContext, patches: dict[str, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery,], cfg : Config, merging_type : str, before_m_flip_prob : float) -> stim.Circuit:

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

    ################################
    # Define initial merging Circuit
    ################################

    merge_init_circuit = stim.Circuit()

    #Indexing of the additional Stabilizers included in the merging/splitting process
    z_stab_index_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "Z-STAB-SURGERY-M"]
    z_stab_boundary_l_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "Z-STAB-SURGERY-L"]
    x_stab_index_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-M"]
    x_stab_boundary_b_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-B"]

    """
    We now redefine the CX implementation which now does Control + Ancilla or Target + Ancilla as one lattice!
    """

    ###################################################
    # Set general global settings for each merging type
    ###################################################

    if merging_type == "AC":
        stab_to_data_curr_merg = stab_to_data_surgery_ac
        stab_to_data_untouched_circ = stab_to_data_target
        x_stab_index_untouched_circ = x_stab_index_target
        z_stab_index_untouched_circ = z_stab_index_target

        #Shortcut for initilized lattice states
        untouched_state_init = target_state_init
        merged_state_init = control_state_init

    elif merging_type == "AT":
        stab_to_data_curr_merg = stab_to_data_surgery_at
        stab_to_data_untouched_circ = stab_to_data_control
        x_stab_index_untouched_circ = x_stab_index_control
        z_stab_index_untouched_circ = z_stab_index_control

        #Shortcut for initilized lattice states
        untouched_state_init = control_state_init
        merged_state_init = target_state_init

    else:
        return ValueError("No valid merging Type in Function selected!")

    ############################
    # CX-GATES-ANCILLA-&-CONTROL
    ############################

    combined_x_stab : list = []

    if merging_type == "AC":
        #Adding h gate for X stabilizers on all lattices -> Filtering out double coords
        for coords in (x_stab_index_ancilla + x_stab_index_control + x_stab_index_target):
            if coords not in combined_x_stab:
                combined_x_stab.append(coords)

    elif merging_type == "AT":
        #Adding h gate for X stabilizers on all lattices -> Filtering out double coords
        for coords in (x_stab_index_ancilla + x_stab_index_control + x_stab_index_target + x_stab_boundary_b_surgery + x_stab_index_surgery):
            if coords not in combined_x_stab:
                combined_x_stab.append(coords)

    merge_init_circuit.append("TICK")
    merge_init_circuit.append("H", combined_x_stab)
    merge_init_circuit.append("TICK")

    #2) CX Operations

    for coord_pairs, order in (stab_to_data_curr_merg |stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    merge_init_circuit.append("TICK")
            
    for coord_pairs, order in (stab_to_data_curr_merg |stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    merge_init_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_curr_merg |stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)
    
    merge_init_circuit.append("TICK")
        
    for coord_pairs, order in (stab_to_data_curr_merg |stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    combined_x_stab_merging_lattices : list = []
    combined_z_stab_merging_lattices : list = []

    if merging_type == "AC":

        #Adding h gate for X stabilizers only on merging lattices -> Filtering out double coords in big lattice
        for coords in (x_stab_index_ancilla + x_stab_index_control):
            if coords not in combined_x_stab_merging_lattices:
                combined_x_stab_merging_lattices.append(coords)

        for coords in (z_stab_index_ancilla + z_stab_index_control + z_stab_boundary_l_surgery + z_stab_index_surgery):
            if coords not in combined_z_stab_merging_lattices:
                combined_z_stab_merging_lattices.append(coords)

    elif merging_type == "AT":

        #Adding h gate for X stabilizers only on merging lattices -> Filtering out double coords in big lattice
        for coords in (x_stab_index_ancilla + x_stab_index_target + x_stab_boundary_b_surgery + x_stab_index_surgery):
            if coords not in combined_x_stab_merging_lattices:
                combined_x_stab_merging_lattices.append(coords)

        for coords in (z_stab_index_ancilla + z_stab_index_target):
            if coords not in combined_z_stab_merging_lattices:
                combined_z_stab_merging_lattices.append(coords)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement):
    merge_init_circuit.append("TICK")
    merge_init_circuit.append("H", combined_x_stab_merging_lattices)
    merge_init_circuit.append("TICK")

    #-------Adding measurement Flip Prob.--------------
    if before_m_flip_prob > 0:
        merge_init_circuit.append("X_ERROR", combined_z_stab_merging_lattices + combined_x_stab_merging_lattices, before_m_flip_prob)
        merge_init_circuit.append("TICK")
    #--------------------------------------------------

    merge_init_circuit.append("MR", combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
    merge_init_circuit.append("TICK")

    ###########################################################################################################
    # Adding Detectors -> Firstly Stabilizers which measurement is already known i.e. outside of merging region
    ###########################################################################################################

    #Determining Position in the measurement Run of only the Ancilla (Shared Stabilizers excluded i.e. shard boundary stabs)
    pos_to_index_ancilla_x_merge : list = []
    pos_to_index_ancilla_z_merge : list = []

    if merging_type == "AC":

        for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
            if index in x_stab_index_ancilla:
                if index not in x_stab_boundary_b_index_ancilla:
                    pos_to_index_ancilla_x_merge.append([pos, index])

            elif index in z_stab_index_ancilla:
                pos_to_index_ancilla_z_merge.append([pos, index])

    elif merging_type == "AT":

        for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
            if index in z_stab_index_ancilla:
                if index not in z_stab_boundary_r_index_ancilla:
                    pos_to_index_ancilla_z_merge.append([pos, index])

            elif index in x_stab_index_ancilla:
                pos_to_index_ancilla_x_merge.append([pos, index])

    #Determining Position in the measurement Run of only the Ancilla
    pos_to_index_ancilla_x : list = []
    pos_to_index_ancilla_z : list = []

    for pos, index in enumerate(x_stab_index_ancilla + z_stab_index_ancilla):
        if index in x_stab_index_ancilla:
            pos_to_index_ancilla_x.append([pos, index])

        elif index in z_stab_index_ancilla:
            pos_to_index_ancilla_z.append([pos, index])

    #Determining Postion in the measurement Run of Target & Control
    pos_to_index_control_x : list = []
    pos_to_index_control_z : list = []
    pos_to_index_target_x : list = []
    pos_to_index_target_z : list = []

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    for pos, index in enumerate(control_target_stabs):
        if index in x_stab_index_control:
            pos_to_index_control_x.append([pos, index])

        elif index in z_stab_index_control:
            pos_to_index_control_z.append([pos, index])

        elif index in x_stab_index_target:
            pos_to_index_target_x.append([pos, index])

        elif index in z_stab_index_target:
            pos_to_index_target_z.append([pos, index])

    #Shifting Coords for valid time-dim.
    merge_init_circuit.append("SHIFT_COORDS", arg=(0,0,1))

    #Adding the needed Detectors
    '''
    Adding x & z stabs in the ancillary lattice
    '''
    
    #X-Stabs
    for index_pos_merge in pos_to_index_ancilla_x_merge:
        current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        q_index = index_pos_merge[1]
        for index_pos in pos_to_index_ancilla_x:
            if q_index == index_pos[1]:
                previous_target = index_pos[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
                merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    #Z-Stabs
    for index_pos_merge in pos_to_index_ancilla_z_merge:
        current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        q_index = index_pos_merge[1]
        for index_pos in pos_to_index_ancilla_z:
            if q_index == index_pos[1]:
                previous_target = index_pos[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
                merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    #Adding the needed Detectors dependent of logical state of the lattice and the current lattice in merging
    pos_to_index_x_merging_lattice : list = []
    pos_to_index_z_merging_lattice : list = []

    if merging_type == "AC":

        #Determining Position in the measurement Run of only the current merging Lattice (Shared Stabilizers excluded)
        for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
            if index in x_stab_index_control:
                #Exclude Shared Ancillas i.e. below stabilizers from Ancilla
                if index not in x_stab_boundary_b_index_ancilla:
                    pos_to_index_x_merging_lattice.append([pos, index])

            elif index in z_stab_index_control:
                pos_to_index_z_merging_lattice.append([pos, index])

        #Z-Stabs
        for index_pos_merge in pos_to_index_z_merging_lattice:
            current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos_merge[1]
            for index_pos in pos_to_index_control_z:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

        #X-Stabs
        for index_pos_merge in pos_to_index_x_merging_lattice:
            current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos_merge[1]
            for index_pos in pos_to_index_control_x:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
    elif merging_type == "AT":
        
        #Determining Position in the measurement Run of only the current merging Lattice (Shared Stabilizers excluded)
        for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
            if index in z_stab_index_target:
                #Exclude Shared Ancillas i.e. right stabilizers from Ancilla
                if index not in z_stab_boundary_r_index_ancilla:
                    pos_to_index_z_merging_lattice.append([pos, index])

            elif index in x_stab_index_target:
                pos_to_index_x_merging_lattice.append([pos, index])

        #Z-Stabs
        for index_pos_merge in pos_to_index_z_merging_lattice:
            current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos_merge[1]
            for index_pos in pos_to_index_target_z:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

        #X-Stabs
        for index_pos_merge in pos_to_index_x_merging_lattice:
            current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos_merge[1]
            for index_pos in pos_to_index_target_x:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    """
    We now initlize the stabilizers newly introduce through the merging; one has to distinguish two types
    1) Boundary stabilizers -> Newly formed weight 4 x stabilizers take measurement record of both old boundary weight 2 stabilizers of control and ancilla
    2) No record history -> Newly formed weight 4 z stabilizers have no old record history and therefore have an non deterministic outcome!
    """
    ###########################
    # Adding shared Stabilizers
    ###########################

    pos_to_index_shared_stabs : list = []

    if merging_type == "AC":

        #Shared X-Stabilizers
        for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
            if index in x_stab_index_control:
                #Include Shared Ancillas i.e. below stabilizers from Ancilla
                if index in x_stab_boundary_b_index_ancilla:
                    pos_to_index_shared_stabs.append([pos, index])

        #Adding the needed Detectors
        for index_pos_merge in pos_to_index_shared_stabs:
            current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos_merge[1]
            for index_pos_a in pos_to_index_ancilla_x:
                if q_index == index_pos_a[1]:
                    previous_target_ancilla = index_pos_a[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
                    for index_pos_c in pos_to_index_control_x:
                        if q_index == index_pos_c[1]:
                            previous_target_control = index_pos_c[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs) 
                            merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target_ancilla), stim.target_rec(previous_target_control)], (i2q[q_index].real, i2q[q_index].imag, 0))

    elif merging_type == "AT":

        #Shared Z-Stabilizer

        for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
            if index in z_stab_index_target:
                #Include Shared Ancillas i.e. below stabilizers from Ancilla
                if index in z_stab_boundary_r_index_ancilla:
                    pos_to_index_shared_stabs.append([pos, index])

        #Adding the needed Detectors
        for index_pos_merge in pos_to_index_shared_stabs:
            current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos_merge[1]
            for index_pos_a in pos_to_index_ancilla_z:
                if q_index == index_pos_a[1]:
                    previous_target_ancilla = index_pos_a[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
                    for index_pos_c in pos_to_index_target_z:
                        if q_index == index_pos_c[1]:
                            previous_target_target = index_pos_c[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs) 
                            merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target_ancilla), stim.target_rec(previous_target_target)], (i2q[q_index].real, i2q[q_index].imag, 0))


    #Continue CX-Implementation for unotuched lattice (As AC/AT-Lattice already has a full run)

    """
    Adding needed H-Gates for X-Stabs which are shared between merged lattice and untouched lattice
    """

    if merging_type == "AT":
        merge_init_circuit.append("H", x_stab_boundary_b_index_ancilla)
        merge_init_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "5-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    merge_init_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "6-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    merge_init_circuit.append("TICK")
    merge_init_circuit.append("H", x_stab_index_untouched_circ)
    merge_init_circuit.append("TICK")

    #-------Adding measurement Flip Prob.--------------
    if before_m_flip_prob > 0:
        merge_init_circuit.append("X_ERROR", x_stab_index_untouched_circ + z_stab_index_untouched_circ, before_m_flip_prob)
        merge_init_circuit.append("TICK")
    #--------------------------------------------------

    merge_init_circuit.append("MR", x_stab_index_untouched_circ + z_stab_index_untouched_circ)

    if merging_type == "AC":

        #Determining Position in the measurement Run of only the exluded Lattice (Excluded from merge -> Normal stabilizer measurement)
        pos_to_index_x_excluded_lattice : list = []
        pos_to_index_z_excluded_lattice : list = []

        for pos, index in enumerate(x_stab_index_target + z_stab_index_target):
            if index in x_stab_index_target:
                pos_to_index_x_excluded_lattice.append([pos, index])

            elif index in z_stab_index_target:
                pos_to_index_z_excluded_lattice.append([pos, index])

        #Z-Stabs
        for index_pos_merge in pos_to_index_z_excluded_lattice:
            current_tar = index_pos_merge[0] - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_merge[1]
            for index_pos in pos_to_index_target_z:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(x_stab_index_target + z_stab_index_target) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

        #X-Stabs
        for index_pos_merge in pos_to_index_x_excluded_lattice:
            current_tar = index_pos_merge[0] - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_merge[1]
            for index_pos in pos_to_index_target_x:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(x_stab_index_target + z_stab_index_target) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

        
    elif merging_type == "AT":

        #Determining Position in the measurement Run of only the exluded Lattice (Excluded from merge -> Normal stabilizer measurement)
        pos_to_index_x_excluded_lattice : list = []
        pos_to_index_z_excluded_lattice : list = []

        for pos, index in enumerate(x_stab_index_control + z_stab_index_control):
            if index in x_stab_index_control:
                pos_to_index_x_excluded_lattice.append([pos, index])

            elif index in z_stab_index_control:
                pos_to_index_z_excluded_lattice.append([pos, index])

        #Z-Stabs
        for index_pos_merge in pos_to_index_z_excluded_lattice:
            current_tar = index_pos_merge[0] - len(x_stab_index_control + z_stab_index_control)
            q_index = index_pos_merge[1]
            for index_pos in pos_to_index_control_z:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(x_stab_index_control + z_stab_index_control) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

        #X-Stabs
        for index_pos_merge in pos_to_index_x_excluded_lattice:
            current_tar = index_pos_merge[0] - len(x_stab_index_control + z_stab_index_control)
            q_index = index_pos_merge[1]
            for index_pos in pos_to_index_control_x:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(x_stab_index_control + z_stab_index_control) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    #########################################
    # Merging-Repeat-Circuit
    #########################################

    #Defining Repeat Circuit
    merge_round_circuit = stim.Circuit()

    #Adding repeat circuit with all stabilizers defined

    merge_round_circuit.append("SHIFT_COORDS", arg=(0,0,1))
    merge_round_circuit.append("TICK")
    merge_round_circuit.append("H", combined_x_stab)
    merge_round_circuit.append("TICK")

    #################
    # CX Operations
    #################

    for coord_pairs, order in (stab_to_data_curr_merg |stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    merge_round_circuit.append("TICK")
            
    for coord_pairs, order in (stab_to_data_curr_merg |stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    merge_round_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_curr_merg |stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)
    
    merge_round_circuit.append("TICK")
        
    for coord_pairs, order in (stab_to_data_curr_merg |stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement):
    merge_round_circuit.append("TICK")
    merge_round_circuit.append("H", combined_x_stab_merging_lattices)
    merge_round_circuit.append("TICK")

    #-------Adding measurement Flip Prob.--------------
    if before_m_flip_prob > 0:
        merge_round_circuit.append("X_ERROR", combined_z_stab_merging_lattices + combined_x_stab_merging_lattices, before_m_flip_prob)
        merge_round_circuit.append("TICK")
    #--------------------------------------------------
    
    merge_round_circuit.append("MR", combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
    merge_round_circuit.append("TICK")

    ###########################################################################################################
    # Adding Detectors -> Firstly Stabilizers which measurement is already known i.e. outside of merging region
    ###########################################################################################################

    #Adding the needed Detectors for Ancilla (X & Z-Stabs)
    for index_pos_merge in pos_to_index_ancilla_x_merge:
        current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        previous_target = index_pos_merge[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
        q_index = index_pos_merge[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    for index_pos_merge in pos_to_index_ancilla_z_merge:
        current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        previous_target = index_pos_merge[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
        q_index = index_pos_merge[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    #Adding the needed Detectors for the merged lattice(Target or Control)

    #Z-Stabs
    for index_pos_merge in pos_to_index_z_merging_lattice:
        current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        previous_target = index_pos_merge[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
        q_index = index_pos_merge[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    #X-Stabs
    for index_pos_merge in pos_to_index_x_merging_lattice:
        current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        previous_target = index_pos_merge[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
        q_index = index_pos_merge[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    #Adding the Ancialla & current merging Lattice shared Detectors
    for index_pos_merge in pos_to_index_shared_stabs:
        current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        previous_target = index_pos_merge[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
        q_index = index_pos_merge[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    #############################################################################
    # Adding newly generated Stabilizers -> non det. measurements from prev round
    #############################################################################

    pos_to_index_newly_gen_stabs : list = []
    first_MM_pos : list = []

    if merging_type == "AC":

        for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
            if index in (z_stab_index_surgery + z_stab_boundary_l_surgery):
                pos_to_index_newly_gen_stabs.append([pos, index])

    elif merging_type == "AT":

        for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
            if index in (x_stab_index_surgery + x_stab_boundary_b_surgery):
                pos_to_index_newly_gen_stabs.append([pos, index])

    #Adding the needed Detectors (Newly Z/X generated Stabs)
    for index_pos_merge in pos_to_index_newly_gen_stabs:
        current_tar = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        previous_target = index_pos_merge[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
        first_MM_pos.append(previous_target - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ))
        q_index = index_pos_merge[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

    ############################################################################################
    # Continue CX-Implementation for untouched lattice (As AC/AT-Lattice already has a full run)
    ############################################################################################

    """
    Adding needed H-Gates for X-Stabs which are shared between merged lattice and untouched lattice
    """

    if merging_type == "AT":
        merge_round_circuit.append("H", x_stab_boundary_b_index_ancilla)
        merge_round_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "5-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    merge_round_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_untouched_circ).items():
   
        #Parallel Implementation of CX
        if order == "6-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    merge_round_circuit.append("TICK")
    merge_round_circuit.append("H", x_stab_index_untouched_circ)
    merge_round_circuit.append("TICK")

    #-------Adding measurement Flip Prob.--------------
    if before_m_flip_prob > 0:
        merge_round_circuit.append("X_ERROR", x_stab_index_untouched_circ + z_stab_index_untouched_circ, before_m_flip_prob)
        merge_round_circuit.append("TICK")
    #--------------------------------------------------

    merge_round_circuit.append("MR", x_stab_index_untouched_circ + z_stab_index_untouched_circ)

    #Adding the needed Detectors for the untouched lattice
    if merging_type == "AC":

        #Z-Stabs
        for index_pos_merge in pos_to_index_z_excluded_lattice:
            current_tar = index_pos_merge[0] - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            previous_target = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - 2 * len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            q_index = index_pos_merge[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

        #X-Stabs
        for index_pos_merge in pos_to_index_x_excluded_lattice:
            current_tar = index_pos_merge[0] - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            previous_target = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - 2 * len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            q_index = index_pos_merge[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
    elif merging_type == "AT":

        #Z-Stabs
        for index_pos_merge in pos_to_index_z_excluded_lattice:
            current_tar = index_pos_merge[0] - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            previous_target = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - 2 * len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            q_index = index_pos_merge[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))

        #X-Stabs
        for index_pos_merge in pos_to_index_x_excluded_lattice:
            current_tar = index_pos_merge[0] - len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            previous_target = index_pos_merge[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - 2 * len(x_stab_index_untouched_circ + z_stab_index_untouched_circ)
            q_index = index_pos_merge[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 0))
  
    ################################
    # Redefining Logical Observables
    ################################

    merge_final_circuit = stim.Circuit()

    merge_final_circuit.append("SHIFT_COORDS", arg=(0,0,1))

    '''
    All the logical Observables who cross the lattice now needs to be updates in order to commute with all stabilizers
    -> Newly introudced stabilizers would else antcommute
    '''
    """
    if merging_type == "AC":

        if control_state_init in {"X+", "X-"}:

            log_x_c_new = []

            for imag in range(1, (distance * 2), 2):
                log_x_c_new.append(q2i[1 + imag*1j])

            #Rewriting in correct form i.e. Z1 Z2 etc...
            paulis_c = [f"X{i}" for i in log_x_c_new]

            merge_final_circuit.append("OBSERVABLE_INCLUDE", paulis_c, 0)
    """

    # AT OBSERVABLE STILL ANTICOMMUTE SOMEHOW -> CHECK DETECTOR COORDINATES, MAYBE WRONG DETECTOR DEFINED!
    
    """
    elif merging_type == "AT":

        if target_state_init in {"Z0", "Z1"}:

            log_z_t_new = []

            for real in range(1, (distance * 2), 2):
                log_z_t_new.append(q2i[real + 1j])

            #Rewriting in correct form i.e. Z1 Z2 etc...
            paulis_t = [f"Z{i}" for i in log_z_t_new]

            merge_final_circuit.append("OBSERVABLE_INCLUDE", paulis_t, 1)
    """
    #####################################################
    # Adding Circuits & Receving the MXX/MZZ Measurements
    #####################################################

    merge_init_circuit += merge_round_circuit * (distance - 1)
    merge_init_circuit += merge_final_circuit

    return merge_init_circuit