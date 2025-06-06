from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext
from .stabilizers import populate_stab_to_data

Coord = complex

def merge(*, lct : LatticeContext, patches: dict[str, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery,], cfg : Config, merging_type : str) -> stim.Circuit:

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

    ##################################
    # Detect what merging type we have
    ##################################

    if merging_type == "AC":
        stab_to_data_curr_merg = stab_to_data_surgery_ac
        stab_to_data_untouched_circ = stab_to_data_target
        x_stab_index_untouched_circ = x_stab_index_target
        z_stab_index_untouched_circ = z_stab_index_target

        #Adding h gate for X stabilizers on all lattices -> Filtering out double coords
        combined_x_stab : list = []
        for coords in (x_stab_index_ancilla + x_stab_index_control + x_stab_index_target):
            if coords not in combined_x_stab:
                combined_x_stab.append(coords)

        #Adding h gate for X stabilizers only on merging lattices -> Filtering out double coords in big lattice
        combined_x_stab_merging_lattices : list = []
        for coords in (x_stab_index_ancilla + x_stab_index_control):
            if coords not in combined_x_stab_merging_lattices:
                combined_x_stab_merging_lattices.append(coords)

        combined_z_stab_merging_lattices : list = []
        for coords in (z_stab_index_ancilla + z_stab_index_control + z_stab_boundary_l_surgery + z_stab_index_surgery):
            if coords not in combined_z_stab_merging_lattices:
                combined_z_stab_merging_lattices.append(coords)

    elif merging_type == "AT":
        stab_to_data_curr_merg = stab_to_data_surgery_at
        stab_to_data_untouched_circ = stab_to_data_control
        x_stab_index_untouched_circ = x_stab_index_control
        z_stab_index_untouched_circ = z_stab_index_control

        #Adding h gate for X stabilizers on all lattices -> Filtering out double coords
        combined_x_stab : list = []
        for coords in (x_stab_index_ancilla + x_stab_index_control + x_stab_index_target + x_stab_boundary_b_surgery + x_stab_index_surgery):
            if coords not in combined_x_stab:
                combined_x_stab.append(coords)

        #Adding h gate for X stabilizers only on merging lattices -> Filtering out double coords in big lattice
        combined_x_stab_merging_lattices : list = []
        for coords in (x_stab_index_ancilla + x_stab_index_target + x_stab_boundary_b_surgery + x_stab_index_surgery):
            if coords not in combined_x_stab_merging_lattices:
                combined_x_stab_merging_lattices.append(coords)

        combined_z_stab_merging_lattices : list = []
        for coords in (z_stab_index_ancilla + z_stab_index_target):
            if coords not in combined_z_stab_merging_lattices:
                combined_z_stab_merging_lattices.append(coords)

    else:
        return ValueError("No valid merging Type in Function selected!")

    #######################################
    # CX-GATES-ANCILLA-&-CONTROL
    #######################################

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

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement):
    merge_init_circuit.append("TICK")
    merge_init_circuit.append("H", combined_x_stab_merging_lattices)
    merge_init_circuit.append("TICK")
    merge_init_circuit.append("MR", combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
    merge_init_circuit.append("TICK")

    ###########################################################################################################
    # Adding Detectors -> Firstly Stabilizers which measurement is already known i.e. outside of merging region
    ###########################################################################################################

    #Determining Position in the measurement Run of only the Ancilla (Shared Stabilizers excluded)
    pos_to_index_ancilla_x_ac : list = []
    pos_to_index_ancilla_z_ac : list = []

    for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
        if index in x_stab_index_ancilla:
            if index not in x_stab_boundary_b_index_ancilla:
                pos_to_index_ancilla_x_ac.append([pos, index])

        elif index in z_stab_index_ancilla:
            pos_to_index_ancilla_z_ac.append([pos, index])

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

    #Adding the needed Detectors
    for index_pos_ac in pos_to_index_ancilla_x_ac:
        current_tar = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        q_index = index_pos_ac[1]
        for index_pos in pos_to_index_ancilla_x:
            if q_index == index_pos[1]:
                previous_target = index_pos[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
                merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #Determining Position in the measurement Run of only the Control (Shared Stabilizers excluded)
    pos_to_index_control_x_ac : list = []
    pos_to_index_control_z_ac : list = []

    for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
        if index in x_stab_index_control:
            #Exclude Shared Ancillas i.e. below stabilizers from Ancilla
            if index not in x_stab_boundary_b_index_ancilla:
                pos_to_index_control_x_ac.append([pos, index])

        elif index in z_stab_index_control:
            pos_to_index_control_z_ac.append([pos, index])

    #Adding the needed Detectors dependent of logical state of the lattice

    #Z-Basis (0/1 - state)
    if control_state_init in {"Z0", "Z1"}:

        for index_pos_ac in pos_to_index_control_z_ac:
            current_tar = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos_ac[1]
            for index_pos in pos_to_index_control_z:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #X-Basis (+/- - state)
    elif control_state_init in {"X-", "X+"}:

        for index_pos_ac in pos_to_index_control_x_ac:
            current_tar = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            q_index = index_pos_ac[1]
            for index_pos in pos_to_index_control_x:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")
    
    """
    We now initlize the stabilizers newly introduce through the merging; one has to distinguish two types
    1) Boundary stabilizers -> Newly formed weight 4 x stabilizers take measurement record of both old boundary weight 2 stabilizers of control and ancilla
    2) No record history -> Newly formed weight 4 z stabilizers have no old record history and therefore have an non deterministic outcome!
    """
    ###########################
    # Adding shared Stabilizers
    ###########################

    pos_to_index_shared_x_stabs : list = []

    for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
        if index in x_stab_index_control:
            #Include Shared Ancillas i.e. below stabilizers from Ancilla
            if index in x_stab_boundary_b_index_ancilla:
                pos_to_index_shared_x_stabs.append([pos, index])

    #Adding the needed Detectors
    for index_pos_ac in pos_to_index_shared_x_stabs:
        current_tar = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        q_index = index_pos_ac[1]
        for index_pos_a in pos_to_index_ancilla_x:
            if q_index == index_pos_a[1]:
                previous_target_ancilla = index_pos_a[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
                for index_pos_c in pos_to_index_control_x:
                    if q_index == index_pos_c[1]:
                        previous_target_control = index_pos_c[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs) 
                        merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target_ancilla), stim.target_rec(previous_target_control)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #Continue CX-Implementation for Target (As AC-Lattice already has a full run)
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
    merge_init_circuit.append("MR", x_stab_index_untouched_circ + z_stab_index_untouched_circ)

    #Determining Position in the measurement Run of only the Target (Excluded from merge -> Normal stabilizer measurement)
    pos_to_index_target_x_ac : list = []
    pos_to_index_target_z_ac : list = []

    for pos, index in enumerate(x_stab_index_target + z_stab_index_target):
        if index in x_stab_index_target:
            pos_to_index_target_x_ac.append([pos, index])

        elif index in z_stab_index_target:
            pos_to_index_target_z_ac.append([pos, index])

    #Adding the needed Detectors dependent of logical state of the lattice

    #Z-Basis (0/1 - state)
    if target_state_init in {"Z0", "Z1"}:

        for index_pos_ac in pos_to_index_target_z_ac:
            current_tar = index_pos_ac[0] - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            for index_pos in pos_to_index_target_z:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(x_stab_index_target + z_stab_index_target) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #X-Basis (+/- - state)
    elif target_state_init in {"X-", "X+"}:

        for index_pos_ac in pos_to_index_target_x_ac:
            current_tar = index_pos_ac[0] - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            for index_pos in pos_to_index_target_x:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(x_stab_index_target + z_stab_index_target) - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")

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
    merge_round_circuit.append("MR", combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
    merge_round_circuit.append("TICK")

    #Adding Detectors -> Firstly Stabilizers which measurement is already known i.e. outside of merging region

    #Adding the needed Detectors for Ancilla
    for index_pos_ac in pos_to_index_ancilla_x_ac:
        current_tar = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_target + z_stab_index_target)
        q_index = index_pos_ac[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #Adding the needed Detectors dependent of logical state of the lattice (Control)

    #Z-Basis (0/1 - state)
    if control_state_init in {"Z0", "Z1"}:

        for index_pos_ac in pos_to_index_control_z_ac:
            current_tar = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #X-Basis (+/- - state)
    elif control_state_init in {"X-", "X+"}:

        for index_pos_ac in pos_to_index_control_x_ac:
            current_tar = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
            previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")
    
    #Adding the Ancialla Control shared Detectors
    for index_pos_ac in pos_to_index_shared_x_stabs:
        current_tar = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_target + z_stab_index_target)
        q_index = index_pos_ac[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))
                        
    #Adding newly generated Stabilizers
    pos_to_index_new_z_stabs : list = []

    for pos, index in enumerate(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices):
        if index in (z_stab_index_surgery + z_stab_boundary_l_surgery):
            pos_to_index_new_z_stabs.append([pos, index])

    #Adding the needed Detectors (Newly Z generated Stabs)
    for index_pos_ac in pos_to_index_new_z_stabs:
        current_tar = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_target + z_stab_index_target)
        q_index = index_pos_ac[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #Continue CX-Implementation for Target (As AC-Lattice already has a full run)
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
    merge_round_circuit.append("MR", x_stab_index_untouched_circ + z_stab_index_untouched_circ)

    #Adding the needed Detectors dependent of logical state of the lattice

    #Z-Basis (0/1 - state)
    if target_state_init in {"Z0", "Z1"}:

        for index_pos_ac in pos_to_index_target_z_ac:
            current_tar = index_pos_ac[0] - len(x_stab_index_target + z_stab_index_target)
            previous_target = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - 2 * len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #X-Basis (+/- - state)
    elif target_state_init in {"X-", "X+"}:

        for index_pos_ac in pos_to_index_target_x_ac:
            current_tar = index_pos_ac[0] - len(x_stab_index_target + z_stab_index_target)
            previous_target = index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - 2 * len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))
    
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")
    
    

    ####################################################################################
    # Adding Conditional X Gate on Target in correspondance to the Logical XX Observable
    ####################################################################################

    """
    The logical ZZ Observable is already defined by the newly implemented Z stabilizers on the merge -> Product of the stabilizers give measurement result
    """

    merge_final_circuit = stim.Circuit()

    observable_z_targets = []

    for index_pos_ac in pos_to_index_new_z_stabs:
        observable_z_targets.append(index_pos_ac[0] - len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices) - len(x_stab_index_target + z_stab_index_target))

    merge_final_circuit.append("TICK")

    for records in observable_z_targets:
        for data in data_control:
            merge_final_circuit.append("CX", [stim.target_rec(records), data])


    #######################################################
    # Adding conditional Z gate on target, if X_L is uneven
    #######################################################

    ###############################################
    # Return Full Merge implementation Circuit
    ###############################################

    merge_init_circuit += merge_round_circuit * (distance - 1)
    merge_init_circuit += merge_final_circuit

    return merge_init_circuit