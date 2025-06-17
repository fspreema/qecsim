from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext

Coord = complex

__all__ = ["initial"]

def initial(*, lct : LatticeContext, patches: dict[str, Patch_Ancilla, Patch_Target, Patch_Control, Patch_Surgery], cfg : Config) -> stim.Circuit:
    
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    ancilla_patch = patches["ancilla"]
    target_patch = patches["target"]
    control_patch = patches["control"]

    #-Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    distance = cfg.distance
    control_state_init = cfg.control_state_init
    target_state_init = cfg.target_state_init
    stab_to_data = lct.stab_to_data

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

    ########################
    # Define Initial Circuit
    ########################

    initial_circuit = stim.Circuit()

    #Appending Coords
    for q, i in q2i.items():
        initial_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

    """
    Looking at every state preperation seperatly seems to be inefficient
    ->  If not all Operators only once used one gets an incorrect formatting in the 
        timeslice view because of the Operations being in different timeslices in each TICK!
    """
    ########################################################################
    # Inilizing Ancilla in Plus (Reset) and Control/ Target in desired State
    ########################################################################

    init_patterns = {
    ("Z0", "Z0"): [("RX", data_ancilla), ("R", data_control + data_target), ("Z", data_ancilla)],
    ("Z0", "Z1"): [("RX", data_ancilla), ("R", data_control + data_target), ("X", data_target), ("Z", data_ancilla)],
    ("Z0", "X+"): [("RX", data_ancilla + data_target), ("R", data_control), ("Z", data_ancilla + data_target)],
    ("Z0", "X-"): [("RX", data_ancilla + data_target), ("R", data_control), ("Z", data_ancilla)],
    ("Z1", "Z0"): [("RX", data_ancilla), ("R", data_control + data_target), ("X", data_control), ("Z", data_ancilla)],
    ("Z1", "Z1"): [("RX", data_ancilla), ("R", data_control + data_target), ("X", data_control + data_target), ("Z", data_ancilla)],
    ("Z1", "X+"): [("RX", data_ancilla + data_target), ("R", data_control), ("X", data_control), ("Z", data_ancilla + data_target)],
    ("Z1", "X-"): [("RX", data_ancilla + data_target), ("R", data_control), ("X", data_control), ("Z", data_ancilla)],
    ("X+", "Z0"): [("RX", data_ancilla + data_control), ("R", data_target), ("Z", data_ancilla + data_control)],
    ("X+", "Z1"): [("RX", data_ancilla + data_control), ("R", data_target), ("X", data_target), ("Z", data_ancilla + data_control)],
    ("X+", "X+"): [("RX", data_ancilla + data_control + data_target), ("Z", data_ancilla + data_control + data_target)],
    ("X+", "X-"): [("RX", data_ancilla + data_control + data_target), ("Z", data_ancilla + data_target)],
    ("X-", "Z0"): [("RX", data_ancilla + data_control), ("R", data_target), ("Z", data_ancilla)],
    ("X-", "Z1"): [("RX", data_ancilla + data_control), ("R", data_target), ("X", data_target), ("Z", data_ancilla)],
    ("X-", "X+"): [("RX", data_ancilla + data_control + data_target), ("Z", data_ancilla + data_control)],
    ("X-", "X-"): [("RX", data_ancilla + data_control + data_target), ("Z", data_ancilla + data_control + data_target)],
    }

    # Apply the initialization pattern
    key = (control_state_init, target_state_init)

    if key not in init_patterns:
        raise ValueError(f"Invalid basis combination: {key}")

    for gate, qubits in init_patterns[key]:
        initial_circuit.append(gate, qubits)


    initial_circuit.append("TICK")

    #Adding h gate for X stabilizers -> Filtering out double coords
    combined_x_stab : list = []
    for coords in (x_stab_index_ancilla + x_stab_index_control + x_stab_index_target):
        if coords not in combined_x_stab:
            combined_x_stab.append(coords)

    initial_circuit.append("H", combined_x_stab)
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

    initial_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)

    initial_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)
    
    initial_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!):
    initial_circuit.append("TICK")
    initial_circuit.append("H", x_stab_index_ancilla)
    initial_circuit.append("TICK")
    initial_circuit.append("MR", x_stab_index_ancilla + z_stab_index_ancilla)
    initial_circuit.append("TICK")
    initial_circuit.append("H", x_stab_boundary_b_index_ancilla)
    initial_circuit.append("TICK")

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
        q_index = index_pos[1]
        initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))


    #Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "5-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)

    initial_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "6-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    initial_circuit.append("TICK")
    initial_circuit.append("H", x_stab_index_control +  x_stab_index_target)
    initial_circuit.append("TICK")
    initial_circuit.append("MR", control_target_stabs)

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
            q_index = index_pos[1]
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    
    #X-Basis (+/- - state)
    elif control_state_init in {"X-", "X+"}:

        for index_pos in pos_to_index_control_x:
            current_tar = index_pos[0] - len(control_target_stabs)
            q_index = index_pos[1]
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")

    ########################################
    # Implementing Detectors for Target
    ########################################

    #Z-Basis (0/1 - state)
    if target_state_init in {"Z0", "Z1"}:
  
        for index_pos in pos_to_index_target_z:
            current_tar = index_pos[0] - len(control_target_stabs)
            q_index = index_pos[1]
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    
    #X-Basis (+/- - state)
    elif target_state_init in {"X-", "X+"}:

        for index_pos in pos_to_index_target_x:
            current_tar = index_pos[0] - len(control_target_stabs)
            q_index = index_pos[1]
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Target Lattice")
    
    ###########################################
    # Adding Repeat Block
    ###########################################

    initial_repeat_circuit = stim.Circuit()

    
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("SHIFT_COORDS", arg=(0,0,1))
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", combined_x_stab)
    initial_repeat_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_repeat_circuit.append("CX", index_pairs)

    initial_repeat_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_repeat_circuit.append("CX", index_pairs)

    initial_repeat_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_repeat_circuit.append("CX", index_pairs)
    
    initial_repeat_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_repeat_circuit.append("CX", index_pairs)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!):
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", x_stab_index_ancilla)
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("MR", x_stab_index_ancilla + z_stab_index_ancilla)
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", x_stab_boundary_b_index_ancilla)
    initial_repeat_circuit.append("TICK")

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

    """
    As intial circuit run is completed, now we define all stabilizers in every basis
    -> Detectors on x and z stabs uncorrelated to the inital state
    """

    #Adding the needed Detectors (X-Basis)
    for index_pos in pos_to_index_ancilla_x:
        current_tar = index_pos[0] - len(x_stab_index_ancilla + z_stab_index_ancilla)
        previous_tar = index_pos[0] - 2 * len(x_stab_index_ancilla + z_stab_index_ancilla) - len(control_target_stabs)
        q_index = index_pos[1]
        initial_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    #Adding the needed Detectors (Z-Basis)
    for index_pos in pos_to_index_ancilla_z:
        current_tar = index_pos[0] - len(x_stab_index_ancilla + z_stab_index_ancilla)
        previous_tar = index_pos[0] - 2 * len(x_stab_index_ancilla + z_stab_index_ancilla) - len(control_target_stabs)
        q_index = index_pos[1]
        initial_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))


    #Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "5-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_repeat_circuit.append("CX", index_pairs)

    initial_repeat_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "6-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_repeat_circuit.append("CX", index_pairs)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", x_stab_index_control +  x_stab_index_target)
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("MR", control_target_stabs)

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

    for index_pos in pos_to_index_control_z:
        current_tar = index_pos[0] - len(control_target_stabs)
        previous_tar = index_pos[0] - 2 * len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
        q_index = index_pos[1]
        initial_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    

    for index_pos in pos_to_index_control_x:
        current_tar = index_pos[0] - len(control_target_stabs)
        previous_tar = index_pos[0] - 2 * len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
        q_index = index_pos[1]
        initial_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    ########################################
    # Implementing Detectors for Target
    ########################################

    for index_pos in pos_to_index_target_z:
        current_tar = index_pos[0] - len(control_target_stabs)
        previous_tar = index_pos[0] - 2 * len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
        q_index = index_pos[1]
        initial_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    
    for index_pos in pos_to_index_target_x:
        current_tar = index_pos[0] - len(control_target_stabs)
        previous_tar = index_pos[0] - 2 * len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
        q_index = index_pos[1]
        initial_repeat_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    initial_circuit += initial_repeat_circuit * (distance - 1)

    return initial_circuit
