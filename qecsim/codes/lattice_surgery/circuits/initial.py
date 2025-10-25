from typing import Mapping, Union
import stim
from qecsim.core.cx_builder import cx_builder
from qecsim.core.data_models import ConfigLatticeSurgery as Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext, NoiseModel

Coord = complex

__all__ = ["initial"]

def initial(*, 
            lct: LatticeContext, 
            patches: Mapping[str, Union[Patch_Ancilla, Patch_Target, Patch_Control, Patch_Surgery]], 
            cfg: Config, 
            noise: NoiseModel) -> stim.Circuit:
    
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

    rounds = distance

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

    """
    Looking at every state preperation seperatly seems to be inefficient
    ->  If not all Operators only once used one gets an incorrect formatting in the 
        timeslice view because of the Operations being in different timeslices in each TICK!
    """
    ########################################################################
    # Inilizing Ancilla in Plus (Reset) and Control/ Target in desired State
    ########################################################################
    
    init_patterns = {
        (a, b): [("RX", data_ancilla)]
        for a in ["Z0", "Z1", "X+", "X-", "Y+", "Y-"]
        for b in ["Z0", "Z1", "X+", "X-", "Y+", "Y-"]
    }

    # Apply the initialization pattern
    key = (control_state_init, target_state_init)

    if key not in init_patterns:
        raise ValueError(f"Invalid basis combination: {key}")

    for gate, qubits in init_patterns[key]:
        initial_circuit.append(gate, qubits)

    #-------Adding Before Round Data Depol.------------
    if noise.before_round_depol > 0:
        initial_circuit.append("DEPOLARIZE1", data_ancilla + data_control + data_target, noise.before_round_depol)
    #--------------------------------------------------

    initial_circuit.append("TICK")
    
    #Adding h gate for X stabilizers -> Filtering out double coords
    combined_x_stab : list = []
    for coords in (x_stab_index_ancilla + x_stab_index_control + x_stab_index_target):
        if coords not in combined_x_stab:
            combined_x_stab.append(coords)

    initial_circuit.append("H", combined_x_stab)

    #-------Adding-After-Clifford-Depol.------------
    if noise.after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", combined_x_stab, noise.after_c_depol_prob)
    #-----------------------------------------------

    initial_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    cx_builder(q2i= q2i,
               stab_to_data= stab_to_data,
               circuit= initial_circuit,
               noise= noise)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!):
    initial_circuit.append("TICK")
    initial_circuit.append("H", x_stab_index_ancilla)

    #-------Adding-After-Clifford-Depol.------------
    if noise.after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", x_stab_index_ancilla, noise.after_c_depol_prob)
    #-----------------------------------------------

    initial_circuit.append("TICK")

    #-------Adding Measurement Flip--------------------
    if noise.before_m_flip_prob > 0:
        initial_circuit.append("X_ERROR", x_stab_index_ancilla + z_stab_index_ancilla, noise.before_m_flip_prob)
    #--------------------------------------------------

    initial_circuit.append("M", x_stab_index_ancilla + z_stab_index_ancilla)
    initial_circuit.append("TICK")
    initial_circuit.append("R", x_stab_index_ancilla + z_stab_index_ancilla)

    #-------Adding-After-Reset-Flip-Prob.------------
    if noise.after_r_flip > 0:
        initial_circuit.append("X_ERROR", x_stab_index_ancilla + z_stab_index_ancilla, noise.after_r_flip)
    #------------------------------------------------

    initial_circuit.append("TICK")
    initial_circuit.append("H", x_stab_boundary_b_index_ancilla)

    #-------Adding-After-Clifford-Depol.------------
    if noise.after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", x_stab_boundary_b_index_ancilla, noise.after_c_depol_prob)
    #-----------------------------------------------

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
    cx_builder(q2i= q2i,
               stab_to_data= stab_to_data,
               circuit= initial_circuit,
               orders= ("5-CX", "6-CX"),
               noise= noise)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    initial_circuit.append("TICK")
    initial_circuit.append("H", x_stab_index_control +  x_stab_index_target)
    initial_circuit.append("SHIFT_COORDS", arg = (0,0,1))

    #-------Adding-After-Clifford-Depol.------------
    if noise.after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", x_stab_index_control +  x_stab_index_target, noise.after_c_depol_prob)
    #-----------------------------------------------

    initial_circuit.append("TICK")

    #-------Adding Measurement Flip--------------------
    if noise.before_m_flip_prob > 0:
        initial_circuit.append("X_ERROR", control_target_stabs, noise.before_m_flip_prob)
    #--------------------------------------------------

    initial_circuit.append("M", control_target_stabs)
    initial_circuit.append("TICK")
    initial_circuit.append("R", control_target_stabs)

    #-------Adding-After-Reset-Flip-Prob.------------
    if noise.after_r_flip > 0:
        initial_circuit.append("X_ERROR", control_target_stabs, noise.after_r_flip)
    #------------------------------------------------

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

    ###########################################
    # Adding Repeat Block
    ###########################################

    initial_repeat_circuit = stim.Circuit()

    initial_repeat_circuit.append("SHIFT_COORDS", arg=(0,0,1))
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", combined_x_stab)

    #-------Adding-After-Clifford-Depol.------------
    if noise.after_c_depol_prob > 0:
        initial_repeat_circuit.append("DEPOLARIZE1", combined_x_stab, noise.after_c_depol_prob)
    #-----------------------------------------------

    initial_repeat_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    cx_builder(q2i= q2i,
               stab_to_data= stab_to_data,
               circuit= initial_repeat_circuit,
               noise= noise)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!):
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", x_stab_index_ancilla)

    #-------Adding-After-Clifford-Depol.------------
    if noise.after_c_depol_prob > 0:
        initial_repeat_circuit.append("DEPOLARIZE1", x_stab_index_ancilla, noise.after_c_depol_prob)
    #-----------------------------------------------

    initial_repeat_circuit.append("TICK")

    #-------Adding Measurement Flip--------------------
    if noise.before_m_flip_prob > 0:
        initial_repeat_circuit.append("X_ERROR", x_stab_index_ancilla + z_stab_index_ancilla, noise.before_m_flip_prob)
    #--------------------------------------------------

    initial_repeat_circuit.append("M", x_stab_index_ancilla + z_stab_index_ancilla)
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("R", x_stab_index_ancilla + z_stab_index_ancilla)

    #-------Adding-After-Reset-Flip-Prob.------------
    if noise.after_r_flip > 0:
        initial_repeat_circuit.append("X_ERROR", x_stab_index_ancilla + z_stab_index_ancilla, noise.after_r_flip)
    #------------------------------------------------

    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", x_stab_boundary_b_index_ancilla)

    #-------Adding-After-Clifford-Depol.------------
    if noise.after_c_depol_prob > 0:
        initial_repeat_circuit.append("DEPOLARIZE1", x_stab_boundary_b_index_ancilla, noise.after_c_depol_prob)
    #-----------------------------------------------

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
    cx_builder(q2i= q2i,
               stab_to_data= stab_to_data,
               circuit= initial_repeat_circuit,
               orders= ("5-CX", "6-CX"),
               noise= noise)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", x_stab_index_control +  x_stab_index_target)
    initial_repeat_circuit.append("SHIFT_COORDS", arg=(0,0,1))

    #-------Adding-After-Clifford-Depol.------------
    if noise.after_c_depol_prob > 0:
        initial_repeat_circuit.append("DEPOLARIZE1", x_stab_index_control +  x_stab_index_target, noise.after_c_depol_prob)
    #-----------------------------------------------

    initial_repeat_circuit.append("TICK")

    #-------Adding Measurement Flip--------------------
    if noise.before_m_flip_prob > 0:
        initial_repeat_circuit.append("X_ERROR", control_target_stabs, noise.before_m_flip_prob)
    #--------------------------------------------------

    initial_repeat_circuit.append("M", control_target_stabs)
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("R", control_target_stabs)

    #-------Adding-After-Reset-Flip-Prob.------------
    if noise.after_r_flip > 0:
        initial_repeat_circuit.append("X_ERROR", control_target_stabs, noise.after_r_flip)
    #------------------------------------------------

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

    initial_circuit += initial_repeat_circuit * (rounds - 1)

    return initial_circuit