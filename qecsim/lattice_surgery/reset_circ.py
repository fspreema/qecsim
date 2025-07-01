from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext

Coord = complex

__all__ = ["reset"]

def reset(*, lct : LatticeContext, patches: dict[str, Patch_Ancilla, Patch_Target, Patch_Control, Patch_Surgery], 
            cfg : Config, flow : str) -> stim.Circuit:
    
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
    z_stab_boundary_r_index_ancilla = ancilla_patch.z_bdyR
    x_stab_index_control = control_patch.x_stab
    z_stab_index_control = control_patch.z_stab
    x_stab_index_target = target_patch.x_stab
    z_stab_index_target = target_patch.z_stab

    #----------------------------------------------
    # Creating List of all Stabilizers (No Double!)
    #----------------------------------------------

    all_stabs_not_double = []
    
    #Setting double counter
    counter_x = 0
    counter_z = 0

    for index in (x_stab_index_ancilla + z_stab_index_ancilla + x_stab_index_control + z_stab_index_control + x_stab_index_target + z_stab_index_target):

        #Double Values possible
        if index in x_stab_boundary_b_index_ancilla:

            #Value already appended?
            if counter_x == 0:
                all_stabs_not_double.append(index)
                counter_x += 1

        elif index in z_stab_boundary_r_index_ancilla:

            #Value already appended?
            if counter_z == 0:
                all_stabs_not_double.append(index)
                counter_z += 1

        else:
            all_stabs_not_double.append(index)

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

    # Target
    t_log_obs_z_index : list[complex] = []
    for real in range((distance * 2) + 1, (distance * 4), 2):
        t_log_obs_z_index.append(q2i[real + 1j])

    t_log_obs_x_index : list[complex] = []
    for imag in range(1, (distance * 2), 2):
        t_log_obs_x_index.append(q2i[((distance * 2) + 1) + imag * 1j])

    # Control
    c_log_obs_z_index : list[complex] = []
    for real in range(1, (distance * 2), 2):
        c_log_obs_z_index.append(q2i[(real + ((distance * 2) + 1) * 1j)])

    c_log_obs_x_index : list[complex] = []
    for imag in range((distance * 2) + 1, (distance * 4), 2):
        c_log_obs_x_index.append(q2i[(1 + imag * 1j)])

    ########################
    # Define Initial Circuit
    ########################

    reset_circuit = stim.Circuit()

    #Appending Coords
    for q, i in q2i.items():
        reset_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

    """
    Looking at every state preperation seperatly seems to be inefficient
    ->  If not all Operators only once used one gets an incorrect formatting in the 
        timeslice view because of the Operations being in different timeslices in each TICK!
    """
    ########################################################################
    # Inilizing Ancilla in Plus (Reset) and Control/ Target in desired State
    ########################################################################

    init_patterns = {
    ("Z0", "Z0"): [("R", data_control + data_target + all_stabs_not_double)],
    ("Z0", "Z1"): [("R", data_control + data_target + all_stabs_not_double), ("X", t_log_obs_x_index)],
    ("Z0", "X+"): [("RX", data_target), ("R", data_control + all_stabs_not_double)],
    ("Z0", "X-"): [("RX", data_target), ("R", data_control + all_stabs_not_double)],
    ("Z1", "Z0"): [("R", data_control + data_target + all_stabs_not_double), ("X", c_log_obs_x_index)],
    ("Z1", "Z1"): [("R", data_control + data_target + all_stabs_not_double), ("X", c_log_obs_x_index + t_log_obs_x_index)],
    ("Z1", "X+"): [("RX", data_target), ("R", data_control + all_stabs_not_double), ("X", c_log_obs_x_index)],
    ("Z1", "X-"): [("RX", data_target), ("R", data_control + all_stabs_not_double), ("X", c_log_obs_x_index), ("Z", t_log_obs_z_index)],
    ("X+", "Z0"): [("RX", data_control), ("R", data_target + all_stabs_not_double)],
    ("X+", "Z1"): [("RX", data_control), ("R", data_target + all_stabs_not_double), ("X", t_log_obs_x_index), ("Z", c_log_obs_z_index)],
    ("X+", "X+"): [("RX", data_control + data_target), ("R", all_stabs_not_double)],
    ("X+", "X-"): [("RX", data_control + data_target), ("R", all_stabs_not_double), ("Z", t_log_obs_z_index)],
    ("X-", "Z0"): [("RX", data_control), ("R", data_target + all_stabs_not_double), ("Z", c_log_obs_z_index)],
    ("X-", "Z1"): [("RX", data_control), ("R", data_target + all_stabs_not_double), ("X", t_log_obs_x_index), ("Z", c_log_obs_z_index)],
    ("X-", "X+"): [("RX", data_control + data_target), ("R", all_stabs_not_double), ("Z", c_log_obs_z_index)],
    ("X-", "X-"): [("RX", data_control + data_target), ("R", all_stabs_not_double), ("Z", c_log_obs_z_index + t_log_obs_z_index)],
    }

    # Apply the initialization pattern
    key = (control_state_init, target_state_init)

    if key not in init_patterns:
        raise ValueError(f"Invalid basis combination: {key}")

    for gate, qubits in init_patterns[key]:
        reset_circuit.append(gate, qubits)

    reset_circuit.append("TICK")

    ##############################
    # Defining Logical Observables
    ##############################
    """
    if flow == "X -> XX":

        if control_state_init in {"X+", "X-"}:
            if target_state_init in {"X+", "X-"}:
                
                # Control stabilized by x logical
                log_x_c = []

                for imag in range(((distance * 2) + 1), (distance * 4), 2):
                    log_x_c.append(q2i[1 + imag*1j])

                #Rewriting in correct form i.e. X1 X2 etc...
                targets_c = [f"X{i}" for i in log_x_c]

                reset_circuit.append("OBSERVABLE_INCLUDE", targets_c, 0)

    elif flow == "X -> X":

        if control_state_init in {"Z0", "Z1", "X+", "X-"}:        
            if target_state_init in {"X+", "X-"}:

                # Control stabilized by x logical
                log_x_t = []

                for imag in range(1, (distance * 2), 2):
                    log_x_t.append(q2i[((distance * 2) + 1) + imag*1j])

                #Rewriting in correct form i.e. X1 X2 etc...
                targets_t = [f"X{i}" for i in log_x_t]

                reset_circuit.append("OBSERVABLE_INCLUDE", targets_t, 0)

    elif flow == "Z -> ZZ":

        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1"}:

                # Control stabilized by z logical
                log_z_t = []

                for real in range(((distance * 2) + 1), (distance * 4), 2):
                    log_z_t.append(q2i[real + 1j])

                #Rewriting in correct form i.e. X1 X2 etc...
                targets_t = [f"Z{i}" for i in log_z_t]

                reset_circuit.append("OBSERVABLE_INCLUDE", targets_t, 0)
        
    elif flow == "Z -> Z":

        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1", "X+", "X-"}:

                # Control stabilized by z logical
                log_z_c = []

                for real in range(1, (distance * 2), 2):
                    log_z_c.append(q2i[real + ((distance * 2) + 1)*1j])

                #Rewriting in correct form i.e. X1 X2 etc...
                targets_c = [f"Z{i}" for i in log_z_c]

                reset_circuit.append("OBSERVABLE_INCLUDE", targets_c, 0)"""

    return reset_circuit