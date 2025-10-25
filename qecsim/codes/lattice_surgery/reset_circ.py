from typing import Mapping, Union
import stim

from qecsim.core.data_models import ConfigLatticeSurgery as Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext, NoiseModel

Coord = complex
Label = str
Index = int
Pair = tuple[Coord, Coord]

__all__ = ["reset"]

def reset(*, 
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
    distance = cfg.distance
    control_state_init = cfg.control_state_init
    target_state_init = cfg.target_state_init

    #-Retrieving Data Coords
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
    ("Z0", "X-"): [("RX", data_target), ("R", data_control + all_stabs_not_double), ("Z", t_log_obs_z_index)],
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

    if key in init_patterns:
        for gate, qubits in init_patterns[key]:
            reset_circuit.append(gate, qubits)

        reset_circuit.append("TICK")
        return reset_circuit
    
    else:
        raise ValueError(f"Invalid control/target state initialization: {control_state_init}, {target_state_init}")


        
