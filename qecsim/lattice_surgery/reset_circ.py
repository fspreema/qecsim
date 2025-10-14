from typing import Dict, Tuple, Mapping
import stim

from tqecd import annotate_detectors_automatically

from .dataclasses import Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext, NoiseModel

# Import Part of Y Surface Code Init
from qecsim.surface_code_rotated.y_initial import y_initial
from qecsim.surface_code_rotated.y_repetition_circ import y_repetition_circ
from qecsim.surface_code_rotated.y_switch import y_switch_circ
from qecsim.surface_code_rotated.stabilizers import populate_stab_to_data
from qecsim.surface_code_rotated.circuit import _add_boundary_labels
from qecsim.surface_code_rotated.data_models import Context, Config, Patch
from qecsim.lattice_surgery.geometry import build_lattice
from qecsim.surface_code_rotated.reset_circ import reset as y_reset

Coord = complex
Label = str
Index = int
Pair = tuple[Coord, Coord]

__all__ = ["reset"]

def reset(*, 
          lct: LatticeContext, 
          patches: dict[str, Patch_Ancilla, Patch_Target, Patch_Control, Patch_Surgery], 
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

    y_patterns = {
    ("Y+", "Z0"): [("R", data_target + x_stab_index_target + z_stab_index_target)],
    ("Y+", "Z1"): [("R", data_target + x_stab_index_target + z_stab_index_target), ("X", t_log_obs_x_index)],
    ("Y+", "X+"): [("RX", data_target), ("R", x_stab_index_target + z_stab_index_target)],
    ("Y+", "X-"): [("RX", data_target), ("R", x_stab_index_target + z_stab_index_target), ("Z", t_log_obs_z_index)],
    ("Y-", "Z0"): [("R", data_target + x_stab_index_target + z_stab_index_target)],
    ("Y-", "Z1"): [("R", data_target + x_stab_index_target + z_stab_index_target), ("X", t_log_obs_x_index)],
    ("Y-", "X+"): [("RX", data_target), ("R", x_stab_index_target + z_stab_index_target)],
    ("Y-", "X-"): [("RX", data_target), ("R", x_stab_index_target + z_stab_index_target), ("Z", t_log_obs_z_index)],
    ("X+", "Y+"): [("RX", data_control), ("R", x_stab_index_control + z_stab_index_control)],
    ("X-", "Y+"): [("RX", data_control), ("R", x_stab_index_control + z_stab_index_control), ("Z", c_log_obs_z_index)],
    ("X+", "Y-"): [("RX", data_control), ("R", x_stab_index_control + z_stab_index_control)],
    ("X-", "Y-"): [("RX", data_control), ("R", x_stab_index_control + z_stab_index_control), ("Z", c_log_obs_z_index)],
    ("Z0", "Y+"): [("R", data_control + x_stab_index_control + z_stab_index_control)],
    ("Z1", "Y+"): [("R", data_control + x_stab_index_control + z_stab_index_control), ("X", c_log_obs_x_index)],
    ("Z0", "Y-"): [("R", data_control + x_stab_index_control + z_stab_index_control)],
    ("Z1", "Y-"): [("R", data_control + x_stab_index_control + z_stab_index_control), ("X", c_log_obs_x_index)],
    }

    # Apply the initialization pattern
    key = (control_state_init, target_state_init)

    if key in init_patterns:
        for gate, qubits in init_patterns[key]:
            reset_circuit.append(gate, qubits)

        reset_circuit.append("TICK")
        return reset_circuit

    if key in y_patterns:

        for gate, qubits in y_patterns[key]:
            reset_circuit.append(gate, qubits)

        # Building qubit Coords for control or target
        if key[1] in {"Y+", "Y-"}:
            qubit_coords: dict[Coord, Label] = build_lattice(distance, 
                                                             offset= (distance*2) + 0j, 
                                                             starting_stabilizer_x=False
                                                             )
            curr_offset = (distance*2) + 0j

        elif key[0] in {"Y+", "Y-"}:
            qubit_coords: dict[Coord, Label] = build_lattice(distance, 
                                                             offset= 0 + (distance*2) * 1j, 
                                                             starting_stabilizer_x=False
                                                             )
            curr_offset = 0 + (distance*2) * 1j

        # Define which fixed run settings should be given to the functions
        if key[1] == "Y+" or key[0] == "Y+":
            cfg_y = Config(distance= distance, 
                           state_init = "+i", 
                           obs = "Y", 
                           rounds = distance
                           )
        else:
            cfg_y = Config(distance= distance, 
                           state_init = "-i", 
                           obs = "Y", 
                           rounds = distance
                           )

        # Adding Boundary Labels
        _add_boundary_labels(distance= distance,
                             offset= curr_offset,
                             qubit_coords= qubit_coords, 
                             y_basis = True
                             )
        
        # Adding Info into Database
        stab_to_data: dict[tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords, 
                                                                             y_basis = True, 
                                                                             offset= curr_offset
                                                                             )
        
        stab_to_data_switch, stab_to_data_xcy = populate_stab_to_data(qubit_coords, 
                                                                      y_basis = True, 
                                                                      y_switch = True, 
                                                                      distance = distance, 
                                                                      offset= curr_offset
                                                                      )
        
        stab_to_data_memory: dict[tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords, 
                                                                                    y_basis = True, 
                                                                                    y_memory= True,
                                                                                    offset= curr_offset
                                                                                    )

        lct_y = Context(q2i= q2i, 
                        i2q= i2q, 
                        stab_to_data = stab_to_data, 
                        stab_to_data_modified = stab_to_data_switch, 
                        stab_to_data_modified2 = stab_to_data_xcy, 
                        stab_to_data_modified3 = stab_to_data_memory
                        )
        
        # Building Patch
        patches : dict[str, Patch] = {"patch": Patch.from_coords(qubit_coords, q2i),}

        build_y_circ = stim.Circuit()

        # Building needed circuits
        build_y_circ += y_reset(lct= lct_y, 
                                patches= patches,
                                cfg = cfg_y,    
                                skip_coords= True,
                                offset= curr_offset
                                ).circuit

        build_y_circ += y_initial(lct= lct_y, 
                        patch= patches["patch"], 
                        cfg = cfg_y,    
                        noise= noise,
                        offset= curr_offset
                        ).circuit
        
        build_y_circ += y_repetition_circ(lct= lct_y, 
                                            patch= patches["patch"], 
                                            cfg = cfg_y,
                                            offset= curr_offset, 
                                            noise= noise
                                            ).circuit
        
        build_y_circ += y_switch_circ(lct= lct_y, 
                                        patch= patches["patch"],
                                        offset= curr_offset, 
                                        cfg = cfg_y, 
                                        noise= noise
                                        ).circuit

        return reset_circuit, build_y_circ