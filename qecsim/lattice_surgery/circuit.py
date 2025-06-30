import stim
from typing import Dict, Tuple, List, Mapping, Any
from dataclasses import dataclass

from .geometry import build_lattice
from .stabilizers import populate_stab_to_data
from .reset_circ import reset
from .initial import initial
from .merging import merge
from .splitting import split
from.final_m_circuit import final_m
from .dataclasses import Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext

Coord = complex
Label = str
Index = int
Pair = Tuple[Coord, Coord]

__all__ = ["surgery_circuit"]

# -------------------------
# Helper Functions
# -------------------------

def _add_boundary_labels(distance: int,
                         ancilla: Dict[Coord, Label],
                         target: Dict[Coord, Label],
                         control: Dict[Coord, Label],
                         surgery: Dict[Coord, Label]) -> None:
    
    """
    Adds the neseccary Boundary and Surgery Stabilizers needed
    """

    max_coord = 2 * distance

    # Z-boundary stabilizers
    for y in range(2, max_coord, 4):
        coord_ancilla = complex(0, y)
        coord_target = complex(distance * 2, y + 2)
        coord_control = complex(0, y + 2 + (distance * 2))
        coord_surgery = complex(y + 2, max_coord)
        ancilla[coord_ancilla] = "Z-STAB-BOUND-L-A"
        target[coord_target] = "Z-STAB-BOUND-L-T"
        control[coord_control] = "Z-STAB-BOUND-L-C"
        surgery[coord_surgery] = "Z-STAB-SURGERY-M"

    for y in range(4, max_coord, 4):
        coord_ancilla = complex(max_coord, y)
        coord_target = complex(max_coord + (distance * 2), y - 2)
        coord_control = complex(max_coord, y - 2 + (distance * 2)) 
        ancilla[coord_ancilla] = "Z-STAB-BOUND-R-A"
        target[coord_target] = "Z-STAB-BOUND-R-T"
        control[coord_control] = "Z-STAB-BOUND-R-C"

    # X-boundary stabilizers
    for y in range(4, max_coord, 4):
        coord_ancilla = complex(y, 0)
        coord_target = complex(y - 2 + (distance * 2), 0)
        coord_control = complex(y - 2, distance * 2)
        coord_surgery = complex(max_coord, y - 2)
        ancilla[coord_ancilla] = "X-STAB-BOUND-A-A"
        target[coord_target] = "X-STAB-BOUND-A-T"
        control[coord_control] = "X-STAB-BOUND-A-C"
        surgery[coord_surgery] = "X-STAB-SURGERY-M"

    for y in range(2, max_coord, 4):
        coord_ancilla = complex(y, max_coord)
        coord_target = complex(y + 2 + (distance * 2), max_coord)
        coord_control = complex(y + 2, max_coord + (distance * 2))
        ancilla[coord_ancilla] = "X-STAB-BOUND-B-A"
        target[coord_target] = "X-STAB-BOUND-B-T"
        control[coord_control] = "X-STAB-BOUND-B-C"

    #Adding X Surgery Stabilizer Between Ancilla & Target
    surgery[complex(max_coord, max_coord)] = "X-STAB-SURGERY-B"

    #Adding Z Surgery Stabilizer Between Ancilla & Control
    surgery[complex(0, max_coord)] = "Z-STAB-SURGERY-L"

# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------

def surgery_circuit(distance: int, *, target_state_init: str, control_state_init: str, flow_observable: str,
                    noise_depol_data_init : float = 0.0, noise_measure_flip : float = 0.0,
                    noise_after_reset : float = 0.0, noise_after_clifford_depol : float = 0.0) -> stim.Circuit:
    """
    Returns the full lattice surgery circuit

    Arguments:
                -> target_state_init: In which basis should the target lattice be initlized?
                -> control_state_init: In which basis should the control lattice be initilized?

    Returns:
                -> Fully implemented CX-Gate in stim.Circuit format
    """

    #########################################
    # Input fixed run settings into dataclass
    #########################################

    cfg = Config(distance= distance, 
                 target_state_init= target_state_init, 
                 control_state_init = control_state_init)

    ###############################################################
    # 1. Build independent square patches (using geometry function)
    ###############################################################
    qubit_coords_ancilla: Dict[Coord, Label] = build_lattice(distance, offset=0+0j, starting_stabilizer_x=True)
    qubit_coords_target: Dict[Coord, Label] = build_lattice(distance, offset= (distance*2) + 0j, starting_stabilizer_x=False)
    qubit_coords_control: Dict[Coord, Label] = build_lattice(distance, offset= 0 + (distance*2) * 1j, starting_stabilizer_x=False)
    qubit_coords_surgery: Dict[Coord, Label] = {}

    #######################################################################################
    # 2. Insert boundary & surgery labels (Only get activated in splitting/merging process)
    #######################################################################################
    _add_boundary_labels(distance, qubit_coords_ancilla, qubit_coords_target, qubit_coords_control, qubit_coords_surgery)

    # Merge into different Patches
    """
    Different Patches are needed, because of different Keywords on identical Coordinates (inside dict.):
    X-Stab-Boundary-Above-Control & X-Stab-Boundary-Below-Ancilla f.ex. get Keywords for surgery stabilizers
    """

    indiv_patches = qubit_coords_ancilla | qubit_coords_control | qubit_coords_target
    full_srgy_ptch = qubit_coords_ancilla | qubit_coords_control | qubit_coords_target | qubit_coords_surgery

    ################################################################################
    # 3. Adding the Mapping from Stabilizer to Data for later CX gate implementation
    ################################################################################
    stab_to_data_ancilla: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords_ancilla, merging = False)
    stab_to_data_target: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords_target, merging = False)
    stab_to_data_control: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords_control, merging = False)
    stab_to_data_surgery_ac : Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(full_srgy_ptch, merging = True, merging_type="AC")
    stab_to_data_surgery_at : Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(full_srgy_ptch, merging = True, merging_type="AT")

    ###############################################
    # 4. Indexing All Qubits From given Coordinates
    ###############################################

    #Indexing Qubits
    q2i: dict[complex, int] = {q: i for i, q in enumerate(
    sorted(full_srgy_ptch, key=lambda v: (v.real, v.imag))
    )}

    #Reverse Indexing
    i2q: dict[int, complex] = {i: q for q, i in q2i.items()}

    ##########################################################
    # Adding Indexes and shared information into lct dataclass
    ##########################################################

    lct = LatticeContext(q2i= q2i, 
                         i2q= i2q, 
                         stab_to_data = stab_to_data_ancilla | stab_to_data_control | stab_to_data_target,
                         stab_to_data_surgery_ac= stab_to_data_surgery_ac,
                         stab_to_data_surgery_at= stab_to_data_surgery_at,
                         surgery_coords= qubit_coords_surgery)
    
    patches : Dict[str, Patch_Ancilla, Patch_Target, Patch_Control, Patch_Surgery] = {
        "ancilla": Patch_Ancilla.from_coords(qubit_coords_ancilla, q2i),
        "target": Patch_Target.from_coords(qubit_coords_target, q2i),
        "control": Patch_Control.from_coords(qubit_coords_control, q2i),
        "surgery": Patch_Surgery.from_coords(qubit_coords_surgery, q2i)
    }

    ###################################
    # 5. Building Initilization Circuit
    ###################################

    initial_circuit = initial(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                              after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)
   
    #############################################
    # 6. Building Merging Ancilla Control Circuit
    #############################################

    merged_circuit_AC = merge(lct = lct, patches = patches, cfg = cfg, merging_type="AC", before_m_flip_prob = noise_measure_flip, after_r_flip = noise_after_reset)

    ###############################################
    # 5. Building Splitting Ancilla Control Circuit
    ###############################################

    split_circuit_AC = split(lct = lct, patches = patches, cfg = cfg, split_type="AC", before_m_flip_prob = noise_measure_flip, after_r_flip = noise_after_reset)

    ############################################
    # 6. Building Merging Ancilla Target Circuit
    ############################################

    merged_circuit_AT = merge(lct = lct, patches = patches, cfg = cfg, merging_type="AT", before_m_flip_prob = noise_measure_flip, after_r_flip = noise_after_reset)

    ##############################################
    # 7. Building splitting Ancilla Target Circuit
    ##############################################

    split_circuit_AT = split(lct = lct, patches = patches, cfg = cfg, split_type="AT", before_m_flip_prob = noise_measure_flip, after_r_flip = noise_after_reset)

    ################################################################################
    # 8. Creating Clipped Circuit (Without State intilization and final measurement)
    ################################################################################

    """
    Still need to add the len(stabs...) in the AT merge if the circuit is runs after AC!!
    -> Currently only constructed to look at the detectors right behind one another in the circuit!
    -> Maybe not??
    """

    initial_circuit += merged_circuit_AC
    initial_circuit += split_circuit_AC
    initial_circuit += merged_circuit_AT
    initial_circuit += split_circuit_AT

    #######################################################################
    # 9. Retrieving final Circuit with postion of parity ZZ XX Measurements
    #######################################################################

    c_log_x = []
    t_log_x = []
    c_log_z = []
    t_log_z = []
    a_log_x = []
    a_log_z = []
    
    for imag in range(((distance * 2) + 1), distance * 4, 2):
        c_log_x.append(q2i[1 + imag * 1j])

    for real in range(1, distance * 2, 2):
        c_log_z.append(q2i[real + ((distance * 2) + 1) * 1j])

    for imag in range(1, distance * 2, 2):
        t_log_x.append(q2i[((distance * 2) + 1) + (imag * 1j)])

    for real in range(((distance * 2) + 1), distance * 4, 2):
        t_log_z.append(q2i[real + 1j])

    for imag in range(1, distance * 2, 2):
        a_log_x.append(q2i[1 + imag * 1j])

    for real in range(1, distance * 2, 2):
        a_log_z.append(q2i[real + 1j])

    ##################################
    # 10. Building logical Observables
    ##################################

    if flow_observable == "X -> XX":
        if control_state_init in {"X+", "X-"}:
            if target_state_init in {"X+", "X-"}:
            
                left = '*'.join(f"X{i}" for i in c_log_x)
                right = '*'.join(f"X{i}" for i in c_log_x + t_log_x)
                result = f"{left} -> {right}"

                (included_measurements,) = initial_circuit.solve_flow_measurements([
                stim.Flow(result),
                ])

            else:
                return ValueError("Wrong target basis for selected flow")
            
        else:
            return ValueError("Wrong control basis for selected flow")

    elif flow_observable == "X -> X":
        if control_state_init in {"X+", "X-", "Z0", "Z1"}:
            if target_state_init in {"X+", "X-"}:

                left = '*'.join(f"X{i}" for i in t_log_x)
                right = '*'.join(f"X{i}" for i in t_log_x)
                result = f"{left} -> {right}"

                (included_measurements,) = initial_circuit.solve_flow_measurements([
                stim.Flow(result),
                ])

            else:
                return ValueError("Invalid target basis for selected flow")
            
        else:
            return ValueError("Invalid control state")
        
    elif flow_observable == "Z -> ZZ":
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1"}:

                left = '*'.join(f"Z{i}" for i in t_log_z)
                right = '*'.join(f"Z{i}" for i in c_log_z + t_log_z)
                result = f"{left} -> {right}"

                (included_measurements,) = initial_circuit.solve_flow_measurements([
                stim.Flow(result),
                ])

            else:
                return ValueError("Wrong target basis for selected flow")
            
        else:
            return ValueError("Wrong control basis for selected flow")

    elif flow_observable == "Z -> Z":   
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1", "X+", "X-"}:

                left = '*'.join(f"Z{i}" for i in c_log_z)
                right = '*'.join(f"Z{i}" for i in c_log_z)
                result = f"{left} -> {right}"

                (included_measurements,) = initial_circuit.solve_flow_measurements([
                stim.Flow(result),
                ])

            else:
                return ValueError("Invalid target state")
            
        else:
            return ValueError("Wrong control basis for selected flow")

    else:
        return ValueError("Invalid Flow selected")

    ###############################
    # 11. Adding State initiliztion
    ###############################

    state_init_circuit = reset(lct = lct, patches = patches, cfg = cfg, flow = flow_observable)

    final_measurement = final_m(lct = lct, patches = patches, cfg = cfg, flow = flow_observable, before_m_flip_prob = noise_measure_flip)

    state_init_circuit += initial_circuit

    ##################################################
    # 12. Adding logical Observable given by stim.Flow
    ##################################################

    # Calculating target rec pos
    rec_pos = []

    for index in included_measurements:
        current_rec_tar = initial_circuit.num_measurements - index
        rec_pos.append(- current_rec_tar)

    state_init_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(k) for k in rec_pos], 0)

    ##########################
    # Adding final measurement
    ##########################

    state_init_circuit += final_measurement

    #return_circuit = state_init.with_inlined_feedback()

    #print(c_log_x, c_log_z)
    #print(t_log_x, t_log_z)
    #print(a_log_x, a_log_z)

    return state_init_circuit

