import stim
from typing import Dict, Tuple, List, Mapping, Any
from dataclasses import dataclass

from qecsim.lattice_surgery.geometry import build_lattice
from .stabilizers import populate_stab_to_data
from .initial import initial
from .repetition_circ import repetition_circ
from .reset_circ import reset
from .final_m_circuit import final_m
from .switched_round import switched_circ
from .switched_init import switched_circ_init
from .dataclasses import Config, Patch, Context

Coord = complex
Label = str
Index = int
Pair = Tuple[Coord, Coord]

__all__ = ["Rotated_Surface_Code"]

# -------------------------
# Helper Functions
# -------------------------

def _add_boundary_labels(distance: int, qubit_coords: Dict[Coord, Label]) -> None:
    
    """
    Adds the neseccary Boundary and Surgery Stabilizers needed
    """

    max_coord = 2 * distance

    # Z-boundary stabilizers
    for y in range(2, max_coord, 4):
        coord_ancilla = complex(0, y)
        qubit_coords[coord_ancilla] = "Z-STAB-BOUND-L"

    for y in range(4, max_coord, 4):
        coord_ancilla = complex(max_coord, y)
        qubit_coords[coord_ancilla] = "Z-STAB-BOUND-R"

    # X-boundary stabilizers
    for y in range(4, max_coord, 4):
        coord_ancilla = complex(y, 0)
        qubit_coords[coord_ancilla] = "X-STAB-BOUND-U"

    for y in range(2, max_coord, 4):
        coord_ancilla = complex(y, max_coord)
        qubit_coords[coord_ancilla] = "X-STAB-BOUND-B"

# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------

def Rotated_Surface_Code(distance: int, rounds : int, *, state_init : str, log_obs : str, flow_observable: str,
                    noise_depol_data_init : float = 0.0, noise_measure_flip : float = 0.0,
                    noise_after_reset : float = 0.0, noise_after_clifford_depol : float = 0.0) -> stim.Circuit:
    
    '''
    Generates Rotated-Surface-Code

    Args: 
        Distance (int): Distance of the surface code i.e. lattice size
        Rounds (int): Number of syndrome measurement rounds per Shot
        noise (float): Probability for x y and z error in Pauli-Channel
    
    Information:
        Logical Operator is Z Operator and pre-Defined!
        Logical State: 0 -> Only z Stabilizer detectors in the first round as x detectors are non deterministc for the first round 
                            (STILL: COMPLETE MEASUREMENT)
                         -> In theory we don not even need to meassure the X stabilizers at all because we do not have phase errors 
                            (Would result in global phases which can be ignored)
        Noise-Model: Analog to Stims Circuit i.e. Full Noise Model implemented
                    -> Before Round Depolarization Data
                    -> Before Measurement Flip Probability
                    -> After Clifford Depolarization
                    -> After Reset Flip Propability

    Returns:
        stim.Circuit: Comiled Circuit in Stim format
    '''

    #########################################
    # Input fixed run settings into dataclass
    #########################################
    cfg = Config(distance= distance, state_init = state_init, obs = log_obs, rounds = rounds)

    ###############################################################
    # 1. Build independent square patches (using geometry function)
    ###############################################################
    qubit_coords: Dict[Coord, Label] = build_lattice(distance, offset=0+0j, starting_stabilizer_x=True)

    #######################################################################################
    # 2. Insert boundary & surgery labels (Only get activated in splitting/merging process)
    #######################################################################################
    _add_boundary_labels(distance, qubit_coords)

    ################################################################################
    # 3. Adding the Mapping from Stabilizer to Data for later CX gate implementation
    ################################################################################
    stab_to_data: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords, is_flipped = False)
    stab_to_data_flipped: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords, is_flipped = True)

    ###############################################
    # 4. Indexing All Qubits From given Coordinates
    ###############################################

    #Indexing Qubits
    q2i: dict[complex, int] = {q: i for i, q in enumerate(
    sorted(qubit_coords, key=lambda v: (v.real, v.imag))
    )}

    #Reverse Indexing
    i2q: dict[int, complex] = {i: q for q, i in q2i.items()}

    ##########################################################
    # Adding Indexes and shared information into lct dataclass
    ##########################################################

    lct = Context(q2i= q2i, i2q= i2q, stab_to_data = stab_to_data, stab_to_data_flipped = stab_to_data_flipped)
    patches : Dict[str, Patch] = {"patch": Patch.from_coords(qubit_coords, q2i),}

    ###################################
    # 5. Building Initilization Circuit
    ###################################

    initial_circuit = initial(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                              after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)
   
    ################################
    # 6. Building repetition Circuit
    ################################

    repet_circ = repetition_circ(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                              after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)

    ########################################################################
    # Implement additional Circuit if Init and Measure Basis is not the same
    ########################################################################

    """
    We now rund d rounds with flipped stabilizer roles
    -> i.e. X stabilizers convert to z stabilizers and x to z
    -> Flipped the stabs_to_data formalism and changed inside the function the role of x and z stab indices
    """

    # Determining if flip is needed
    if state_init in {"0", "1"} and log_obs in {"X"}:
        flip_needed = True
    elif state_init in {"+", "-"} and log_obs in {"Z"}:
        flip_needed = True
    else:
        flip_needed = False

    # Adding the needed circuits
    if flip_needed is True: 

        repet_switch_init = switched_circ_init(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                                after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)

        repet_switched = switched_circ(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                                after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)
        
        initial_circuit += repet_circ
        initial_circuit += repet_switch_init
        initial_circuit += repet_switched
        
    else:
        initial_circuit += repet_circ

    #######################################################################
    # 9. Retrieving final Circuit with postion of parity ZZ XX Measurements
    #######################################################################

    log_x = []
    log_z = []

    for imag in range(1, distance * 2, 2):
        log_x.append(q2i[1 + imag * 1j])

    for real in range(1, distance * 2, 2):
        log_z.append(q2i[real + 1j])

    ##################################
    # 10. Building logical Observables
    ##################################

    if flow_observable == "X -> Z":
        if state_init in {"+", "-"}:
            if log_obs in {"Z"}:
            
                left = '*'.join(f"X{i}" for i in log_x)
                right = '*'.join(f"Z{i}" for i in log_z)
                result = f"{left} -> {right}"

                #(included_measurements,) = initial_circuit.solve_flow_measurements([stim.Flow(result),])

            else:
                return ValueError("Wrong target basis for selected flow")
            
        else:
            return ValueError("Wrong control basis for selected flow")

    elif flow_observable == "X -> X":
        if state_init in {"+", "-"}:
            if log_obs in {"X"}:

                left = '*'.join(f"X{i}" for i in log_x)
                right = '*'.join(f"X{i}" for i in log_x)
                result = f"{left} -> {right}"

                (included_measurements,) = initial_circuit.solve_flow_measurements([
                stim.Flow(result),
                ])

            else:
                return ValueError("Invalid target basis for selected flow")
            
        else:
            return ValueError("Invalid control state")
        
    elif flow_observable == "Z -> X":
        if state_init in {"0", "1"}:
            if log_obs in {"X"}:

                left = '*'.join(f"Z{i}" for i in log_z)
                right = '*'.join(f"X{i}" for i in log_x)
                result = f"{left} -> {right}"

                (included_measurements,) = initial_circuit.solve_flow_measurements([
                stim.Flow(result),
                ])

            else:
                return ValueError("Wrong target basis for selected flow")
            
        else:
            return ValueError("Wrong control basis for selected flow")

    elif flow_observable == "Z -> Z":   
        if state_init in {"0", "1"}:
            if log_obs in {"Z"}:

                left = '*'.join(f"Z{i}" for i in log_z)
                right = '*'.join(f"Z{i}" for i in log_z)
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

    state_init_circuit = reset(lct = lct, patches = patches, cfg = cfg)

    final_measurement = final_m(lct = lct, patches = patches, cfg = cfg, flow = flow_observable, before_m_flip_prob = noise_measure_flip, is_flipped = flip_needed)

    state_init_circuit += initial_circuit

    ##################################################
    # 12. Adding logical Observable given by stim.Flow
    ##################################################

    # Calculating target rec pos
    rec_pos = []

    #for index in included_measurements:
    #    current_rec_tar = initial_circuit.num_measurements - index
    #    rec_pos.append(- current_rec_tar)

    #state_init_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(k) for k in rec_pos], 0)

    ##########################
    # Adding final measurement
    ##########################

    state_init_circuit += final_measurement

    return state_init_circuit

