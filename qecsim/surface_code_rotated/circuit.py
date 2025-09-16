import stim
from typing import Dict, Tuple, List, Mapping, Any
from dataclasses import dataclass

from qecsim.lattice_surgery.geometry import build_lattice
from .stabilizers import populate_stab_to_data
from .initial import initial
from .repetition_circ import repetition_circ
from .reset_circ import reset
from .final_m_circuit import final_m
from .h_switched_round import h_switched_circ
from .h_switched_init import h_switched_circ_init
from .y_initial import y_initial
from .y_repetition_circ import y_repetition_circ
from .y_switch import y_switch_circ
from .dataclasses import Config, Patch, Context

Coord = complex
Label = str
Index = int
Pair = Tuple[Coord, Coord]

__all__ = ["Rotated_Surface_Code"]

# -------------------------
# Helper Functions
# -------------------------

def _add_boundary_labels(distance: int, qubit_coords: Dict[Coord, Label], y_basis : bool = False) -> None:
    
    """
    Adds the neseccary Boundary and Surgery Stabilizers needed
    """

    max_coord = 2 * distance

    if not y_basis:
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

    else:
        # Z-boundary stabilizers
        for y in range(4, max_coord, 4):
            coord_ancilla = complex(y, max_coord)
            qubit_coords[coord_ancilla] = "Z-STAB-BOUND-B"

        for y in range(4, max_coord, 4):
            coord_ancilla = complex(max_coord, y)
            qubit_coords[coord_ancilla] = "Z-STAB-BOUND-R"

        # Additional after H
        for y in range(2, max_coord, 4):
            coord_ancilla = complex(max_coord, y)
            qubit_coords[coord_ancilla] = "Z-STAB-BOUND-R-H"

        # X-boundary stabilizers
        for y in range(4, max_coord, 4):
            coord_ancilla = complex(y, 0)
            qubit_coords[coord_ancilla] = "X-STAB-BOUND-U"

        for y in range(4, max_coord, 4):
            coord_ancilla = complex(0, y)
            qubit_coords[coord_ancilla] = "X-STAB-BOUND-L"

        # Additional after H
        for y in range(2, max_coord, 4):
            coord_ancilla = complex(y, 0)
            qubit_coords[coord_ancilla] = "X-STAB-BOUND-U-H"

# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------

def Rotated_Surface_Code(distance: int, rounds : int, *, state_init : str, log_obs : str, logical_H : bool = False,
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

    ####################
    # 2. Insert boundary
    ####################

    """
    As we need one x and one z edge for y init, we need to differentiate the boundray labels
    """
    
    if state_init in {"+i", "-i"}:
        _add_boundary_labels(distance, qubit_coords, y_basis = True)
    else:
        _add_boundary_labels(distance, qubit_coords, y_basis = False)

    ###############################################
    # 3. Indexing All Qubits From given Coordinates
    ###############################################

    #Indexing Qubits
    q2i: dict[complex, int] = {q: i for i, q in enumerate(
    sorted(qubit_coords, key=lambda v: (v.real, v.imag))
    )}

    #Reverse Indexing
    i2q: dict[int, complex] = {i: q for q, i in q2i.items()}

    #######################################################################################################
    # 4. Adding the Mapping from Stabilizer to Data for later CX gate implementation and add into dataclass
    #######################################################################################################

    """
    stab_to_data_switch: For the implementation fo the XCY gates for the Y basis init
    stab_to-data_flipped: For the logical H gate implementation -> Switch of X and Z stabilizers
    """

    if state_init in {"+i", "-i"}:
        stab_to_data: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords, y_basis = True)
        stab_to_data_switch, stab_to_data_xcy = populate_stab_to_data(qubit_coords, y_basis = True, y_switch = True , distance = distance)
        lct = Context(q2i= q2i, i2q= i2q, stab_to_data = stab_to_data, stab_to_data_modified = stab_to_data_switch, 
                      stab_to_data_modified2 = stab_to_data_xcy)

    elif logical_H:
        stab_to_data: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords)
        stab_to_data_flipped: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords, is_flipped = True)
        lct = Context(q2i= q2i, i2q= i2q, stab_to_data = stab_to_data, stab_to_data_modified = stab_to_data_flipped)
    else:
        stab_to_data: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords)
        lct = Context(q2i= q2i, i2q= i2q, stab_to_data = stab_to_data)

    ##########################################################
    # Adding Indexes and shared information into lct dataclass
    ##########################################################

    patches : Dict[str, Patch] = {"patch": Patch.from_coords(qubit_coords, q2i),}

    ###################################
    # 5. Building Initilization Circuit
    ###################################

    # Check whether we need Y basis initilization
    if state_init in {"+i", "-i"}:
        initial_circuit = y_initial(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                              after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)

    else:
        initial_circuit = initial(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                              after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)

    ################################
    # 6. Building repetition Circuit
    ################################

    if state_init in {"+i", "-i"}:
        repet_circ = y_repetition_circ(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                              after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)
        
        switch_circ = y_switch_circ(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                              after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)

        initial_circuit += repet_circ
        initial_circuit += switch_circ       
        
    else:
        repet_circ = repetition_circ(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                              after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)
        
        initial_circuit += repet_circ

    #############################################################
    # Implement additional Circuit if Logical H gate was selected
    #############################################################

    """
    We now rund d rounds with flipped stabilizer roles
    -> i.e. X stabilizers convert to z stabilizers and x to z
    -> Flipped the stabs_to_data formalism and changed inside the function the role of x and z stab indices
    """

    # Determining if flip is needed
    if state_init in {"0", "1"} and log_obs in {"X"} and logical_H == True:
        flip_needed = True
    elif state_init in {"+", "-"} and log_obs in {"Z"} and logical_H == True:
        flip_needed = True
    else:
        flip_needed = False

    # Adding the needed circuits
    if flip_needed is True: 

        repet_switch_init = h_switched_circ_init(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                                after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)

        repet_switched = h_switched_circ(lct = lct, patches = patches, cfg = cfg, before_round_depol = noise_depol_data_init, before_m_flip_prob = noise_measure_flip, 
                                after_r_flip = noise_after_reset, after_c_depol_prob = noise_after_clifford_depol)
        
        initial_circuit += repet_circ
        initial_circuit += repet_switch_init
        initial_circuit += repet_switched

    ###############################
    # 11. Adding State initiliztion
    ###############################

    state_init_circuit = reset(lct = lct, patches = patches, cfg = cfg)

    final_measurement = final_m(lct = lct, patches = patches, cfg = cfg, before_m_flip_prob = noise_measure_flip, is_flipped = flip_needed)

    state_init_circuit += initial_circuit

    return(state_init_circuit)

    ##########################
    # Adding final measurement
    ##########################

    if isinstance(final_measurement, tuple):

        state_init_circuit += final_measurement[0]
        return state_init_circuit, final_measurement[1]
    
    else:
        state_init_circuit += final_measurement
        return state_init_circuit


