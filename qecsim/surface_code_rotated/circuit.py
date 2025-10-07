import stim

from tqecd import annotate_detectors_automatically

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
from .y_rev_switch_circ import y_rev_switch_circ
from .data_models import Config, Patch, Context, NoiseModel, CircuitResult

Coord = complex
Label = str
Index = int
Pair = tuple[Coord, Coord]

__all__ = ["rotated_surface_code"]

# -------------------------
# Helper Functions
# -------------------------

def _add_boundary_labels(distance: int, 
                         qubit_coords: dict[Coord, Label], 
                         y_basis: bool = False) -> None:
    
    """
    Adds the neseccary Boundary and Surgery Stabilizers needed for the code
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
            qubit_coords[coord_ancilla] = "X-STAB-BOUND-B"

        for y in range(4, max_coord, 4):
            coord_ancilla = complex(max_coord, y)
            qubit_coords[coord_ancilla] = "X-STAB-BOUND-R"

        # Additional after H
        for y in range(2, max_coord, 4):
            coord_ancilla = complex(max_coord, y)
            qubit_coords[coord_ancilla] = "X-STAB-BOUND-R-H"

        # X-boundary stabilizers
        for y in range(4, max_coord, 4):
            coord_ancilla = complex(y, 0)
            qubit_coords[coord_ancilla] = "Z-STAB-BOUND-U"

        for y in range(4, max_coord, 4):
            coord_ancilla = complex(0, y)
            qubit_coords[coord_ancilla] = "Z-STAB-BOUND-L"

        # Additional after H
        for y in range(2, max_coord, 4):
            coord_ancilla = complex(y, 0)
            qubit_coords[coord_ancilla] = "Z-STAB-BOUND-U-H"

def _needs_flip(state_init: str, 
                log_obs: str, 
                logical_h: bool) -> bool:
    
    # Determining if flip is needed
    return (
        (state_init in {"0", "1"} and log_obs == "X" and logical_h)
        or (state_init in {"+", "-"} and log_obs == "Z" and logical_h)
    )

# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------

def rotated_surface_code(distance: int, 
                        rounds: int, *, 
                        state_init: str, 
                        log_obs: str, 
                        logical_h: bool = False,
                        noise_depol_data_init: float = 0.0, 
                        noise_measure_flip: float = 0.0,
                        noise_after_reset: float = 0.0, 
                        noise_after_clifford_depol: float = 0.0,
                        noise_h_flip_prob: float = 0.0) -> stim.Circuit:
    
    """
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
    """

    #########################################
    # Input fixed run settings into dataclass
    #########################################
    cfg = Config(distance= distance, state_init = state_init, obs = log_obs, rounds = rounds)

    noise = NoiseModel(
        before_round_depol= noise_depol_data_init,
        before_m_flip_prob= noise_measure_flip,
        after_r_flip= noise_after_reset,
        after_c_depol_prob= noise_after_clifford_depol
    )

    ###############################################################
    # 1. Build independent square patches (using geometry function)
    ###############################################################

    # Is logical Basis y?
    is_y = state_init in {"+i", "-i"}

    if is_y:
        qubit_coords: dict[Coord, Label] = build_lattice(distance, offset=0+0j, starting_stabilizer_x=False)
    else: 
        qubit_coords: dict[Coord, Label] = build_lattice(distance, offset=0+0j, starting_stabilizer_x=True)

    ####################
    # 2. Insert boundary
    ####################

    """
    As we need one x and one z edge for y init, we need to differentiate the boundray labels
    """
    
    if is_y:
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

    if is_y:
        stab_to_data: dict[tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords, y_basis = True)
        stab_to_data_switch, stab_to_data_xcy = populate_stab_to_data(qubit_coords, y_basis = True, y_switch = True, distance = distance)
        stab_to_data_memory: dict[tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords, y_basis = True, y_memory= True)
        lct = Context(q2i= q2i, i2q= i2q, stab_to_data = stab_to_data, stab_to_data_modified = stab_to_data_switch, 
                      stab_to_data_modified2 = stab_to_data_xcy, stab_to_data_modified3 = stab_to_data_memory)

    elif logical_h:
        stab_to_data: dict[tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords)
        stab_to_data_flipped: dict[tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords, is_flipped = True)
        lct = Context(q2i= q2i, i2q= i2q, stab_to_data = stab_to_data, stab_to_data_modified = stab_to_data_flipped)
    else:
        stab_to_data: dict[tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords)
        lct = Context(q2i= q2i, i2q= i2q, stab_to_data = stab_to_data)

    ##########################################################
    # Adding Indexes and shared information into lct dataclass
    ##########################################################

    patches : dict[str, Patch] = {"patch": Patch.from_coords(qubit_coords, q2i),}

    ###################################
    # 5. Building Initilization Circuit
    ###################################

    # Check whether we need Y basis initilization
    if is_y:
        initial_circuit = y_initial(lct = lct, patches = patches, cfg = cfg, noise = noise)

    else:
        initial_circuit = initial(lct = lct, patches = patches, cfg = cfg, noise = noise)

    ################################
    # 6. Building repetition Circuit
    ################################

    if is_y:
        repeat_circ = y_repetition_circ(lct = lct, patches = patches, cfg = cfg, noise = noise)
        
        switch_circ = y_switch_circ(lct = lct, patches = patches, cfg = cfg, noise = noise)

        initial_circuit += repeat_circ
        initial_circuit += switch_circ       
        
    else:
        repeat_circ = repetition_circ(lct = lct, patches = patches, cfg = cfg, noise = noise)
        
        initial_circuit += repeat_circ

    #############################################################
    # Implement additional Circuit if Logical H gate was selected
    #############################################################

    """
    We now rund d rounds with flipped stabilizer roles
    -> i.e. X stabilizers convert to z stabilizers and x to z
    -> Flipped the stabs_to_data formalism and changed inside the function the role of x and z stab indices
    """

    #Determine if flip needed by helper
    flip_needed = _needs_flip(state_init= state_init,
                              log_obs= log_obs,
                              logical_h= logical_h)

    # Adding the needed circuits
    if flip_needed is True: 

        repeat_switch_init = h_switched_circ_init(lct = lct, patches = patches, cfg = cfg, noise = noise)

        repeat_switched = h_switched_circ(lct = lct, patches = patches, cfg = cfg, noise = noise)
        
        initial_circuit += repeat_switch_init
        initial_circuit += repeat_switched

    ###############################
    # 11. Adding State initiliztion
    ###############################

    state_init_circuit = reset(lct = lct, patches = patches, cfg = cfg, logical_h = flip_needed)

    if not is_y:
        final_measurement = final_m(lct = lct, patches = patches, cfg = cfg, noise = noise, is_flipped = flip_needed)
        state_init_circuit += initial_circuit
        state_init_circuit += final_measurement

    ##################################
    # Adding Actual Y basis Memory run
    ##################################

    else:
        y_memory = y_repetition_circ(lct = lct, patches = patches, cfg = cfg, noise = noise, memory_round= True)

        state_init_circuit += initial_circuit
        state_init_circuit += y_memory

        # Adding the basis reverse
        reverse_switch = y_rev_switch_circ(lct = lct, patches = patches, cfg = cfg, noise = noise)
        
        #############################
        # Adding logical Y Observable
        #############################

        if log_obs == "Y":

            logical_circ_contraction = reverse_switch
            logical_circ_creation = switch_circ

            logical_x_string = []
            logical_z_string = []

            fixed_coord = (distance * 2) - 1

            # Finding logical Strings for x and z
            for imag in range(1, (distance * 2) - 1, 2):
                logical_x_string.append(q2i[fixed_coord + 1j * imag])

            for real in range(1, (distance * 2) - 1, 2):
                logical_z_string.append(q2i[real + fixed_coord * 1j])


            # Adding logical z string
            logical_xyz_string = '*'.join([f"Z{idz}" for j, idz in enumerate(logical_z_string)] + 
                                    [f"Y{q2i[fixed_coord + fixed_coord * 1j]}"] + 
                                    [f"X{idx}" for j, idx in enumerate(logical_x_string)])
            
            logical_creation = f"{1} -> {logical_xyz_string}"
            logical_contraction = f"{logical_xyz_string} -> {1}"

            (logical_creation_rec,) = logical_circ_creation.circuit.solve_flow_measurements([stim.Flow(logical_creation)])
            (logical_contraction_rec,) = logical_circ_contraction.circuit.solve_flow_measurements([stim.Flow(logical_contraction)])

            # Adding the final Measurement Round & Missing Detectors
            state_init_circuit += reverse_switch

            # Calculating target rec pos
            rec_pos = []

            contraction_records = reverse_switch.circuit.num_measurements
            creation_records= reverse_switch.circuit.num_measurements + y_memory.circuit.num_measurements + switch_circ.circuit.num_measurements

            # Adding the Observable
            for index_creation in logical_creation_rec:
                current_rec_crea = creation_records - index_creation
                rec_pos.append(- current_rec_crea)
            
            for index_contraction in logical_contraction_rec:
                current_rec_cont = contraction_records - index_contraction
                rec_pos.append(- current_rec_cont)

            state_init_circuit.circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(k) for k in rec_pos], 0)

            #############################################################################
            # Y ONLY: Adding needed y_inital rounds in order top guarentee faul tolerance
            #############################################################################

            final_measurement = y_repetition_circ(lct = lct, patches = patches, cfg = cfg, noise = noise, ft_round = True)
            state_init_circuit += final_measurement

        # Calc the measurement rec pos for the X/Z Basis measruement
        if log_obs in {"X", "Z"}:

            """
            Adding ft round can be skipped as state was already measured (This is not FT either way)
            """

            # Adding Empty Meas Circ
            final_measurement = CircuitResult(circuit= stim.Circuit())

            # Finding full target history within full circuit
            recs_after = reverse_switch.circuit.num_measurements
            obs_index = [i - recs_after for i in y_memory.obs_indices]

            final_measurement.obs_indices = obs_index

    ############################################################
    # Return Circuit and measurement rec postitions for logicals
    ############################################################

    """
    If another Basis for measurement then init was chosen one needs the 
    measurement observable inices in order to determine the current measurement 
    result per sample
    """

    if state_init in {"+i", "-i"}:

        if final_measurement.obs_indices is not None:

            annotate_circuit = state_init_circuit.circuit

            # Add all Detectors with package
            circ_with_dets = annotate_detectors_automatically(annotate_circuit)

            return circ_with_dets , final_measurement.obs_indices

        else:   

            annotate_circuit = state_init_circuit.circuit

            # Add all Detectors with package
            circ_with_dets = annotate_detectors_automatically(annotate_circuit)

            return circ_with_dets
        
    else:

        if final_measurement.obs_indices is not None:

            return state_init_circuit.circuit, final_measurement.obs_indices

        else:   

            return state_init_circuit.circuit


