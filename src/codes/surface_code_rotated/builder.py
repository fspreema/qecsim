import stim
from tqecd import annotate_detectors_automatically

from src.codes.surface_code_rotated.get_stab_pairings import SurfacePairings
from src.codes.surface_code_rotated.surface_geom import SurfaceGeometry
from src.core.data_models import (
    CircuitResult,
    ConfigSurface as Config,
    Context,
    NoiseModel,
    Patch,
)
from src.core.noise_models import CircuitNoise

from .circuits.final_measure import final_m
from .circuits.h_switched_init import h_switched_circ_init
from .circuits.h_switched_round import h_switched_circ
from .circuits.initial import initial
from .circuits.repetition import repetition_circ
from .circuits.reset import reset
from .circuits.y_initial import y_initial
from .circuits.y_repetition import y_repetition_circ
from .circuits.y_rev_switch import y_rev_switch_circ
from .circuits.y_switch import y_switch_circ

Coord = complex
Label = str
Index = int
Pair = tuple[Coord, Coord]

__all__ = ["rotated_surface_code"]

# -------------------------
# Helper Functions
# -------------------------


def _needs_flip(state_init: str, log_obs: str, logical_h: bool) -> bool:
    # Determining if flip is needed
    return (state_init in {"0", "1"} and log_obs == "X" and logical_h) or (
        state_init in {"+", "-"} and log_obs == "Z" and logical_h
    )


# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------


def rotated_surface_code(
    distance: int,
    rounds: int,
    *,
    state_init: str,
    log_obs: str,
    logical_h: bool = False,
    noise_depol_data_init: float = 0.0,
    noise_measure_flip: float = 0.0,
    noise_after_reset: float = 0.0,
    noise_after_clifford_depol: float = 0.0,
) -> stim.Circuit:
    """
    Generates Rotated-Surface-Code

    Args:
        Distance (int): Distance of the surface code i.e. lattice size
        Rounds (int): Number of syndrome measurement rounds per Shot
        noise (float): Probability for x y and z error in Pauli-Channel

    Information:
        Logical Operator is Z Operator and pre-Defined!
        Logical State: 0 -> Only z Stabilizer detectors in the first round as x detectors
                            are non deterministc for the first round
                            (STILL: COMPLETE MEASUREMENT)
                         -> In theory we don not even need to meassure the X stabilizers
                            at all because we do not have phase errors
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
    cfg = Config(distance=distance, state_init=state_init, obs=log_obs, rounds=rounds)

    noise = NoiseModel(
        before_round_depol=noise_depol_data_init,
        before_m_flip_prob=noise_measure_flip,
        after_r_flip=noise_after_reset,
        after_c_depol_prob=noise_after_clifford_depol,
    )

    ###############################################################
    # 1. Build independent square patches (using geometry function)
    ###############################################################

    """
    As we need one x and one z edge for y init, we need to differentiate the boundray labels
    """

    # Is logical Basis y?
    is_y = state_init in {"+i", "-i"}

    if is_y:
        geometry = SurfaceGeometry(distance=distance, y_basis=True, starting_stabilizer_x=False)
        qubit_coords: dict[Coord, Label] = geometry.coords
    else:
        geometry = SurfaceGeometry(distance=distance, y_basis=False, starting_stabilizer_x=True)
        qubit_coords: dict[Coord, Label] = geometry.coords

    ###############################################
    # 3. Indexing All Qubits From given Coordinates
    ###############################################

    # Indexing Qubits
    q2i = geometry._get_q2i()

    # Reverse Indexing
    i2q = geometry._get_i2q()

    ############################################################
    # 4. Adding the Mapping from Stabilizer to Data for later CX
    # gate implementation and add into dataclass
    ############################################################

    """
    stab_to_data_switch: For the implementation fo the XCY gates for the Y basis init
    stab_to-data_flipped: For the logical H gate implementation -> Switch of X and Z stabilizers
    """

    if is_y:  # Logical Y Basis
        surface_pairings1 = SurfacePairings(
            patch=qubit_coords,
            distance=distance,
            y_basis=True,
        )
        stab_to_data = surface_pairings1.get_schedule()

        surface_pairings2 = SurfacePairings(
            patch=qubit_coords,
            distance=distance,
            y_basis=True,
            y_switch=True,
        )
        stab_to_data_switch, stab_to_data_xcy = surface_pairings2.get_schedule()

        surface_pairings3 = SurfacePairings(
            patch=qubit_coords,
            distance=distance,
            y_basis=True,
            y_memory=True,
        )
        stab_to_data_memory: dict[tuple[Coord, Coord], str] = surface_pairings3.get_schedule()

        lct = Context(
            q2i=q2i,
            i2q=i2q,
            stab_to_data=stab_to_data,
            stab_to_data_modified=stab_to_data_switch,
            stab_to_data_modified2=stab_to_data_xcy,
            stab_to_data_modified3=stab_to_data_memory,
        )

    elif logical_h:  # Logical H Gate
        surface_pairings = SurfacePairings(
            patch=qubit_coords,
            distance=distance,
        )

        stab_to_data: dict[tuple[Coord, Coord], str] = surface_pairings.get_schedule()

        surface_pairings_flipped = SurfacePairings(
            patch=qubit_coords,
            distance=distance,
            is_flipped=True,
        )

        stab_to_data_flipped: dict[tuple[Coord, Coord], str] = (
            surface_pairings_flipped.get_schedule()
        )

        lct = Context(
            q2i=q2i,
            i2q=i2q,
            stab_to_data=stab_to_data,
            stab_to_data_modified=stab_to_data_flipped,
        )
    else:  # Regular X or Z Basis
        surface_pairings_reg = SurfacePairings(
            patch=qubit_coords,
            distance=distance,
        )

        stab_to_data: dict[tuple[Coord, Coord], str] = surface_pairings_reg.get_schedule()
        lct = Context(q2i=q2i, i2q=i2q, stab_to_data=stab_to_data)

    ##########################################################
    # Adding Indexes and shared information into lct dataclass
    ##########################################################

    patches: dict[str, Patch] = {"patch": Patch.from_coords(qubit_coords, q2i)}

    ###################################
    # 5. Building Initilization Circuit
    ###################################

    # Check whether we need Y basis initilization
    if is_y:
        initial_circuit = y_initial(lct=lct, patch=patches["patch"])

    else:
        initial_circuit = initial(lct=lct, patches=patches, cfg=cfg)

    ################################
    # 6. Building repetition Circuit
    ################################

    if is_y:
        repeat_circ = y_repetition_circ(lct=lct, patch=patches["patch"], cfg=cfg, noise=noise)

        switch_circ = y_switch_circ(lct=lct, patch=patches["patch"], cfg=cfg)

        initial_circuit += repeat_circ
        initial_circuit += switch_circ

    else:
        repeat_circ = repetition_circ(lct=lct, patches=patches, cfg=cfg, noise=noise)

        initial_circuit += repeat_circ

    #############################################################
    # Implement additional Circuit if Logical H gate was selected
    #############################################################

    """
    We now rund d rounds with flipped stabilizer roles
    -> i.e. X stabilizers convert to z stabilizers and x to z
    -> Flipped the stabs_to_data formalism and changed inside 
       the function the role of x and z stab indices
    """

    # Determine if flip needed by helper
    flip_needed = _needs_flip(state_init=state_init, log_obs=log_obs, logical_h=logical_h)

    # Adding the needed circuits
    if flip_needed is True:
        repeat_switch_init = h_switched_circ_init(lct=lct, patches=patches, noise=noise)

        repeat_switched = h_switched_circ(lct=lct, patches=patches, cfg=cfg)

        initial_circuit += repeat_switch_init
        initial_circuit += repeat_switched

    ###############################
    # 11. Adding State initiliztion
    ###############################

    state_init_circuit = reset(lct=lct, patches=patches, cfg=cfg, logical_h=flip_needed)

    if not is_y:
        final_measurement = final_m(
            lct=lct,
            patches=patches,
            cfg=cfg,
            is_flipped=flip_needed,
        )

        state_init_circuit += initial_circuit
        state_init_circuit += final_measurement

    ##################################
    # Adding Actual Y basis Memory run
    ##################################

    else:
        y_memory = y_repetition_circ(
            lct=lct,
            patch=patches["patch"],
            cfg=cfg,
            noise=noise,
            memory_round=True,
        )

        state_init_circuit += initial_circuit
        state_init_circuit += y_memory

        # Adding the basis reverse
        reverse_switch = y_rev_switch_circ(
            lct=lct,
            patches=patches,
        )

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
            logical_xyz_string = "*".join(
                [f"Z{idz}" for j, idz in enumerate(logical_z_string)]
                + [f"Y{q2i[fixed_coord + fixed_coord * 1j]}"]
                + [f"X{idx}" for j, idx in enumerate(logical_x_string)],
            )

            logical_creation = f"{1} -> {logical_xyz_string}"
            logical_contraction = f"{logical_xyz_string} -> {1}"

            (logical_creation_rec,) = logical_circ_creation.circuit.solve_flow_measurements(
                [stim.Flow(logical_creation)],
            )
            (logical_contraction_rec,) = logical_circ_contraction.circuit.solve_flow_measurements(
                [stim.Flow(logical_contraction)],
            )

            # Adding the final Measurement Round & Missing Detectors
            state_init_circuit += reverse_switch

            # Calculating target rec pos
            rec_pos = []

            contraction_records = reverse_switch.circuit.num_measurements
            creation_records = (
                reverse_switch.circuit.num_measurements
                + y_memory.circuit.num_measurements
                + switch_circ.circuit.num_measurements
            )

            # Adding the Observable
            for index_creation in logical_creation_rec:
                current_rec_crea = creation_records - index_creation
                rec_pos.append(-current_rec_crea)

            for index_contraction in logical_contraction_rec:
                current_rec_cont = contraction_records - index_contraction
                rec_pos.append(-current_rec_cont)

            state_init_circuit.circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_rec(k) for k in rec_pos],
                0,
            )

            #############################################################################
            # Y ONLY: Adding needed y_inital rounds in order top guarentee faul tolerance
            #############################################################################

            final_measurement = y_repetition_circ(
                lct=lct,
                patch=patches["patch"],
                cfg=cfg,
                noise=noise,
            )

            state_init_circuit += final_measurement

        # Calc the measurement rec pos for the X/Z Basis measruement
        if log_obs in {"X", "Z"}:
            """
            Adding ft round can be skipped as state was already measured (This is not FT either way)
            """

            # Adding Empty Meas Circ
            final_measurement = CircuitResult(circuit=stim.Circuit())

            # Finding full target history within full circuit
            recs_after = reverse_switch.circuit.num_measurements
            obs_index = [i - recs_after for i in y_memory.obs_indices]

            final_measurement.obs_indices = obs_index

    ###############################
    # Adding Noise to final Circuit
    ###############################

    if (
        noise.before_m_flip_prob > 0.0
        or noise.after_r_flip > 0.0
        or noise.after_c_depol_prob > 0.0
        or noise.before_round_depol > 0.0
    ):
        # Create noise dict
        noise_dict = {
            "before_round_depol": noise.before_round_depol,
            "before_m_flip_prob": noise.before_m_flip_prob,
            "after_r_flip": noise.after_r_flip,
            "after_c_depol_prob": noise.after_c_depol_prob,
        }

        circuit_noise_builder = CircuitNoise(circuit=state_init_circuit.circuit, noise=noise_dict)

        # Building final circuit with noise
        return_circuit = circuit_noise_builder.apply()

    else:
        return_circuit = state_init_circuit.circuit

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
            # Add all Detectors with package
            circ_with_dets = annotate_detectors_automatically(return_circuit)

            return circ_with_dets, final_measurement.obs_indices

        else:
            # Add all Detectors with package
            circ_with_dets = annotate_detectors_automatically(return_circuit)

            return circ_with_dets

    else:
        if final_measurement.obs_indices is not None:
            return return_circuit, final_measurement.obs_indices

        else:
            return return_circuit
