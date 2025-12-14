import stim

from src.codes.lattice_surgery.get_stab_pairings import LatticeSurgeryPairings
from src.codes.lattice_surgery.logical_strings import get_logical_strings
from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.core.data_models import (
    ConfigLatticeSurgery as Config,
    LatticeContext,
    NoiseModel,
    PatchAncilla,
    PatchControl,
    PatchSurgery,
    PatchTarget,
)
from src.core.noise_models import CircuitNoise

from .circuits.final_measure import final_m
from .circuits.initial import initial
from .circuits.merge import merge
from .circuits.reset import reset
from .circuits.split import split

Coord = complex
Label = str
Index = int
Pair = tuple[Coord, Coord]

__all__ = ["surgery_circuit"]

# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------


def surgery_circuit(
    distance: int,
    *,
    round_num: int = 0,
    round_merge: int = 0,
    round_split: int = 0,
    target_state_init: str,
    control_state_init: str,
    flow_observable: str,
    noise_depol_data_init: float = 0.0,
    noise_measure_flip: float = 0.0,
    noise_after_reset: float = 0.0,
    noise_after_clifford_depol: float = 0.0,
) -> stim.Circuit:
    """
    Returns the full lattice surgery circuit

    Arguments:
        -> target_state_init:
        In which basis should the target lattice be initlized?
        -> control_state_init:
        In which basis should the control lattice be initilized?
        -> flow_observable:
        What Observable should stim track?
        -> distance/rounds:
        Distance and rounds of the circuit (If not given any value distance = rounds)
        -> noise_depol_data_init:
        Depolaization noise after initlization of the data qubits with given probability
        -> noise_measure_flip:
        X_Error before measurement to simulate faulty measurement with given probability
        -> noise_after_reset:
        X_Error after reset of ancilla in order so simulate faulty reset with given probability
        -> noise_after_clifford_depol:
        Depolarization gate after each Clifford (CX and H) with given probability

    Returns:
        -> Fully implemented CX-Gate in stim.Circuit format
    """
    ####################
    # Check Round number
    ####################

    if round_num == 0:
        round_num = distance

    if round_merge == 0:
        round_merge = distance

    if round_split == 0:
        round_split = distance

    #########################################
    # Input fixed run settings into dataclass
    #########################################

    cfg = Config(
        distance=distance,
        flow_observable=flow_observable,
        target_state_init=target_state_init,
        control_state_init=control_state_init,
    )

    noise = NoiseModel(
        before_round_depol=noise_depol_data_init,
        before_m_flip_prob=noise_measure_flip,
        after_r_flip=noise_after_reset,
        after_c_depol_prob=noise_after_clifford_depol,
    )

    #######################################################################################
    # 1. Build independent square patches (using geometry function) and add Boundary Labels
    #######################################################################################
    full_geometry = SurgeryGeometry(distance=distance)

    qubit_coords_ancilla: dict[Coord, Label] = full_geometry.coords_ancilla
    qubit_coords_target: dict[Coord, Label] = full_geometry.coords_target
    qubit_coords_control: dict[Coord, Label] = full_geometry.coords_control
    qubit_coords_surgery: dict[Coord, Label] = full_geometry.coords_surgery

    # Merge into different Patches
    """
    Different Patches are needed, because of different Keywords on 
    identical Coordinates (inside dict.):

    -> X-Stab-Boundary-Above-Control & X-Stab-Boundary-Below-Ancilla f.ex. 
       get Keywords for surgery stabilizers
    """

    ################################################################################
    # 3. Adding the Mapping from Stabilizer to Data for later CX gate implementation
    ################################################################################

    pairings_ancilla = LatticeSurgeryPairings(
        qubit_coords=qubit_coords_ancilla,
        merging=False,
    )
    stab_to_data_ancilla: dict[tuple[Coord, Coord], str] = pairings_ancilla.get_schedule()

    pairings_target = LatticeSurgeryPairings(
        qubit_coords=qubit_coords_target,
        merging=False,
    )
    stab_to_data_target: dict[tuple[Coord, Coord], str] = pairings_target.get_schedule()

    pairings_control = LatticeSurgeryPairings(
        qubit_coords=qubit_coords_control,
        merging=False,
    )
    stab_to_data_control: dict[tuple[Coord, Coord], str] = pairings_control.get_schedule()

    pairings_surgery_ac = LatticeSurgeryPairings(
        qubit_coords=qubit_coords_surgery,
        merging=True,
        merging_type="AC",
    )
    stab_to_data_surgery_ac: dict[tuple[Coord, Coord], str] = pairings_surgery_ac.get_schedule()

    pairings_surgery_at = LatticeSurgeryPairings(
        qubit_coords=qubit_coords_surgery,
        merging=True,
        merging_type="AT",
    )
    stab_to_data_surgery_at: dict[tuple[Coord, Coord], str] = pairings_surgery_at.get_schedule()

    ###############################################
    # 4. Indexing All Qubits From given Coordinates
    ###############################################

    # Indexing Qubits
    q2i: dict[complex, int] = full_geometry._get_q2i()

    # Reverse Indexing
    i2q: dict[int, complex] = full_geometry._get_i2q()

    ##########################################################
    # Adding Indexes and shared information into lct dataclass
    ##########################################################

    lct = LatticeContext(
        q2i=q2i,
        i2q=i2q,
        stab_to_data=stab_to_data_ancilla | stab_to_data_control | stab_to_data_target,
        stab_to_data_surgery_ac=stab_to_data_surgery_ac,
        stab_to_data_surgery_at=stab_to_data_surgery_at,
        surgery_coords=qubit_coords_surgery,
    )

    patches: dict[str, PatchAncilla | PatchTarget | PatchControl | PatchSurgery] = {
        "ancilla": PatchAncilla.from_coords(qubit_coords_ancilla, q2i),
        "target": PatchTarget.from_coords(qubit_coords_target, q2i),
        "control": PatchControl.from_coords(qubit_coords_control, q2i),
        "surgery": PatchSurgery.from_coords(qubit_coords_surgery, q2i),
    }

    #################################
    # Add Y Init FT Proc. if selected
    #################################

    if target_state_init in {"Y+", "Y-"} or control_state_init in {"Y+", "Y-"}:
        reset_circ, y_circ = reset(
            lct=lct,
            patches=patches,
            cfg=cfg,
            noise=noise,
        )

    else:
        reset_circ = reset(
            lct=lct,
            patches=patches,
            cfg=cfg,
            noise=noise,
        )

    ###################################
    # 5. Building Initilization Circuit
    ###################################

    if target_state_init in {"Y+", "Y-"} or control_state_init in {"Y+", "Y-"}:
        # Adding y basis initilization circuit before normal initilization
        flow_circuit = stim.Circuit()
        flow_circuit += y_circ

        flow_circuit += initial(
            lct=lct,
            patches=patches,
            cfg=cfg,
            noise=noise,
        )

    else:
        flow_circuit = stim.Circuit()

        flow_circuit += initial(
            lct=lct,
            patches=patches,
            cfg=cfg,
            noise=noise,
        )

    #############################################
    # 6. Building Merging Ancilla Control Circuit
    #############################################

    """
    Depending on the proposed flows either standard XX or ZZ parity measurements
    are used, or new YY measurements which need additional s gates
    -> Therefore different merging or splitting procedures for certain flows
    """

    merged_circuit_ac = merge(
        lct=lct,
        patches=patches,
        cfg=cfg,
        merging_type="AC",
    )

    ###############################################
    # 5. Building Splitting Ancilla Control Circuit
    ###############################################

    split_circuit_ac = split(
        lct=lct,
        patches=patches,
        cfg=cfg,
        split_type="AC",
    )

    ############################################
    # 6. Building Merging Ancilla Target Circuit
    ############################################

    merged_circuit_at = merge(
        lct=lct,
        patches=patches,
        cfg=cfg,
        merging_type="AT",
    )

    ##############################################
    # 7. Building splitting Ancilla Target Circuit
    ##############################################

    split_circuit_at = split(
        lct=lct,
        patches=patches,
        cfg=cfg,
        split_type="AT",
    )

    ################################################################################
    # 8. Creating Clipped Circuit (Without State intilization and final measurement)
    ################################################################################

    flow_circuit += merged_circuit_ac
    flow_circuit += split_circuit_ac
    flow_circuit += merged_circuit_at
    flow_circuit += split_circuit_at

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

    # Y logical components (using shared helper for correctness)
    log_strings = get_logical_strings(q2i, distance)
    # Shifted variants for Y-including flows (changes c_x at x=2d-1, t_z at y=2d-1)
    log_strings_shift = get_logical_strings(
        q2i,
        distance,
        shift_cx_for_y=True,
        shift_tz_for_y=True,
        shift_tx_for_y=True,
        shift_cz_for_y=True,
    )
    c_y = log_strings["c_y"]  # dict with keys: z_string, y_corner, x_string
    t_y = log_strings["t_y"]

    ##################################
    # 10. Building logical Observables
    ##################################

    # -----------------ALL-X-----------------------

    if flow_observable == "X -> XX":
        if control_state_init in {"X+", "X-"}:
            if target_state_init in {"X+", "X-"}:
                left = "*".join(f"X{i}" for i in c_log_x)
                right = "*".join(f"X{i}" for i in c_log_x + t_log_x)
                result = f"{left} -> {right}"

                (included_measurements,) = flow_circuit.solve_flow_measurements(
                    [
                        stim.Flow(result),
                    ],
                )

            else:
                raise ValueError("Wrong target basis for selected flow")

        else:
            raise ValueError("Wrong control basis for selected flow")

    elif flow_observable == "XX -> X":
        if control_state_init in {"X+", "X-"}:
            if target_state_init in {"X+", "X-"}:
                left = "*".join(f"X{i}" for i in c_log_x + t_log_x)
                right = "*".join(f"X{i}" for i in c_log_x)
                result = f"{left} -> {right}"

                (included_measurements,) = flow_circuit.solve_flow_measurements(
                    [
                        stim.Flow(result),
                    ],
                )

            else:
                raise ValueError("Invalid target state")

        else:
            raise ValueError("Wrong control basis for selected flow")

    elif flow_observable == "X -> X":
        if control_state_init in {"X+", "X-", "Z0", "Z1"}:
            if target_state_init in {"X+", "X-"}:
                left = "*".join(f"X{i}" for i in t_log_x)
                right = "*".join(f"X{i}" for i in t_log_x)
                result = f"{left} -> {right}"

                (included_measurements,) = flow_circuit.solve_flow_measurements(
                    [
                        stim.Flow(result),
                    ],
                )

            else:
                raise ValueError("Invalid target basis for selected flow")

        else:
            raise ValueError("Invalid control state")

    # -------------ALL-Z--------------------------

    elif flow_observable == "Z -> ZZ":
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1"}:
                left = "*".join(f"Z{i}" for i in t_log_z)
                right = "*".join(f"Z{i}" for i in c_log_z + t_log_z)
                result = f"{left} -> {right}"

                (included_measurements,) = flow_circuit.solve_flow_measurements(
                    [
                        stim.Flow(result),
                    ],
                )

            else:
                raise ValueError("Wrong target basis for selected flow")

        else:
            raise ValueError("Wrong control basis for selected flow")

    elif flow_observable == "ZZ -> Z":
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1"}:
                left = "*".join(f"Z{i}" for i in c_log_z + t_log_z)
                right = "*".join(f"Z{i}" for i in t_log_z)
                result = f"{left} -> {right}"

                (included_measurements,) = flow_circuit.solve_flow_measurements(
                    [
                        stim.Flow(result),
                    ],
                )

            else:
                raise ValueError("Wrong target basis for selected flow")

        else:
            raise ValueError("Wrong control basis for selected flow")

    elif flow_observable == "Z -> Z":
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1", "X+", "X-"}:
                left = "*".join(f"Z{i}" for i in c_log_z)
                right = "*".join(f"Z{i}" for i in c_log_z)
                result = f"{left} -> {right}"

                (included_measurements,) = flow_circuit.solve_flow_measurements(
                    [
                        stim.Flow(result),
                    ],
                )

            else:
                raise ValueError("Invalid target state")

        else:
            raise ValueError("Wrong control basis for selected flow")

    # -------------------XZ-MIX-------------------------------------

    elif flow_observable == "ZX -> ZX":
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"X+", "X-"}:
                left = "*".join([f"Z{i}" for i in c_log_z] + [f"X{i}" for i in t_log_x])
                right = "*".join([f"Z{i}" for i in c_log_z] + [f"X{i}" for i in t_log_x])
                result = f"{left} -> {right}"

                (included_measurements,) = flow_circuit.solve_flow_measurements(
                    [
                        stim.Flow(result),
                    ],
                )

            else:
                raise ValueError("Invalid target state")

        else:
            raise ValueError("Wrong control basis for selected flow")

    # -----------------------Y-Included-measurements----------------

    # As the y observable is initlizized inside the circuit we do not propose a flow with Y
    # -> So YZ -> XY is proposed as Z -> XY as Y is created iniside the circuit
    # by the y basis initilization procedure

    elif flow_observable == "YZ -> XY":
        # Left: control Y, target Z; Right: control X, target Y
        if control_state_init in {"Y+", "Y-"} and target_state_init in {"Z0", "Z1"}:
            left_terms = [f"Z{i}" for i in log_strings_shift["t_z"]]
            right_terms = (
                [f"X{i}" for i in log_strings_shift["c_x"]]
                + [f"Z{i}" for i in t_y["z_string"]]
                + [f"Y{i}" for i in t_y["y_corner"]]
                + [f"X{i}" for i in t_y["x_string"]]
            )
            result = f"{'*'.join(left_terms)} -> {'*'.join(right_terms)}"

            (included_measurements,) = flow_circuit.solve_flow_measurements([stim.Flow(result)])

        else:
            raise ValueError("Invalid basis for selected flow: YZ -> XY")

    elif flow_observable == "YX -> YI":
        # Left: control Y, target X; Right: control Y, target I
        if control_state_init in {"Y+", "Y-"} and target_state_init in {"X+", "X-"}:
            left_terms = [f"X{i}" for i in log_strings_shift["t_x"]]
            right_terms = (
                [f"Z{i}" for i in c_y["z_string"]]
                + [f"Y{i}" for i in c_y["y_corner"]]
                + [f"X{i}" for i in c_y["x_string"]]
            )
            result = f"{'*'.join(left_terms)} -> {'*'.join(right_terms)}"

            (included_measurements,) = flow_circuit.solve_flow_measurements([stim.Flow(result)])
        else:
            raise ValueError("Invalid basis for selected flow: YX -> YI")

    elif flow_observable == "YI -> YX":
        # Left: control Y, target I; Right: control Y, target X
        if control_state_init in {"Y+", "Y-"} and target_state_init in {"Z0", "Z1", "X+", "X-"}:
            right_terms = (
                [f"X{i}" for i in log_strings_shift["t_x"]]
                + [f"Z{i}" for i in c_y["z_string"]]
                + [f"Y{i}" for i in c_y["y_corner"]]
                + [f"X{i}" for i in c_y["x_string"]]
            )
            result = f"{1} -> {'*'.join(right_terms)}"

            (included_measurements,) = flow_circuit.solve_flow_measurements([stim.Flow(result)])
        else:
            raise ValueError("Invalid basis for selected flow: YI -> YX")

    elif flow_observable == "YY -> XZ":
        # Left: control Y, target Y; Right: control X, target Z
        if control_state_init in {"Y+", "Y-"} and target_state_init in {"Y+", "Y-"}:
            right_terms = [f"X{i}" for i in log_strings_shift["c_x"]] + [
                f"Z{i}" for i in log_strings_shift["t_z"]
            ]
            result = f"{1} -> {'*'.join(right_terms)}"

            (included_measurements,) = flow_circuit.solve_flow_measurements([stim.Flow(result)])
        else:
            raise ValueError("Invalid basis for selected flow: YY -> XZ")

    elif flow_observable == "IY -> ZY":
        # Left: control I, target Y; Right: control Z, target Y
        if target_state_init in {"Y+", "Y-"}:
            right_terms = (
                [f"Z{i}" for i in log_strings_shift["c_z"]]
                + [f"Z{i}" for i in t_y["z_string"]]
                + [f"Y{i}" for i in t_y["y_corner"]]
                + [f"X{i}" for i in t_y["x_string"]]
            )
            result = f"{1} -> {'*'.join(right_terms)}"

            (included_measurements,) = flow_circuit.solve_flow_measurements([stim.Flow(result)])
        else:
            raise ValueError("Invalid basis for selected flow: IY -> ZY")

    elif flow_observable == "XY -> YZ":
        # Left: control X, target Y; Right: control Y, target Z
        if control_state_init in {"X+", "X-"} and target_state_init in {"Y+", "Y-"}:
            left_terms = [f"X{i}" for i in log_strings_shift["c_x"]]
            right_terms = (
                [f"Z{i}" for i in c_y["z_string"]]
                + [f"Y{i}" for i in c_y["y_corner"]]
                + [f"X{i}" for i in c_y["x_string"]]
                + [f"Z{i}" for i in log_strings_shift["t_z"]]
            )
            result = f"{'*'.join(left_terms)} -> {'*'.join(right_terms)}"

            (included_measurements,) = flow_circuit.solve_flow_measurements([stim.Flow(result)])
        else:
            raise ValueError("Invalid basis for selected flow: XY -> YZ")

    else:
        raise ValueError("Invalid Flow selected")

    #################################
    # 11. Adding State initialization
    #################################

    final_measurement = final_m(
        lct=lct,
        patches=patches,
        cfg=cfg,
        flow=flow_observable,
    )

    reset_circ += flow_circuit

    ##################################################
    # 12. Adding logical Observable given by stim.Flow
    ##################################################

    # Calculating target rec pos
    rec_pos = []

    for index in included_measurements:
        current_rec_tar = flow_circuit.num_measurements - index
        rec_pos.append(-current_rec_tar)

    reset_circ.append("OBSERVABLE_INCLUDE", [stim.target_rec(k) for k in rec_pos], 0)

    ##########################
    # Adding final measurement
    ##########################

    reset_circ += final_measurement

    ##########################################
    # Adding noise depending on noise selected
    ##########################################

    # 1) Circuit Noise Model
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

        circuit_noise_builder = CircuitNoise(circuit=reset_circ, noise=noise_dict)

        # Building final circuit with noise
        return_circuit = circuit_noise_builder.apply()

    else:
        return_circuit = reset_circ

    # Normalize to avoid double TICKs after composing subcircuits
    return return_circuit
