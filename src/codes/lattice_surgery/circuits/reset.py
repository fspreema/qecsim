from collections.abc import Mapping

import stim

from src.codes.lattice_surgery.logical_strings import get_logical_strings
from src.codes.surface_code_rotated.circuits.reset import reset as y_reset

# Import Part of Y Surface Code Init
from src.codes.surface_code_rotated.circuits.y_initial import y_initial
from src.codes.surface_code_rotated.circuits.y_repetition import y_repetition_circ
from src.codes.surface_code_rotated.circuits.y_switch import y_switch_circ
from src.codes.surface_code_rotated.get_stab_pairings import SurfacePairings
from src.codes.surface_code_rotated.surface_geom import SurfaceGeometry
from src.core.data_models import (
    ConfigLatticeSurgery as Config,
    ConfigSurface as SurfaceConfig,
    Context,
    LatticeContext,
    NoiseModel,
    Patch,
    PatchAncilla,
    PatchControl,
    PatchSurgery,
    PatchTarget,
)

Coord = complex
Label = str
Index = int
Pair = tuple[Coord, Coord]

__all__ = ["reset"]


def reset(
    *,
    lct: LatticeContext,
    patches: Mapping[str, PatchAncilla | PatchTarget | PatchControl | PatchSurgery],
    cfg: Config,
    noise: NoiseModel,
) -> stim.Circuit:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # -Loading in Patches
    ancilla_patch = patches["ancilla"]
    target_patch = patches["target"]
    control_patch = patches["control"]

    # -Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    distance = cfg.distance
    control_state_init = cfg.control_state_init
    target_state_init = cfg.target_state_init
    stab_to_data = lct.stab_to_data

    # -Retrieving Data Coords
    data_control = control_patch.data
    data_target = target_patch.data

    # -Retrieving Index from Stabilizers of the Lattices
    x_stab_index_ancilla = ancilla_patch.x_stab
    z_stab_index_ancilla = ancilla_patch.z_stab
    x_stab_boundary_b_index_ancilla = ancilla_patch.x_bdy_b
    z_stab_boundary_r_index_ancilla = ancilla_patch.z_bdy_r
    x_stab_index_control = control_patch.x_stab
    z_stab_index_control = control_patch.z_stab
    x_stab_index_target = target_patch.x_stab
    z_stab_index_target = target_patch.z_stab

    # ----------------------------------------------
    # Creating List of all Stabilizers (No Double!)
    # ----------------------------------------------

    all_stabs_not_double = []

    # Setting double counter
    counter_x = 0
    counter_z = 0

    for index in (
        x_stab_index_ancilla
        + z_stab_index_ancilla
        + x_stab_index_control
        + z_stab_index_control
        + x_stab_index_target
        + z_stab_index_target
    ):
        # Double Values possible
        if index in x_stab_boundary_b_index_ancilla:
            # Value already appended?
            if counter_x == 0:
                all_stabs_not_double.append(index)
                counter_x += 1

        elif index in z_stab_boundary_r_index_ancilla:
            # Value already appended?
            if counter_z == 0:
                all_stabs_not_double.append(index)
                counter_z += 1

        else:
            all_stabs_not_double.append(index)

    # ------------------------------------------------------
    # Creating list of Logical X/Z string and their indices
    # ------------------------------------------------------
    """
    -> Used for swithcing of the state in a given basis
    """

    log_strings = get_logical_strings(q2i, distance)
    t_log_obs_z_index = log_strings["t_z"]
    t_log_obs_x_index = log_strings["t_x"]
    c_log_obs_z_index = log_strings["c_z"]
    c_log_obs_x_index = log_strings["c_x"]

    ########################
    # Define Initial Circuit
    ########################

    reset_circuit = stim.Circuit()

    # Appending Coords
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
        ("Z0", "Z1"): [
            ("R", data_control + data_target + all_stabs_not_double),
            ("X", t_log_obs_x_index),
        ],
        ("Z0", "X+"): [("RX", data_target), ("R", data_control + all_stabs_not_double)],
        ("Z0", "X-"): [
            ("RX", data_target),
            ("R", data_control + all_stabs_not_double),
            ("Z", t_log_obs_z_index),
        ],
        ("Z1", "Z0"): [
            ("R", data_control + data_target + all_stabs_not_double),
            ("X", c_log_obs_x_index),
        ],
        ("Z1", "Z1"): [
            ("R", data_control + data_target + all_stabs_not_double),
            ("X", c_log_obs_x_index + t_log_obs_x_index),
        ],
        ("Z1", "X+"): [
            ("RX", data_target),
            ("R", data_control + all_stabs_not_double),
            ("X", c_log_obs_x_index),
        ],
        ("Z1", "X-"): [
            ("RX", data_target),
            ("R", data_control + all_stabs_not_double),
            ("X", c_log_obs_x_index),
            ("Z", t_log_obs_z_index),
        ],
        ("X+", "Z0"): [("RX", data_control), ("R", data_target + all_stabs_not_double)],
        ("X+", "Z1"): [
            ("RX", data_control),
            ("R", data_target + all_stabs_not_double),
            ("X", t_log_obs_x_index),
            ("Z", c_log_obs_z_index),
        ],
        ("X+", "X+"): [("RX", data_control + data_target), ("R", all_stabs_not_double)],
        ("X+", "X-"): [
            ("RX", data_control + data_target),
            ("R", all_stabs_not_double),
            ("Z", t_log_obs_z_index),
        ],
        ("X-", "Z0"): [
            ("RX", data_control),
            ("R", data_target + all_stabs_not_double),
            ("Z", c_log_obs_z_index),
        ],
        ("X-", "Z1"): [
            ("RX", data_control),
            ("R", data_target + all_stabs_not_double),
            ("X", t_log_obs_x_index),
            ("Z", c_log_obs_z_index),
        ],
        ("X-", "X+"): [
            ("RX", data_control + data_target),
            ("R", all_stabs_not_double),
            ("Z", c_log_obs_z_index),
        ],
        ("X-", "X-"): [
            ("RX", data_control + data_target),
            ("R", all_stabs_not_double),
            ("Z", c_log_obs_z_index + t_log_obs_z_index),
        ],
    }

    y_patterns = {
        ("Y+", "Z0"): [("R", data_target + x_stab_index_target + z_stab_index_target)],
        ("Y+", "Z1"): [
            ("R", data_target + x_stab_index_target + z_stab_index_target),
            ("X", t_log_obs_x_index),
        ],
        ("Y+", "X+"): [("RX", data_target), ("R", x_stab_index_target + z_stab_index_target)],
        ("Y+", "X-"): [
            ("RX", data_target),
            ("R", x_stab_index_target + z_stab_index_target),
            ("Z", t_log_obs_z_index),
        ],
        ("Y-", "Z0"): [("R", data_target + x_stab_index_target + z_stab_index_target)],
        ("Y-", "Z1"): [
            ("R", data_target + x_stab_index_target + z_stab_index_target),
            ("X", t_log_obs_x_index),
        ],
        ("Y-", "X+"): [("RX", data_target), ("R", x_stab_index_target + z_stab_index_target)],
        ("Y-", "X-"): [
            ("RX", data_target),
            ("R", x_stab_index_target + z_stab_index_target),
            ("Z", t_log_obs_z_index),
        ],
        ("X+", "Y+"): [("RX", data_control), ("R", x_stab_index_control + z_stab_index_control)],
        ("X-", "Y+"): [
            ("RX", data_control),
            ("R", x_stab_index_control + z_stab_index_control),
            ("Z", c_log_obs_z_index),
        ],
        ("X+", "Y-"): [("RX", data_control), ("R", x_stab_index_control + z_stab_index_control)],
        ("X-", "Y-"): [
            ("RX", data_control),
            ("R", x_stab_index_control + z_stab_index_control),
            ("Z", c_log_obs_z_index),
        ],
        ("Z0", "Y+"): [("R", data_control + x_stab_index_control + z_stab_index_control)],
        ("Z1", "Y+"): [
            ("R", data_control + x_stab_index_control + z_stab_index_control),
            ("X", c_log_obs_x_index),
        ],
        ("Z0", "Y-"): [("R", data_control + x_stab_index_control + z_stab_index_control)],
        ("Z1", "Y-"): [
            ("R", data_control + x_stab_index_control + z_stab_index_control),
            ("X", c_log_obs_x_index),
        ],
    }

    # Apply the initialization pattern
    key = (control_state_init, target_state_init)

    if key in init_patterns:
        for gate, qubits in init_patterns[key]:
            reset_circuit.append(gate, qubits)

        reset_circuit.append("TICK")
        return reset_circuit

    elif key in y_patterns or key in {("Y+", "Y-"), ("Y-", "Y+"), ("Y+", "Y+"), ("Y-", "Y-")}:
        if key in y_patterns:
            for gate, qubits in y_patterns[key]:
                reset_circuit.append(gate, qubits)

        # Check if both control and target are Y-basis
        both_y = key[0] in {"Y+", "Y-"} and key[1] in {"Y+", "Y-"}

        if both_y:
            # Handle both control and target in Y basis
            # Build lattices for both control and target

            y_geometry_control = SurfaceGeometry(
                distance=distance,
                offset=0 + (distance * 2) * 1j,
                starting_stabilizer_x=False,
                y_basis=True,
            )
            qubit_coords_control: dict[Coord, Label] = y_geometry_control.coords

            y_geometry_target = SurfaceGeometry(
                distance=distance,
                offset=(distance * 2) + 0j,
                starting_stabilizer_x=False,
                y_basis=True,
            )
            qubit_coords_target: dict[Coord, Label] = y_geometry_target.coords

            # Define configs for control and target
            cfg_y_control = SurfaceConfig(
                distance=distance,
                state_init="+i" if key[0] == "Y+" else "-i",
                obs="Y",
                rounds=distance,
            )
            cfg_y_target = SurfaceConfig(
                distance=distance,
                state_init="+i" if key[1] == "Y+" else "-i",
                obs="Y",
                rounds=distance,
            )

        elif key[1] in {"Y+", "Y-"}:
            # Only target is Y-basis
            y_geometry_target = SurfaceGeometry(
                distance=distance,
                offset=(distance * 2) + 0j,
                starting_stabilizer_x=False,
                y_basis=True,
            )

            qubit_coords: dict[Coord, Label] = y_geometry_target.coords
            curr_offset = (distance * 2) + 0j

            cfg_y = SurfaceConfig(
                distance=distance,
                state_init="+i" if key[1] == "Y+" else "-i",
                obs="Y",
                rounds=distance,
            )

        elif key[0] in {"Y+", "Y-"}:
            # Only control is Y-basis
            y_geometry_control = SurfaceGeometry(
                distance=distance,
                offset=0 + (distance * 2) * 1j,
                starting_stabilizer_x=False,
                y_basis=True,
            )

            qubit_coords: dict[Coord, Label] = y_geometry_control.coords
            curr_offset = 0 + (distance * 2) * 1j

            cfg_y = SurfaceConfig(
                distance=distance,
                state_init="+i" if key[0] == "Y+" else "-i",
                obs="Y",
                rounds=distance,
            )

        if both_y:
            # Build circuits for both control and target
            build_y_circ = stim.Circuit()

            # Process control patch
            # Pairing Information for normal y basis rounds
            offset_control = 0 + (distance * 2) * 1j
            control_pairings = SurfacePairings(
                patch=qubit_coords_control,
                distance=distance,
                y_basis=True,
                offset=offset_control,
            )
            stab_to_data_control = control_pairings.get_schedule()

            # Pairing Information for switch rounds
            control_pairings_switch_xzy = SurfacePairings(
                patch=qubit_coords_control,
                distance=distance,
                y_basis=True,
                y_switch=True,
                offset=offset_control,
            )
            stab_to_data_switch_control, stab_to_data_xcy_control = (
                control_pairings_switch_xzy.get_schedule()
            )

            # Pairing Information for memory rounds
            control_pairings_memory = SurfacePairings(
                patch=qubit_coords_control,
                distance=distance,
                y_basis=True,
                y_memory=True,
                offset=offset_control,
            )
            stab_to_data_memory_control = control_pairings_memory.get_schedule()

            lct_y_control = Context(
                q2i=q2i,
                i2q=i2q,
                stab_to_data=stab_to_data_control,
                stab_to_data_modified=stab_to_data_switch_control,
                stab_to_data_modified2=stab_to_data_xcy_control,
                stab_to_data_modified3=stab_to_data_memory_control,
            )

            patches_control: dict[str, Patch] = {
                "patch": Patch.from_coords(qubit_coords_control, q2i),
            }

            build_y_circ += y_reset(
                lct=lct_y_control,
                patches=patches_control,
                cfg=cfg_y_control,
                skip_coords=True,
                offset=offset_control,
            ).circuit

            build_y_circ += y_initial(
                lct=lct_y_control,
                patch=patches_control["patch"],
                offset=offset_control,
            ).circuit

            rep_circ_control = y_repetition_circ(
                lct=lct_y_control,
                patch=patches_control["patch"],
                cfg=cfg_y_control,
                offset=offset_control,
                noise=noise,
            ).circuit

            switch_circ_control = y_switch_circ(
                lct=lct_y_control,
                patch=patches_control["patch"],
                offset=offset_control,
                cfg=cfg_y_control,
            ).circuit

            build_y_circ += rep_circ_control
            build_y_circ += switch_circ_control

            # Process target patch
            # Pairing Information for normal y basis rounds
            offset_target = (distance * 2) + 0j
            target_pairings = SurfacePairings(
                patch=qubit_coords_target,
                distance=distance,
                y_basis=True,
                offset=offset_target,
            )
            stab_to_data_target = target_pairings.get_schedule()

            # Pairing Information for switch rounds
            target_pairings_switch_xzy = SurfacePairings(
                patch=qubit_coords_target,
                distance=distance,
                y_basis=True,
                y_switch=True,
                offset=offset_target,
            )
            stab_to_data_switch_target, stab_to_data_xcy_target = (
                target_pairings_switch_xzy.get_schedule()
            )

            # Pairing Information for memory rounds
            target_pairings_memory = SurfacePairings(
                patch=qubit_coords_target,
                distance=distance,
                y_basis=True,
                y_memory=True,
                offset=offset_target,
            )
            stab_to_data_memory_target = target_pairings_memory.get_schedule()

            lct_y_target = Context(
                q2i=q2i,
                i2q=i2q,
                stab_to_data=stab_to_data_target,
                stab_to_data_modified=stab_to_data_switch_target,
                stab_to_data_modified2=stab_to_data_xcy_target,
                stab_to_data_modified3=stab_to_data_memory_target,
            )

            patches_target: dict[str, Patch] = {
                "patch": Patch.from_coords(qubit_coords_target, q2i),
            }

            build_y_circ += y_reset(
                lct=lct_y_target,
                patches=patches_target,
                cfg=cfg_y_target,
                skip_coords=True,
                offset=offset_target,
            ).circuit

            build_y_circ += y_initial(
                lct=lct_y_target,
                patch=patches_target["patch"],
                offset=offset_target,
            ).circuit

            rep_circ_target = y_repetition_circ(
                lct=lct_y_target,
                patch=patches_target["patch"],
                cfg=cfg_y_target,
                offset=offset_target,
                noise=noise,
            ).circuit

            switch_circ_target = y_switch_circ(
                lct=lct_y_target,
                patch=patches_target["patch"],
                offset=offset_target,
                cfg=cfg_y_target,
            ).circuit

            build_y_circ += rep_circ_target
            build_y_circ += switch_circ_target

            flow_circ = stim.Circuit()
            flow_circ += rep_circ_control
            flow_circ += switch_circ_control
            flow_circ += rep_circ_target
            flow_circ += switch_circ_target

            return reset_circuit, build_y_circ

        else:
            # Single Y-basis patch (either control or target)
            pairings = SurfacePairings(
                patch=qubit_coords,
                y_basis=True,
                distance=distance,
                offset=curr_offset,
            )
            stab_to_data: dict[tuple[Coord, Coord], str] = pairings.get_schedule()

            # Get Scheudle for switch rounds
            pairings_switch_xzy = SurfacePairings(
                patch=qubit_coords,
                y_basis=True,
                y_switch=True,
                distance=distance,
                offset=curr_offset,
            )
            stab_to_data_switch, stab_to_data_xcy = pairings_switch_xzy.get_schedule()

            # Get Schedule for memory rounds
            pairings_memory = SurfacePairings(
                patch=qubit_coords,
                y_basis=True,
                y_memory=True,
                distance=distance,
                offset=curr_offset,
            )
            stab_to_data_memory: dict[tuple[Coord, Coord], str] = pairings_memory.get_schedule()

            lct_y = Context(
                q2i=q2i,
                i2q=i2q,
                stab_to_data=stab_to_data,
                stab_to_data_modified=stab_to_data_switch,
                stab_to_data_modified2=stab_to_data_xcy,
                stab_to_data_modified3=stab_to_data_memory,
            )

            # Building Patch
            patches: dict[str, Patch] = {
                "patch": Patch.from_coords(qubit_coords, q2i),
            }

            build_y_circ = stim.Circuit()

            # Building needed circuits
            build_y_circ += y_reset(
                lct=lct_y,
                patches=patches,
                cfg=cfg_y,
                skip_coords=True,
                offset=curr_offset,
            ).circuit

            build_y_circ += y_initial(
                lct=lct_y,
                patch=patches["patch"],
                offset=curr_offset,
            ).circuit

            rep_circ = y_repetition_circ(
                lct=lct_y,
                patch=patches["patch"],
                cfg=cfg_y,
                offset=curr_offset,
                noise=noise,
            ).circuit

            switch_circ = y_switch_circ(
                lct=lct_y,
                patch=patches["patch"],
                offset=curr_offset,
                cfg=cfg_y,
            ).circuit

            build_y_circ += rep_circ
            build_y_circ += switch_circ

            flow_circ = stim.Circuit()
            flow_circ += rep_circ
            flow_circ += switch_circ

        #############################################
        # Adding Logical Y-Observables-Creation-Flow:
        #############################################

        """
        In the following the Y Observable is created depending on which qubit is in Y basis
        """

        # 1) Control Flow:
        if key[0] in {"Y+", "Y-"}:
            logical_x_string = []
            logical_z_string = []
            logical_y_string = distance * 2 - 1 + (distance * 4 - 1) * 1j

            # Finding logical Strings for x and z
            for imag in range(1, (distance * 2) - 1, 2):
                logical_x_string.append(q2i[distance * 2 - 1 + (imag + distance * 2) * 1j])

            for real in range(1, (distance * 2) - 1, 2):
                logical_z_string.append(q2i[real + (distance * 4 - 1) * 1j])

            # Adding logical z string
            logical_xyz_string = "*".join(
                [f"Z{idz}" for j, idz in enumerate(logical_z_string)]
                + [f"Y{q2i[logical_y_string]}"]
                + [f"X{idx}" for j, idx in enumerate(logical_x_string)],
            )

            logical_creation = f"{1} -> {logical_xyz_string}"

            (logical_creation_rec,) = flow_circ.solve_flow_measurements(
                [stim.Flow(logical_creation)],
            )

            # Adding the Observable
            rec_pos = []

            for index_creation in logical_creation_rec:
                current_rec_crea = flow_circ.num_measurements - index_creation
                rec_pos.append(-current_rec_crea)

            build_y_circ.append("OBSERVABLE_INCLUDE", [stim.target_rec(k) for k in rec_pos], 0)

        # 2) Target Flow:
        elif key[1] in {"Y+", "Y-"}:
            logical_x_string = []
            logical_z_string = []
            logical_y_string = distance * 4 - 1 + (distance * 2 - 1) * 1j

            # Finding logical Strings for x and z
            for imag in range(1, (distance * 2) - 1, 2):
                logical_x_string.append(q2i[distance * 4 - 1 + imag * 1j])

            for real in range(1, (distance * 2) - 1, 2):
                logical_z_string.append(q2i[real + distance * 2 + (distance * 2 - 1) * 1j])

            # Adding logical z string
            logical_xyz_string = "*".join(
                [f"Z{idz}" for j, idz in enumerate(logical_z_string)]
                + [f"Y{q2i[logical_y_string]}"]
                + [f"X{idx}" for j, idx in enumerate(logical_x_string)],
            )

            logical_creation = f"{1} -> {logical_xyz_string}"

            (logical_creation_rec,) = flow_circ.solve_flow_measurements(
                [stim.Flow(logical_creation)],
            )

            # Adding the Observable
            rec_pos = []

            for index_creation in logical_creation_rec:
                current_rec_crea = flow_circ.num_measurements - index_creation
                rec_pos.append(-current_rec_crea)

            build_y_circ.append("OBSERVABLE_INCLUDE", [stim.target_rec(k) for k in rec_pos], 0)

        return reset_circuit, build_y_circ
