import stim

from src.codes.lattice_surgery.logical_strings import get_logical_strings
from src.core.data_models import (
    ConfigLatticeSurgery as Config,
    LatticeContext,
    PatchAncilla,
    PatchControl,
    PatchSurgery,
    PatchTarget,
)


def final_m(
    *,
    lct: LatticeContext,
    patches: dict[str, PatchAncilla, PatchControl, PatchTarget, PatchSurgery],
    cfg: Config,
    flow: str,
) -> stim.Circuit:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # -Loading in Patches
    target_patch = patches["target"]
    control_patch = patches["control"]

    # -Retrieving Global Infomration
    distance = cfg.distance
    q2i = lct.q2i
    control_state_init = cfg.control_state_init
    target_state_init = cfg.target_state_init

    # -Retrieving Data Coords
    data_control = control_patch.data
    data_target = target_patch.data

    # -Get all logical operator strings
    log_strings_shift = get_logical_strings(
        q2i,
        distance,
        shift_cx_for_y=True,
        shift_tz_for_y=True,
        shift_tx_for_y=True,
        shift_cz_for_y=True,
    )
    log_strings = get_logical_strings(q2i, distance)

    ##############################
    # Initlize Measurement Circuit
    ##############################

    measure_circuit = stim.Circuit()

    #############################
    # Meassuring all Data Qubits:
    #############################

    if control_state_init in {"Y+", "Y-"} or target_state_init in {"Y+", "Y-"}:
        # Y-including flows are handled per-flow below (measure only logical strings)
        pass
    else:
        if control_state_init in {"X+", "X-"}:
            measure_circuit.append("TICK")
            measure_circuit.append("MX", data_control)

        elif control_state_init in {"Z0", "Z1"}:
            measure_circuit.append("TICK")
            measure_circuit.append("MZ", data_control)

        if target_state_init in {"X+", "X-"}:
            measure_circuit.append("TICK")
            measure_circuit.append("MX", data_target)

        elif target_state_init in {"Z0", "Z1"}:
            measure_circuit.append("TICK")
            measure_circuit.append("MZ", data_target)

    ##############################
    # Defining Logical Observables
    ##############################

    # -----------------------ONLY-X-----------------------------

    if flow == "X -> XX":
        if control_state_init in {"X+", "X-"}:
            if target_state_init in {"X+", "X-"}:
                # Control stabilized by x logical
                log_x_ct = log_strings["c_x"] + log_strings["t_x"]

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_x_ct:
                        tar_rec.append(rec_pos)

                measure_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec],
                    0,
                )

    if flow == "XX -> X":
        if control_state_init in {"X+", "X-"}:
            if target_state_init in {"X+", "X-"}:
                # Control stabilized by x logical
                log_x_c = log_strings["c_x"]

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_x_c:
                        tar_rec.append(rec_pos)

                measure_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec],
                    0,
                )

    elif flow == "X -> X":
        if control_state_init in {"X+", "X-", "Z0", "Z1"}:
            if target_state_init in {"X+", "X-"}:
                # Control stabilized by x logical
                log_x_t = log_strings["t_x"]

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_x_t:
                        tar_rec.append(rec_pos)

                measure_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec],
                    0,
                )

    # --------------------------ONLY-Z----------------------------

    elif flow == "Z -> ZZ":
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1"}:
                # Control stabilized by x logical
                log_z_ct = log_strings["c_z"] + log_strings["t_z"]

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_z_ct:
                        tar_rec.append(rec_pos)

                measure_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec],
                    0,
                )

    elif flow == "ZZ -> Z":
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1"}:
                # Control stabilized by x logical
                log_z_t = log_strings["t_z"]

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_z_t:
                        tar_rec.append(rec_pos)

                measure_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec],
                    0,
                )

    elif flow == "Z -> Z":
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1", "X+", "X-"}:
                # Control stabilized by x logical
                log_z_c = log_strings["c_z"]

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_z_c:
                        tar_rec.append(rec_pos)

                measure_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec],
                    0,
                )

    # -----------------------ONLY-ZX-MIX-------------------------

    elif flow == "ZX -> ZX":
        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"X+", "X-"}:
                # Control stabilized by z logical, Target stabilized by x logical
                combined_log_xz = log_strings["c_z"] + log_strings["t_x"]

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in combined_log_xz:
                        tar_rec.append(rec_pos)

                measure_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec],
                    0,
                )

    # -----------------------Y-Included-measurements----------------

    elif flow == "YZ -> XY":
        # RHS: control X, target Y
        k_count = 0

        # Control X logical -> measure MX on c_x indices
        c_x = log_strings_shift["c_x"]
        measure_circuit.append("MX", c_x)
        k_count += len(c_x)

        # Target Y logical -> measure MZ on z_string, MY on corner, MX on x_string
        t_y = log_strings["t_y"]
        measure_circuit.append("MZ", t_y["z_string"])
        k_count += len(t_y["z_string"])

        measure_circuit.append("MY", t_y["y_corner"])
        k_count += len(t_y["y_corner"])

        measure_circuit.append("MX", t_y["x_string"])
        k_count += len(t_y["x_string"])

        # Include last k_count measurements as the observable
        measure_circuit.append(
            "OBSERVABLE_INCLUDE",
            [stim.target_rec(-i) for i in range(1, k_count + 1)],
            0,
        )

    elif flow == "YI -> YX":
        # RHS: control Y, target X
        k_count = 0

        # Control Y logical
        c_y = log_strings["c_y"]
        measure_circuit.append("MZ", c_y["z_string"])
        k_count += len(c_y["z_string"])

        measure_circuit.append("MY", c_y["y_corner"])
        k_count += len(c_y["y_corner"])

        measure_circuit.append("MX", c_y["x_string"])
        k_count += len(c_y["x_string"])

        # Target X logical
        t_x = log_strings["t_x"]
        measure_circuit.append("MX", t_x)
        k_count += len(t_x)

        measure_circuit.append(
            "OBSERVABLE_INCLUDE",
            [stim.target_rec(-i) for i in range(1, k_count + 1)],
            0,
        )

    elif flow == "YX -> YI":
        # RHS: control Y, target I (no target logical measurement)
        k_count = 0

        # Control Y logical
        c_y = log_strings["c_y"]
        measure_circuit.append("MZ", c_y["z_string"])
        k_count += len(c_y["z_string"])

        measure_circuit.append("MY", c_y["y_corner"])
        k_count += len(c_y["y_corner"])

        measure_circuit.append("MX", c_y["x_string"])
        k_count += len(c_y["x_string"])

        measure_circuit.append(
            "OBSERVABLE_INCLUDE",
            [stim.target_rec(-i) for i in range(1, k_count + 1)],
            0,
        )

    elif flow == "YY -> XZ":
        # RHS: control X, target Z
        k_count = 0

        # Control X logical
        c_x = log_strings["c_x"]
        measure_circuit.append("MX", c_x)
        k_count += len(c_x)

        # Target Z logical
        t_z = log_strings["t_z"]
        measure_circuit.append("MZ", t_z)
        k_count += len(t_z)

        measure_circuit.append(
            "OBSERVABLE_INCLUDE",
            [stim.target_rec(-i) for i in range(1, k_count + 1)],
            0,
        )

    elif flow == "IY -> ZY":
        # RHS: control Z, target Y
        k_count = 0

        # Control Z logical
        c_z = log_strings["c_z"]
        measure_circuit.append("MZ", c_z)
        k_count += len(c_z)

        # Target Y logical
        t_y = log_strings["t_y"]
        measure_circuit.append("MZ", t_y["z_string"])
        k_count += len(t_y["z_string"])

        measure_circuit.append("MY", t_y["y_corner"])
        k_count += len(t_y["y_corner"])

        measure_circuit.append("MX", t_y["x_string"])
        k_count += len(t_y["x_string"])

        measure_circuit.append(
            "OBSERVABLE_INCLUDE",
            [stim.target_rec(-i) for i in range(1, k_count + 1)],
            0,
        )

    elif flow == "XY -> YZ":
        # RHS: control Y, target Z
        k_count = 0

        # Control Y logical
        c_y = log_strings["c_y"]
        measure_circuit.append("MZ", c_y["z_string"])
        k_count += len(c_y["z_string"])

        measure_circuit.append("MY", c_y["y_corner"])
        k_count += len(c_y["y_corner"])

        measure_circuit.append("MX", c_y["x_string"])
        k_count += len(c_y["x_string"])

        # Target Z logical
        t_z = log_strings["t_z"]
        measure_circuit.append("MZ", t_z)
        k_count += len(t_z)

        measure_circuit.append(
            "OBSERVABLE_INCLUDE",
            [stim.target_rec(-i) for i in range(1, k_count + 1)],
            0,
        )

    return measure_circuit
