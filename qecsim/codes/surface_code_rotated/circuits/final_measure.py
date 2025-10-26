import stim

from qecsim.core.data_models import (
    CircuitResult,
    ConfigSurface as Config,
    Context,
    NoiseModel,
    Patch,
)

__all__ = ["final_m"]

# Determine the neighbouring Data qubits to each Stabilizer
OFFSETS = {
    "Z-STAB": [(-1, -1), (+1, -1), (-1, +1), (+1, +1)],
    "Z-STAB-BOUND-L": [(+1, -1), (+1, +1)],
    "Z-STAB-BOUND-R": [(-1, -1), (-1, +1)],
    "X-STAB": [(-1, -1), (+1, -1), (-1, +1), (+1, +1)],
    "X-STAB-BOUND-U": [(-1, +1), (+1, +1)],
    "X-STAB-BOUND-B": [(-1, -1), (+1, -1)],
}


def final_m(
    *,
    lct: Context,
    patches: dict[str, Patch],
    cfg: Config,
    noise: NoiseModel,
    is_flipped: bool,
) -> CircuitResult:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    """
    The final measurement needs to be either in the normal basis or the swapped stab basis
    -> We do an if condition corresponding to if a flip is needed due to the needed basis
    measurements
    """

    # -Loading in Patches
    patch = patches["patch"]

    # -Retrieving Global Infomration
    q2i = lct.q2i
    distance = cfg.distance
    init_state = cfg.state_init
    log_obs = cfg.obs

    # -Retrieving Data Coords
    data = patch.data
    qubit_coords = patch.coords

    # -Retrieving Index from Stabilizers of the Lattices
    if not is_flipped:
        x_stab_index = patch.x_stab
        z_stab_index = patch.z_stab
    else:
        x_stab_index = patch.z_stab
        z_stab_index = patch.x_stab

    ######################
    # Define Final Circuit
    ######################

    final_circuit = stim.Circuit()

    # -------Adding-Before-Measurement-Flip-Prob.----------
    if noise.before_m_flip_prob > 0:
        final_circuit.append("X_ERROR", data, noise.before_m_flip_prob)

    # -------Continue-Circuit----------

    ########################################################
    # Adding mutliple measurement rounds for fault tolerance
    ########################################################

    """
    If measuremnt in same basis as initlized, only one readout is neccessary
    -> Now we implement d rounds of data measurement for non inital basis measurements
    -> Guarentees faul tolerance (hopefully)
    """

    def _logical_x_indices() -> list[int]:
        # Vertical string at x=1 (odd grid), along imag axis
        return [q2i[1 + imag * 1j] for imag in range(1, 2 * distance, 2)]

    def _logical_z_indices() -> list[int]:
        # Horizontal string at y=1, along real axis
        return [q2i[real + 1j] for real in range(1, 2 * distance, 2)]

    def _logical_y_indices() -> list[int]:
        # Both
        z_string = [q2i[real + 1j] for real in range(3, 2 * distance, 2)]
        x_string = [q2i[1 + imag * 1j] for imag in range(3, 2 * distance, 2)]
        y_string = [q2i[1 + 1j]]

        return (x_string, y_string, z_string)

    # Adding final measurement of all Data qubits
    m_circ_x = stim.Circuit()
    m_circ_x.append("MX", data)

    m_circ_z = stim.Circuit()
    m_circ_z.append("MZ", data)

    # Getting corresponding logical string and rec for Y Meas (Non Det)

    m_circ_y = stim.Circuit()
    m_circ_y.append("MY", _logical_y_indices()[1])
    m_circ_y.append("MX", _logical_y_indices()[0])
    m_circ_y.append("MZ", _logical_y_indices()[2])

    # Adding final M round for deterministic measurement basis
    if not is_flipped:
        if init_state in {"1", "0"} and log_obs == "Z":
            final_circuit += m_circ_z

        elif init_state in {"+", "-"} and log_obs == "X":
            final_circuit += m_circ_x

    # Defining Data to measurement indexing
    index_to_rec_data: dict[int, int] = {q: i for i, q in enumerate(reversed(data))}
    index_to_rec_ancilla: dict[int, int] = {
        q: i
        for i, q in enumerate(
            reversed(x_stab_index + z_stab_index),
        )
    }

    # Select which stabilizers form the final detectors
    if is_flipped:
        # flipped: X-basis => use Z stabs; Z-basis => use X stabs
        commuting_stabs = (
            {"Z-STAB", "Z-STAB-BOUND-L", "Z-STAB-BOUND-R"}
            if log_obs in {"X"}
            else {"X-STAB", "X-STAB-BOUND-U", "X-STAB-BOUND-B"}
        )

    else:
        # not flipped: same family as the measurement basis
        commuting_stabs = (
            {"Z-STAB", "Z-STAB-BOUND-L", "Z-STAB-BOUND-R"}
            if log_obs in {"Z"}
            else {"X-STAB", "X-STAB-BOUND-U", "X-STAB-BOUND-B"}
        )

    ###############################
    # Unified detector construction
    ###############################
    """
    Here we construct the stabilizers from the last mr round, 
    i.e. measurement of Data is copared with measurement of the ancilla
    -> Results in 5 components that the detector get build up from

    Non-Determinstic measurements:
        * Cannot build these detectors as they are not determinstic 
        (Due to emasurement in different basis) + They are fundamentally not fault tolerant!!
    """

    # 1) Det circ for comparing stabilizers to data redout -> determinstic measurement basis

    det_circ1 = stim.Circuit()

    for q, qtype in qubit_coords.items():
        if qtype not in commuting_stabs:
            continue

        # Needed Data Qubits
        neighbor_coords = [(q.real + dx) + (q.imag + dy) * 1j for dx, dy in OFFSETS[qtype]]

        # Finding the Correct Data Index
        data_indices = [q2i[c] for c in neighbor_coords]

        # Defining the current record targets
        current_record = [-index_to_rec_data[idx] - 1 for idx in data_indices]

        # Defining the last record targets (Normal Detectors from last round)
        ancilla_index = q2i[q]
        last_record = [-index_to_rec_ancilla[ancilla_index] - 1 - len(data)]

        # Combining the record targets
        final_record = current_record + last_record

        # Appending Detector
        det_circ1.append(
            "DETECTOR", [stim.target_rec(i) for i in final_record], arg=(q.real, q.imag, 1),
        )

    if not is_flipped:
        # Det Observable
        if (init_state in {"1", "0"} and log_obs == "Z") or (
            init_state in {"+", "-"} and log_obs == "X"
        ):
            final_circuit += det_circ1

        # Non-Det Observable
        elif init_state in {"1", "0"} and log_obs == "X":
            final_circuit += m_circ_x

        elif init_state in {"+", "-"} and log_obs == "Z":
            final_circuit += m_circ_z

        elif init_state in {"+", "-", "0", "1"} and log_obs == "Y":
            final_circuit += m_circ_y

    ##############################
    # Defining Logical Observables
    ##############################

    if is_flipped:
        if init_state in {"+", "-"}:
            if log_obs == "Z":
                # Getting corresponding logical string and rec
                log_x = _logical_x_indices()

                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_x:
                        tar_rec.append(rec_pos)

                final_circuit.append(
                    "OBSERVABLE_INCLUDE", [stim.target_rec(-len(data) + k) for k in tar_rec], 0,
                )

        elif init_state in {"0", "1"}:
            if log_obs == "X":
                # Getting corresponding logical string and rec
                log_z = _logical_z_indices()

                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_z:
                        tar_rec.append(rec_pos)

                final_circuit.append(
                    "OBSERVABLE_INCLUDE", [stim.target_rec(-len(data) + k) for k in tar_rec], 0,
                )

    elif not is_flipped:
        if init_state in {"+", "-"}:
            if log_obs == "X":
                # Getting corresponding logical string and rec
                log_x = _logical_x_indices()

                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_x:
                        tar_rec.append(rec_pos)

                final_circuit.append(
                    "OBSERVABLE_INCLUDE", [stim.target_rec(-len(data) + k) for k in tar_rec], 0,
                )

            if log_obs == "Z":
                # Getting corresponding logical string and rec
                log_z = _logical_z_indices()

                final_circuit.append("OBSERVABLE_INCLUDE", [f"Z{index}" for index in log_z], 0)

                # For later decoding we need the measurement record postiions of the logical
                # operator
                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_z:
                        tar_rec.append(rec_pos)

                rec_list = [-len(data) + k for k in tar_rec]

                return CircuitResult(circuit=final_circuit, obs_indices=rec_list)

            if log_obs in {"Y"}:
                final_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [f"X{index}" for index in _logical_y_indices()[0]]
                    + [f"Y{index}" for index in _logical_y_indices()[1]]
                    + [f"Z{index}" for index in _logical_y_indices()[2]],
                    0,
                )

                # For later decoding we need the measurement record postiions of the logical
                # operator
                tar_rec = []

                for rec_pos, _index in enumerate(
                    _logical_y_indices()[0] + _logical_y_indices()[1] + _logical_y_indices()[2],
                ):
                    tar_rec.append(rec_pos)

                rec_list = [-k - 1 for k in tar_rec]

                return CircuitResult(circuit=final_circuit, obs_indices=rec_list)

        elif init_state in {"0", "1"}:
            if log_obs == "Z":
                # Getting corresponding logical string and rec
                log_z = _logical_z_indices()

                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_z:
                        tar_rec.append(rec_pos)

                final_circuit.append(
                    "OBSERVABLE_INCLUDE", [stim.target_rec(-len(data) + k) for k in tar_rec], 0,
                )

            if log_obs in {"X"}:
                # Getting corresponding logical string and rec
                log_x = _logical_x_indices()

                final_circuit.append("OBSERVABLE_INCLUDE", [f"X{index}" for index in log_x], 0)

                # For later decoding we need the measurement record postiions of the logical
                # operator
                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_x:
                        tar_rec.append(rec_pos)

                rec_list = [-len(data) + k for k in tar_rec]

                return CircuitResult(circuit=final_circuit, obs_indices=rec_list)

            if log_obs in {"Y"}:
                final_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [f"X{index}" for index in _logical_y_indices()[0]]
                    + [f"Y{index}" for index in _logical_y_indices()[1]]
                    + [f"Z{index}" for index in _logical_y_indices()[2]],
                    0,
                )

                # For later decoding we need the measurement record postiions of the logical
                # operator
                tar_rec = []

                for rec_pos, _index in enumerate(
                    _logical_y_indices()[0] + _logical_y_indices()[1] + _logical_y_indices()[2],
                ):
                    tar_rec.append(rec_pos)

                rec_list = [-k - 1 for k in tar_rec]

                return CircuitResult(circuit=final_circuit, obs_indices=rec_list)

    return CircuitResult(circuit=final_circuit)
