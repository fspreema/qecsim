import stim

__all__ = ["repetition_code"]


def repetition_code(
    distance: int,
    rounds: int,
    before_measurement_flip_prob: float = 0.0,
    before_round_data_depol: float = 0.0,
    before_round_data_flip: float = 0.0,
    before_round_data_y: float = 0.0,
    before_round_data_z: float = 0.0,
) -> stim.Circuit:
    """Created manual Repetition Code implementation

    Args:
    -> Distance (Int): Gives Distance of the Code i.e. the number of data qubits
    -> rounds (Int): Gives the number of Repeating rounds
    -> Noise (Float): Gives the probability of an initilization-Error at the beginning for the data qubits

    Returns:
    -> Stim.Circuit() with logical Z Operator
    """
    # -------Priliminary-Stup--------

    # Check Distance
    if distance % 2 == 0:
        return IndexError

    # Defining ancilla positions:
    index_ancilla = []
    for i in range(1, distance * 2, 2):
        if len(index_ancilla) != distance - 1:
            index_ancilla.append(i)

    # Defining Data positions:
    index_data = []
    for i in range(0, distance * 2, 2):
        if len(index_data) != distance:
            index_data.append(i)

    # Defining CX-Pairs
    cx_pairs: dict[list, str] = {}
    for i in index_ancilla:
        cx_pairs[i - 1, i] = "UPPER"
        cx_pairs[i + 1, i] = "LOWER"

    # Initlize Circuit
    round_circuit = stim.Circuit()
    final_circuit = stim.Circuit()

    # -------Initilization-Block---------

    # Qubit Coords
    for i, _q in enumerate(index_ancilla + index_data):
        final_circuit.append("QUBIT_COORDS", [i], [i, 0])

    final_circuit.append("R", index_data)
    final_circuit.append("R", index_ancilla)

    final_circuit.append("TICK")

    # Adding Depolarize Noise after init
    if before_round_data_depol > 0.0:
        final_circuit.append("DEPOLARIZE1", index_data, before_round_data_depol)
        final_circuit.append("TICK")

    # Adding Flip Noise after init
    if before_round_data_flip > 0.0:
        final_circuit.append("X_ERROR", index_data, before_round_data_flip)
        final_circuit.append("TICK")

    # Adding Y Noise after init
    if before_round_data_y > 0.0:
        final_circuit.append("Y_ERROR", index_data, before_round_data_y)
        final_circuit.append("TICK")

    # Adding Z Noise after init
    if before_round_data_z > 0.0:
        final_circuit.append("Z_ERROR", index_data, before_round_data_z)
        final_circuit.append("TICK")

    # Adding Upper CX
    for index, qtype in cx_pairs.items():
        if qtype == "UPPER":
            final_circuit.append("CX", index)

    final_circuit.append("TICK")

    # Adding Lower CX
    for index, qtype in cx_pairs.items():
        if qtype == "LOWER":
            final_circuit.append("CX", index)

    final_circuit.append("TICK")

    # Adding before Measurement flip Noise
    if before_measurement_flip_prob > 0.0:
        final_circuit.append("X_ERROR", index_ancilla, before_measurement_flip_prob)

    # Adding Measurements/Detectors
    final_circuit.append("MR", index_ancilla)

    for i, q in enumerate(index_ancilla):
        current_tar = -len(index_ancilla) + i
        final_circuit.append("DETECTOR", [stim.target_rec(current_tar)], arg=(q, 0, 0))

    final_circuit.append("TICK")

    # -----------Repeat-Block-----------

    # Adding Depolarize Noise after init
    if before_round_data_depol > 0.0:
        round_circuit.append("DEPOLARIZE1", index_data, before_round_data_depol)
        round_circuit.append("TICK")

    # Adding Flip Noise after init
    if before_round_data_flip > 0.0:
        round_circuit.append("X_ERROR", index_data, before_round_data_flip)
        round_circuit.append("TICK")

    # Adding Y Noise after init
    if before_round_data_y > 0.0:
        round_circuit.append("Y_ERROR", index_data, before_round_data_y)
        round_circuit.append("TICK")

    # Adding Z Noise after init
    if before_round_data_z > 0.0:
        round_circuit.append("Z_ERROR", index_data, before_round_data_z)
        round_circuit.append("TICK")

    # Adding Upper CX
    for index, qtype in cx_pairs.items():
        if qtype == "UPPER":
            round_circuit.append("CX", index)

    round_circuit.append("TICK")

    # Adding Lower CX
    for index, qtype in cx_pairs.items():
        if qtype == "LOWER":
            round_circuit.append("CX", index)

    round_circuit.append("TICK")

    # Adding before Measurement flip Noise
    if before_measurement_flip_prob > 0.0:
        round_circuit.append("X_ERROR", index_ancilla, before_measurement_flip_prob)

    # Adding Measurements
    round_circuit.append("MR", index_ancilla)

    # Adding Detectors
    # -> Shifting Detectors to create Detector Graph Volume for temporal detection
    round_circuit.append("SHIFT_COORDS", arg=(0, 0, 1))

    for i, q in enumerate(index_ancilla):
        previous_tar = -2 * len(index_ancilla) + i
        current_tar = -len(index_ancilla) + i
        round_circuit.append(
            "DETECTOR",
            [stim.target_rec(current_tar), stim.target_rec(previous_tar)],
            arg=(q, 0, 0),
        )

    round_circuit.append("TICK")

    # Construct final repetition Block
    final_circuit += round_circuit * (rounds - 1)

    # -------Logical-Observables--------

    # Adding before Measurement flip Noise
    if before_measurement_flip_prob > 0.0:
        final_circuit.append("X_ERROR", index_data, before_measurement_flip_prob)

    final_circuit.append("M", index_data)

    # Implementing last Detector round
    for i, q in enumerate(index_data):
        if i != 0:
            prev_dec = -len(index_data) + i - len(index_ancilla) - 1
            upper_meas = -len(index_data) + i - 1
            lower_meas = -len(index_data) + i
            final_circuit.append(
                "DETECTOR",
                [stim.target_rec(prev_dec), stim.target_rec(upper_meas), stim.target_rec(lower_meas)],
                arg=(q, 0, 1),
            )

    final_circuit.append("TICK")

    # Adding logical Observables
    final_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-i - 1) for i, q in enumerate(index_data)], 0)

    return final_circuit
