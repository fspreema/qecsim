import stim

from src.core.data_models import ConfigXZZX, PatchXZZX, XZZXContext

__all__ = ["final_measure"]


def final_measure(
    *,
    lct: XZZXContext,
    cfg: ConfigXZZX,
    patch: PatchXZZX,
) -> stim.Circuit:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # Retrieving Global Infomration
    q2i = lct.q2i
    distance = cfg.distance
    state_init = cfg.state_init

    # Retrieving Stabilizer Indices
    stab_index = patch.stab_index

    # Retireiving qubit coords
    qubit_coords = patch.coords

    # Retrieving Data Indices
    data_z = patch.data_z
    data_x = patch.data_x

    ##########################
    # Adding final measurement
    ##########################

    final_circuit = stim.Circuit()

    final_circuit.append("SHIFT_COORDS", arg=(0, 0, 1))

    ###########################################
    # 6) Implementing Final Measurement-Round
    # -> Detectors compromised by last round stab. measurements
    #    & Data parity checks from final measurement
    ###########################################

    final_circuit.append("MZ", data_z)
    final_circuit.append("MX", data_x)

    # -> Defining Data to measurement indexing
    index_to_rec_data: dict[int, int] = {q: i for i, q in enumerate(reversed(data_z + data_x))}
    index_to_rec_ancilla: dict[int, int] = {q: i for i, q in enumerate(reversed(stab_index))}

    """
    Why do we only check the Z stabilizers in the last measurement round and 
    not also the x stabilizers as we did before
    -> We need to meassure the X stabilizers in the X basis but we already 
       meassured in the z basis (XZ do not commute)
    """

    if state_init in {"XZZX-VER"}:
        for q, qtype in qubit_coords.items():
            if qtype == "STAB-Ver":
                # Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                # Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                # Defining the current record targets
                current_record = [
                    -index_to_rec_data[index_upper_left] - 1,
                    -index_to_rec_data[index_upper_right] - 1,
                    -index_to_rec_data[index_lower_left] - 1,
                    -index_to_rec_data[index_lower_right] - 1,
                ]

                # Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [-index_to_rec_ancilla[ancilla_index] - 1 - len(data_z + data_x)]

                # Combining the record targets
                final_record = current_record + last_record

                # Appending Detector
                final_circuit.append(
                    "DETECTOR",
                    [stim.target_rec(i) for i in final_record],
                    arg=(q.real, q.imag, 1),
                )

            elif qtype == "STAB-BOUND-A-Ver":
                # Needed Data Qubits
                lower_right = q.real + 1 + (q.imag + 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j

                # Finding the Correct Data Index
                index_lower_right = q2i[lower_right]
                index_lower_left = q2i[lower_left]

                # Defining the current record targets
                current_record = [
                    -index_to_rec_data[index_lower_right] - 1,
                    -index_to_rec_data[index_lower_left] - 1,
                ]

                # Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [-index_to_rec_ancilla[ancilla_index] - 1 - len(data_z + data_x)]

                # Combining the record targets
                final_record = current_record + last_record

                # Appending Detector
                final_circuit.append(
                    "DETECTOR",
                    [stim.target_rec(i) for i in final_record],
                    arg=(q.real, q.imag, 1),
                )

            elif qtype == "STAB-BOUND-B-Ver":
                # Needed Data Qubits
                upper_right = (q.real + 1) + (q.imag - 1) * 1j
                upper_left = q.real - 1 + (q.imag - 1) * 1j

                # Finding the Correct Data Index
                index_upper_right = q2i[upper_right]
                index_upper_left = q2i[upper_left]

                # Defining the current record targets
                current_record = [
                    -index_to_rec_data[index_upper_right] - 1,
                    -index_to_rec_data[index_upper_left] - 1,
                ]

                # Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [-index_to_rec_ancilla[ancilla_index] - 1 - len(data_z + data_x)]

                # Combining the record targets
                final_record = current_record + last_record

                # Appending Detector
                final_circuit.append(
                    "DETECTOR",
                    [stim.target_rec(i) for i in final_record],
                    arg=(q.real, q.imag, 1),
                )

    elif state_init in {"XZZX-HOR"}:
        for q, qtype in qubit_coords.items():
            if qtype == "STAB-Hor":
                # Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                # Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                # Defining the current record targets
                current_record = [
                    -index_to_rec_data[index_upper_left] - 1,
                    -index_to_rec_data[index_upper_right] - 1,
                    -index_to_rec_data[index_lower_left] - 1,
                    -index_to_rec_data[index_lower_right] - 1,
                ]

                # Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [-index_to_rec_ancilla[ancilla_index] - 1 - len(data_z + data_x)]

                # Combining the record targets
                final_record = current_record + last_record

                # Appending Detector
                final_circuit.append(
                    "DETECTOR",
                    [stim.target_rec(i) for i in final_record],
                    arg=(q.real, q.imag, 1),
                )

            elif qtype == "STAB-BOUND-L-Hor":
                # Needed Data Qubits
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                # Finding the Correct Data Index
                index_upper_right = q2i[upper_right]
                index_lower_right = q2i[lower_right]

                # Defining the current record targets
                current_record = [
                    -index_to_rec_data[index_upper_right] - 1,
                    -index_to_rec_data[index_lower_right] - 1,
                ]

                # Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [-index_to_rec_ancilla[ancilla_index] - 1 - len(data_z + data_x)]

                # Combining the record targets
                final_record = current_record + last_record

                # Appending Detector
                final_circuit.append(
                    "DETECTOR",
                    [stim.target_rec(i) for i in final_record],
                    arg=(q.real, q.imag, 1),
                )

            elif qtype == "STAB-BOUND-R-Hor":
                # Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j

                # Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_lower_left = q2i[lower_left]

                # Defining the current record targets
                current_record = [
                    -index_to_rec_data[index_upper_left] - 1,
                    -index_to_rec_data[index_lower_left] - 1,
                ]

                # Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [-index_to_rec_ancilla[ancilla_index] - 1 - len(data_z + data_x)]

                # Combining the record targets
                final_record = current_record + last_record

                # Appending Detector
                final_circuit.append(
                    "DETECTOR",
                    [stim.target_rec(i) for i in final_record],
                    arg=(q.real, q.imag, 1),
                )

    ########################################################
    # 7) Defining logical Operators (Final Measurement-Round)
    ########################################################

    if state_init in {"XZZX-VER"}:
        # Control stabilized by x logical
        log_ver = []

        for imag in range(1, (distance * 2), 2):
            log_ver.append(q2i[1 + imag * 1j])

        tar_rec = []

        for rec_pos, index in enumerate(data_z + data_x):
            if index in log_ver:
                tar_rec.append(rec_pos)

        final_circuit.append(
            "OBSERVABLE_INCLUDE",
            [stim.target_rec(-len(data_z + data_x) + k) for k in tar_rec],
            0,
        )

    elif state_init in {"XZZX-HOR"}:
        # Control stabilized by x logical
        log_hor = []

        for real in range(1, (distance * 2), 2):
            log_hor.append(q2i[real + 1j])

        tar_rec = []

        for rec_pos, index in enumerate(data_z + data_x):
            if index in log_hor:
                tar_rec.append(rec_pos)

        final_circuit.append(
            "OBSERVABLE_INCLUDE",
            [stim.target_rec(-len(data_z + data_x) + k) for k in tar_rec],
            0,
        )

    return final_circuit
