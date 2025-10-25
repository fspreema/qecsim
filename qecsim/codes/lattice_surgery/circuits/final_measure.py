import stim

from qecsim.core.data_models import (
    ConfigLatticeSurgery as Config,
    LatticeContext,
    Patch_Ancilla,
    Patch_Control,
    Patch_Surgery,
    Patch_Target,
)


def final_m(*, lct : LatticeContext, patches: dict[str, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery], cfg : Config, flow : str,
            before_m_flip_prob : float) -> stim.Circuit:

    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    ancilla_patch = patches["ancilla"]
    target_patch = patches["target"]
    control_patch = patches["control"]

    #-Retrieving Global Infomration
    distance = cfg.distance
    q2i = lct.q2i
    control_state_init = cfg.control_state_init
    target_state_init = cfg.target_state_init

    #-Retrieving Lattice Coords
    qubit_coords_control = control_patch.coords
    qubit_coords_target = target_patch.coords

    #-Retrieving Data Coords
    data_ancilla = ancilla_patch.data
    data_control = control_patch.data
    data_target = target_patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index_control = control_patch.x_stab
    z_stab_index_control = control_patch.z_stab
    x_stab_index_target = target_patch.x_stab
    z_stab_index_target = target_patch.z_stab

    ##############################
    # Initlize Measurement Circuit
    ##############################

    measure_circuit = stim.Circuit()

    #############################
    # Meassuring all Data Qubits:
    #############################

    measure_circuit.append("TICK")

    if control_state_init in {"X+", "X-"}:
        #-------Adding measurement Flip Prob.--------------
        if before_m_flip_prob > 0:
            measure_circuit.append("X_ERROR", data_control, before_m_flip_prob)
        #--------------------------------------------------
        measure_circuit.append("MX", data_control)

    elif control_state_init in {"Z0", "Z1"}:
        #-------Adding measurement Flip Prob.--------------
        if before_m_flip_prob > 0:
            measure_circuit.append("X_ERROR", data_control, before_m_flip_prob)
        #--------------------------------------------------
        measure_circuit.append("MZ", data_control)

    if target_state_init in {"X+", "X-"}:
        #-------Adding measurement Flip Prob.--------------
        if before_m_flip_prob > 0:
            measure_circuit.append("X_ERROR", data_target, before_m_flip_prob)
        #--------------------------------------------------
        measure_circuit.append("MX", data_target)

    elif target_state_init in {"Z0", "Z1"}:
        #-------Adding measurement Flip Prob.--------------
        if before_m_flip_prob > 0:
            measure_circuit.append("X_ERROR", data_target, before_m_flip_prob)
        #--------------------------------------------------
        measure_circuit.append("MZ", data_target)

    ########################
    # Adding final Detectors
    ########################

    # -> Defining Data to measurement indexing
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    Index_to_rec_data : dict[int, int] = {q: i for i, q in enumerate(reversed(data_control + data_target))}
    Index_to_rec_ancilla: dict[int, int] = {q: i for i, q in enumerate(reversed(control_target_stabs))}

    # Target
    for q, qtype in qubit_coords_target.items():

        if target_state_init in {"Z0", "Z1"}:

            if qtype == "Z-STAB":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_upper_right] - 1,
                                -Index_to_rec_data[index_lower_left] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "Z-STAB-BOUND-L-T":
                #Needed Data Qubits
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_right = q2i[upper_right]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_right] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "Z-STAB-BOUND-R-T":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_lower_left = q2i[lower_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_lower_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

        if target_state_init in {"X+", "X-"}:

            if qtype == "X-STAB":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_upper_right] - 1,
                                -Index_to_rec_data[index_lower_left] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "X-STAB-BOUND-A-T":
                #Needed Data Qubits
                lower_right = q.real + 1 + (q.imag + 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_lower_right = q2i[lower_right]
                index_lower_left = q2i[lower_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_lower_right] - 1, -Index_to_rec_data[index_lower_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "X-STAB-BOUND-B-T":
                #Needed Data Qubits
                upper_right = (q.real + 1) + (q.imag - 1) * 1j
                upper_left = q.real - 1 + (q.imag - 1) * 1j

                #Finding the Correct Data Index
                index_upper_right = q2i[upper_right]
                index_upper_left = q2i[upper_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_right] - 1, -Index_to_rec_data[index_upper_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

    # Control
    for q, qtype in qubit_coords_control.items():

        if control_state_init in {"Z0", "Z1"}:

            if qtype == "Z-STAB":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_upper_right] - 1,
                                -Index_to_rec_data[index_lower_left] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "Z-STAB-BOUND-L-C":
                #Needed Data Qubits
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_right = q2i[upper_right]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_right] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "Z-STAB-BOUND-R-C":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_lower_left = q2i[lower_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_lower_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

        if control_state_init in {"X+", "X-"}:

            if qtype == "X-STAB":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = q.real + 1 + (q.imag - 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j
                lower_right = q.real + 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_upper_right] - 1,
                                -Index_to_rec_data[index_lower_left] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "X-STAB-BOUND-A-C":
                #Needed Data Qubits
                lower_right = q.real + 1 + (q.imag + 1) * 1j
                lower_left = q.real - 1 + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_lower_right = q2i[lower_right]
                index_lower_left = q2i[lower_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_lower_right] - 1, -Index_to_rec_data[index_lower_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "X-STAB-BOUND-B-C":
                #Needed Data Qubits
                upper_right = (q.real + 1) + (q.imag - 1) * 1j
                upper_left = q.real - 1 + (q.imag - 1) * 1j

                #Finding the Correct Data Index
                index_upper_right = q2i[upper_right]
                index_upper_left = q2i[upper_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_right] - 1, -Index_to_rec_data[index_upper_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                ancilla_index = q2i[q]
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data_control) - len(data_target) - len(data_ancilla)]

                #Combining the record targets
                final_record = current_record + last_record

                #Appending Detector
                measure_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

    ##############################
    # Defining Logical Observables
    ##############################

    #-----------------------ONLY-X-----------------------------

    if flow == "X -> XX":

        if control_state_init in {"X+", "X-"}:
            if target_state_init in {"X+", "X-"}:

                # Control stabilized by x logical
                log_x_ct = []

                for imag in range(((distance * 2) + 1), (distance * 4), 2):
                    log_x_ct.append(q2i[1 + imag*1j])

                for imag in range(1, (distance * 2), 2):
                    log_x_ct.append(q2i[((distance * 2) + 1) + imag*1j])

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_x_ct:
                        tar_rec.append(rec_pos)

                measure_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec], 0)

    if flow == "XX -> X":

        if control_state_init in {"X+", "X-"}:
            if target_state_init in {"X+", "X-"}:

                # Control stabilized by x logical
                log_x_c = []

                for imag in range(((distance * 2) + 1), (distance * 4), 2):
                    log_x_c.append(q2i[1 + imag*1j])

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_x_c:
                        tar_rec.append(rec_pos)

                measure_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec], 0)

    elif flow == "X -> X":

        if control_state_init in {"X+", "X-", "Z0", "Z1"}:
            if target_state_init in {"X+", "X-"}:

                # Control stabilized by x logical
                log_x_t = []

                for imag in range(1, (distance * 2), 2):
                    log_x_t.append(q2i[((distance * 2) + 1) + imag*1j])

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_x_t:
                        tar_rec.append(rec_pos)

                measure_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec], 0)

    #--------------------------ONLY-Z----------------------------

    elif flow == "Z -> ZZ":

        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1"}:

                # Control stabilized by x logical
                log_z_ct = []

                for real in range(1, (distance * 2), 2):
                    log_z_ct.append(q2i[real + ((distance * 2) + 1) * 1j])

                for real in range(((distance * 2) + 1), (distance * 4), 2):
                    log_z_ct.append(q2i[real + 1j])

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_z_ct:
                        tar_rec.append(rec_pos)

                measure_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(- len(data_control + data_target) + k) for k in tar_rec], 0)

    elif flow == "ZZ -> Z":

        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1"}:

                # Control stabilized by x logical
                log_z_t = []

                for real in range(((distance * 2) + 1), (distance * 4), 2):
                    log_z_t.append(q2i[real + 1j])

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_z_t:
                        tar_rec.append(rec_pos)

                measure_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(- len(data_control + data_target) + k) for k in tar_rec], 0)

    elif flow == "Z -> Z":

        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"Z0", "Z1", "X+", "X-"}:

                # Control stabilized by x logical
                log_z_c = []

                for real in range(1, (distance * 2), 2):
                    log_z_c.append(q2i[real + ((distance * 2) + 1) * 1j])

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in log_z_c:
                        tar_rec.append(rec_pos)

                measure_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(- len(data_control + data_target) + k) for k in tar_rec], 0)

    #-----------------------ONLY-ZX-MIX-------------------------

    elif flow == "ZX -> ZX":

        if control_state_init in {"Z0", "Z1"}:
            if target_state_init in {"X+", "X-"}:

                # Control stabilized by x logical
                combined_log_xz = []

                for real in range(1, (distance * 2), 2):
                    combined_log_xz.append(q2i[real + ((distance * 2) + 1) * 1j])

                # Target stabilized by z logical
                for imag in range(1, (distance * 2), 2):
                    combined_log_xz.append(q2i[((distance * 2) + 1) + imag*1j])

                tar_rec = []

                for rec_pos, index in enumerate(data_control + data_target):
                    if index in combined_log_xz:
                        tar_rec.append(rec_pos)

                measure_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(- len(data_control + data_target) + k) for k in tar_rec], 0)

    return measure_circuit
