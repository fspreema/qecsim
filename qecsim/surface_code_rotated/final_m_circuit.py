from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch, Context

Coord = complex

__all__ = ["final_m"]

def final_m(*, lct : Context, patches: dict[str, Patch], flow ,
            cfg : Config, before_m_flip_prob : float, is_flipped : bool ) -> stim.Circuit:
    
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    """
    The final measurement needs to be either in the normal basis or the swapped stab basis
    -> We do an if condition corresponding to if a flip is needed due to the needed basis measurements
    """

    #-Loading in Patches
    patch = patches["patch"]

    #-Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    distance = cfg.distance
    init_state = cfg.state_init
    log_obs = cfg.obs

    if is_flipped == False:
        stab_to_data = lct.stab_to_data
    else:
        stab_to_data = lct.stab_to_data_flipped

    #-Retrieving Data Coords
    data = patch.data
    qubit_coords = patch.coords

    #-Retrieving Index from Stabilizers of the Lattices
    if is_flipped == False:
        x_stab_index = patch.x_stab
        z_stab_index = patch.z_stab
    else:
        x_stab_index = patch.z_stab
        z_stab_index = patch.x_stab

    ######################
    # Define Final Circuit
    ######################

    final_circuit = stim.Circuit()

    #-------Adding-Before-Measurement-Flip-Prob.----------

    if before_m_flip_prob > 0:
        final_circuit.append("X_ERROR", data, before_m_flip_prob)

    #-------Continue-Circuit----------

    if log_obs in {"X"}:

        #6) Implementing Final Measurement-Round -> Detectors compromised by last round stab. measurements & Data parity checks from final measurement
        final_circuit.append("M", data)

        # -> Defining Data to measurement indexing
        Index_to_rec_data : dict[int, int] = {q: i for i, q in enumerate(reversed(data))}
        Index_to_rec_ancilla: dict[int, int] = {q: i for i, q in enumerate(reversed(x_stab_index + z_stab_index))}

        """
        Why do we only check the Z stabilizers in the last measurement round and not also the x stabilizers as we did before
        -> We need to meassure the X stabilizers in the X basis but we already meassured in the z basis (XZ do not commute)
        """
        
        for q, qtype in qubit_coords.items():

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
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "Z-STAB-BOUND-L":
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
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "Z-STAB-BOUND-R":
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
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))
    
    elif log_obs in {"Z"}:

        #6) Implementing Final Measurement-Round -> Detectors compromised by last round stab. measurements & Data parity checks from final measurement
        final_circuit.append("MX", data)

        # -> Defining Data to measurement indexing
        Index_to_rec_data : dict[int, int] = {q: i for i, q in enumerate(reversed(data))}
        Index_to_rec_ancilla: dict[int, int] = {q: i for i, q in enumerate(reversed(x_stab_index + z_stab_index))}

        for q, qtype in qubit_coords.items():

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
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "X-STAB-BOUND-U":
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
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

            elif qtype == "X-STAB-BOUND-B":
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
                last_record = [- Index_to_rec_ancilla[ancilla_index] - 1 - len(data)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 1))

    ##############################
    # Defining Logical Observables
    ##############################

    if flow == "X -> Z":
            
        if init_state in {"+", "-"}:
            if log_obs in {"Z"}:

                # Getting corresponding logical string and rec
                log_z = []

                for real in range(1, (distance * 2), 2):
                    log_z.append(q2i[real + 1j])

                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_z:
                        tar_rec.append(rec_pos)

                #final_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data) + k) for k in tar_rec], 0)

    elif flow == "X -> X":
    
        if init_state in {"+", "-"}:
            if log_obs in {"X"}:

                # Getting corresponding logical string and rec
                log_x = []

                for imag in range(1, (distance * 2), 2):
                    log_x.append(q2i[1 + imag*1j])

                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_x:
                        tar_rec.append(rec_pos)

                #final_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data) + k) for k in tar_rec], 0)

    elif flow == "Z -> X":

        if init_state in {"0", "1"}:
            if log_obs in {"X"}:

                # Getting corresponding logical string and rec
                log_x = []

                for imag in range(1, (distance * 2), 2):
                    log_x.append(q2i[1 + imag*1j])

                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_x:
                        tar_rec.append(rec_pos)

                #final_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data) + k) for k in tar_rec], 0)

    elif flow == "Z -> Z":

        if init_state in {"0", "1"}:
            if log_obs in {"Z"}:

                # Getting corresponding logical string and rec
                log_z = []

                for real in range(1, (distance * 2), 2):
                    log_z.append(q2i[real + 1j])

                tar_rec = []

                for rec_pos, index in enumerate(data):
                    if index in log_z:
                        tar_rec.append(rec_pos)

                #final_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data) + k) for k in tar_rec], 0)

    return(final_circuit)