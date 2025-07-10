import stim
import sinter
import os
import matplotlib.pyplot as plt
import pickle

__all__ = ["Rotated_Surface_Code"]

def Rotated_Surface_Code(distance: int, rounds : int, logical_z : int, before_measurement_flip_prob:float = 0.0, 
                         before_round_depol_data:float = 0.0, after_reset_flip_prob:float = 0.0,
                         after_clifford_depol:float = 0.0) -> stim.Circuit:
    '''
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
    '''

    #------------Preliminary-Setup------------

    #Checking Distance
    if distance <= 2:
        raise ValueError("Length must be a minimum of 3")
    
    #Checking Odd Distance
    if distance % 2 == 0:
        raise ValueError("Distance must be odd")

    #Implementing Coordinate for Qubits
    qubit_coords: dict[complex, str] = {}
    stab_to_data: dict[list[complex], str] = {}

    #----------Coordinate-Setup----------------

    #Data/Stab Coords -> WRONG!!!!!
    start_with_x = True

    for real in range(distance * 2):
        stab_counter = 0  # Reset for each row

        for imag in range(distance * 2):

        # Create Data Coords
            if real % 2 != 0 and imag % 2 != 0:
                coord = complex(real, imag)
                qubit_coords[coord] = "DATA"

        # Check Stabilizer Positions
            elif real % 2 == 0 and imag % 2 == 0:
                # Exclude Boundary
                if real != 0 and imag != 0:
                    coord = complex(real, imag)

                    # Decide which stabilizer to use
                    if start_with_x:
                        # Even positions are X, odd are Z
                        if stab_counter % 2 == 0:
                            qubit_coords[coord] = "X-STAB"
                        else:
                            qubit_coords[coord] = "Z-STAB"
                    else:
                        # Even positions are Z, odd are X
                        if stab_counter % 2 == 0:
                            qubit_coords[coord] = "Z-STAB"
                        else:
                            qubit_coords[coord] = "X-STAB"

                    stab_counter += 1

        # After finishing a row, flip the starting stabilizer for the next row
        if real % 2 == 0 and real != 0:
            start_with_x = not start_with_x

    #Create Boundary Stab Coords
    max_coord = 2 * distance

    if distance % 2 == 0:
        # Z-boundary stabilizers
        for y in range(4, max_coord, 4):
            coord = complex(0, y)
            qubit_coords[coord] = "Z-STAB-BOUND-L"

        for y in range(2, max_coord, 4):
            coord = complex(max_coord, y)
            qubit_coords[coord] = "Z-STAB-BOUND-R"

        # X-boundary stabilizers
        for x in range(2, max_coord, 4):
            coord = complex(x, 0)
            qubit_coords[coord] = "X-STAB-BOUND-U"

        for x in range(4, max_coord, 4):
            coord = complex(x, max_coord)
            qubit_coords[coord] = "X-STAB-BOUND-B"

    else:
        # Z-boundary stabilizers
        for y in range(2, max_coord, 4):
            coord = complex(0, y)
            qubit_coords[coord] = "Z-STAB-BOUND-L"

        for y in range(4, max_coord, 4):
            coord = complex(max_coord, y)
            qubit_coords[coord] = "Z-STAB-BOUND-R"

        # X-boundary stabilizers
        for x in range(4, max_coord, 4):
            coord = complex(x, 0)
            qubit_coords[coord] = "X-STAB-BOUND-U"

        for x in range(2, max_coord, 4):
            coord = complex(x, max_coord)
            qubit_coords[coord] = "X-STAB-BOUND-B"

    """
        BOUNDARY CX STABILIZER IMPLEMENTATION BELOW:
        -> very specific implementation! 
        -> If the ordering changes Stabilizers get non deterministic result 
           (Look at Det. Timeslice in order to see why this is the correct implementation -> CX all in the correct "pattern")
    """
    
    #Coordinates of Data qubits for each Stabilizer
    for coords,string in qubit_coords.items():
        #Already in Correct Orientation for Measurement of CX
        if string == "X-STAB":
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
            new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "1-CX"
            stab_to_data[new_cord2, coords] = "2-CX"
            stab_to_data[new_cord3, coords] = "3-CX"
            stab_to_data[new_cord4, coords] = "4-CX"

        elif string == "Z-STAB":
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
            new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "1-CX"
            stab_to_data[coords, new_cord2] = "2-CX"
            stab_to_data[coords, new_cord3] = "3-CX"
            stab_to_data[coords, new_cord4] = "4-CX"

        elif string == "Z-STAB-BOUND-L":
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "1-CX"
            stab_to_data[coords, new_cord2] = "2-CX"

        elif string == "Z-STAB-BOUND-R":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
            stab_to_data[coords, new_cord1] = "3-CX"
            stab_to_data[coords, new_cord2] = "4-CX"
        
        elif string == "X-STAB-BOUND-U":
            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data[new_cord1, coords] = "4-CX"
            stab_to_data[new_cord2, coords] = "3-CX"

        elif string == "X-STAB-BOUND-B":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
            stab_to_data[new_cord1, coords] = "2-CX"
            stab_to_data[new_cord2, coords] = "1-CX"

    #Indexing Qubits
    q2i: dict[complex, int] = {q: i for i, q in enumerate(
    sorted(qubit_coords, key=lambda v: (v.real, v.imag))
    )}

    #Reverse Indexing
    i2q: dict[int, complex] = {i: q for q, i in q2i.items()}

    #Indexing Z and X Stabilizers
    x_stab_index = [q2i[q] for q, qtype in qubit_coords.items() if qtype in {"X-STAB","X-STAB-BOUND-U", "X-STAB-BOUND-B"}]
    z_stab_index = [q2i[q] for q, qtype in qubit_coords.items() if qtype in {"Z-STAB","Z-STAB-BOUND-L", "Z-STAB-BOUND-R"}]

    #Indexing Data-Qubits
    data = [q2i[q] for q, qtype in qubit_coords.items() if qtype == "DATA"]
    
    #Building Circuit
    inital_circuit = stim.Circuit()
    round_circuit = stim.Circuit()
    final_circuit = stim.Circuit()

    #-----BUILDING-INITILIZATION-CIRCUIT----

    #Appending Coords
    for q, i in q2i.items():
        inital_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

    #Appending Reset (In Theory not needed right now -> Later: Init. in X-Basis)
    inital_circuit.append("R", data + x_stab_index + z_stab_index)

    #Create logical X Data String:
    data_log = []

    for imag in range(1,(distance * 2),2):
        data_log.append(q2i[3 + imag * 1j])

    if logical_z == 1:
        inital_circuit.append("X", data_log)

    #-------Adding-After-Reset-Flip-Prob.------------

    if after_reset_flip_prob > 0:
        inital_circuit.append("X_ERROR", data + x_stab_index + z_stab_index, after_reset_flip_prob)

    #-------Continue-Circuit------------

    inital_circuit.append("TICK")

    #-------Adding-Before-Round-Depol.-Data------------

    if before_round_depol_data > 0:
        inital_circuit.append("DEPOLARIZE1", data, before_round_depol_data)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    inital_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
        inital_circuit.append("DEPOLARIZE1", x_stab_index, after_clifford_depol)
        
    #-------Continue-Circuit------------

    inital_circuit.append("TICK")

    #2) CX Operations

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "1-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                inital_circuit.append("DEPOLARIZE2", index_pairs, after_clifford_depol)

    #-------Continue-Circuit------------

    inital_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "2-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                inital_circuit.append("DEPOLARIZE2", index_pairs, after_clifford_depol)

    #-------Continue-Circuit------------

    inital_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "3-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                inital_circuit.append("DEPOLARIZE2", index_pairs, after_clifford_depol)

    #-------Continue-Circuit------------

    inital_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)

    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "4-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                inital_circuit.append("DEPOLARIZE2", index_pairs, after_clifford_depol)

    #-------Continue-Circuit------------
    
    inital_circuit.append("TICK")

    #3) Basis/ Measurement
    inital_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
        inital_circuit.append("DEPOLARIZE1", x_stab_index, after_clifford_depol)

    #-------Continue-Circuit------------

    inital_circuit.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.------------

    if before_measurement_flip_prob > 0:
        inital_circuit.append("X_ERROR", x_stab_index + z_stab_index, before_measurement_flip_prob)

    #-------Continue-Circuit------------
    
    inital_circuit.append("MR", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if after_reset_flip_prob > 0:
        inital_circuit.append("X_ERROR", x_stab_index + z_stab_index, after_reset_flip_prob)

    #-------Continue-Circuit------------

    inital_circuit.append("TICK")

    #4) DETECTORS -> Measure only Z-Stabilizers!
    num_measurements_initial = len(z_stab_index)
    
    for index, q_index in enumerate(z_stab_index):
        current_tar = -1 * num_measurements_initial + index
        inital_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    #-----BUILDING-REPETITION-CIRC------

    #-------Adding-Before-Round-Depol.-Data------------

    if before_round_depol_data > 0:
        round_circuit.append("DEPOLARIZE1", data, before_round_depol_data)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    round_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
        round_circuit.append("DEPOLARIZE1", x_stab_index, after_clifford_depol)

    #-------Continue-Circuit------------

    round_circuit.append("TICK")

    #2) CX Operations

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            round_circuit.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "1-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                round_circuit.append("DEPOLARIZE2", index_pairs, after_clifford_depol)

    #-------Continue-Circuit------------

    round_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            round_circuit.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "2-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                round_circuit.append("DEPOLARIZE2", index_pairs, after_clifford_depol)

    #-------Continue-Circuit------------

    round_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            round_circuit.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "3-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                round_circuit.append("DEPOLARIZE2", index_pairs, after_clifford_depol)

    #-------Continue-Circuit------------

    round_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            round_circuit.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "4-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                round_circuit.append("DEPOLARIZE2", index_pairs, after_clifford_depol)

    #-------Continue-Circuit------------
    
    round_circuit.append("TICK")

    #3) Basis/ Measurement
    round_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_clifford_depol > 0:
        round_circuit.append("DEPOLARIZE1", x_stab_index, after_clifford_depol)

    #-------Continue-Circuit------------

    round_circuit.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.-------

    if before_measurement_flip_prob > 0:
        round_circuit.append("X_ERROR", x_stab_index + z_stab_index, before_measurement_flip_prob)

    #-------Continue-Circuit----------

    round_circuit.append("MR", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if after_reset_flip_prob > 0:
        round_circuit.append("X_ERROR", x_stab_index + z_stab_index, after_reset_flip_prob)

    #-------Continue-Circuit------------

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    round_circuit.append("SHIFT_COORDS", arg= (0,0,1))

    #4) Detectors
    num_measurements_repeat = len(x_stab_index + z_stab_index)

    for index, q_index in enumerate(x_stab_index + z_stab_index):
        prev_tar = -2 * num_measurements_repeat + index
        current_tar = -1 * num_measurements_repeat + index
        round_circuit.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                             (i2q[q_index].real, i2q[q_index].imag, 0))
        
    round_circuit.append("TICK")

    #5) Building Repeat Block
    final_circuit += round_circuit * (rounds - 1)

    #-------Adding-Before-Measurement-Flip-Prob.----------

    if before_measurement_flip_prob > 0:
        final_circuit.append("X_ERROR", data, before_measurement_flip_prob)

    #-------Continue-Circuit----------

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

    #7) Defining logical Operators (Final Measurement-Round)
    # -> F.ex Distance: 5
    # -> Data-Qubits [2, 3, 4, 5, 6, 12, 13, 14, 15, 16, 22, 23, 24, 25, 26, 32, 33, 34, 35, 36, 42, 43, 44, 45, 46] 
    # -> Z-Logical-Qubits [2, 12, 22, 32, 42]

    observable_targets_z = [stim.target_rec(-i - 1) for i in range(len(data)-1, 0, -distance)]
    final_circuit.append("OBSERVABLE_INCLUDE", observable_targets_z, [0])
    
    #8) Outputting final Circuit
    inital_circuit += final_circuit

    return inital_circuit