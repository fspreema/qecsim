from itertools import chain
import stim
from .dataclasses import Config, Patch, Context, NoiseModel, CircuitResult

Coord = complex

__all__ = ["y_final_circ"]

def y_final_circ(*, 
                 lct: Context, 
                 patches: dict[str, Patch], 
                 cfg: Config, 
                 noise: NoiseModel) -> CircuitResult:
    
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    patch = patches["patch"]

    #-Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    rounds = cfg.rounds
    distance = cfg.distance
    stab_to_data = lct.stab_to_data
    stab_to_data_switch = lct.stab_to_data_modified
    stab_to_data_switch_xcy = lct.stab_to_data_modified2

    #-Retrieving Data Coords
    data = patch.data

    # Getting infromation about the additional newly formed boundary operators:
    r_h_stabs = patch.right_h
    u_h_stabs = patch.upper_h

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab
    switch_stab_apply_h = patch.stab_switch_apply_h

    # Finding Upper right qubit index -> need to look in 2-CX
    y_coords = 1 + 1j
    y_index = q2i[y_coords]

    #################################
    # Creating final MPP measurements
    #################################

    final_measurement = stim.Circuit()
    final_measurement.append("TICK")

    #####################################
    # Define First Stab X MPP measurement
    #####################################

    reversed_switch_circ = stim.Circuit()

    """
    The following procedure will be done
    -> Along the mirrored diagonal one half h gate one half not
    -> Boundary on the site of the h gate gets expanded
        -> Every postion is now a boundary! (No 2 coords distance between them)
    -> Along the Y digaonal we do SQRT_X_DAG on all ancilla qubits
    """

    # Getting infromation about the additional newly formed boundary operators:
    r_h_stabs = patch.right_h
    u_h_stabs = patch.upper_h

    # Defining new diagonal and the corresponding stabilizers after the switch:
    
    index_h = []
    index_x_deg = []
    index_nh = []

    for cords, qtype in patch.coords.items():

        if cords != y_coords:

            # Diagonal Cut
            if cords.real > cords.imag:
                index_h.append(q2i[cords])

            # Filtering out the X_DAG -> Not on Data
            elif cords.real == cords.imag:
                if qtype != "DATA":
                    index_x_deg.append(q2i[cords])

            else:
                index_nh.append(q2i[cords])

    h_gates_rep : list = []

    for qubit, q_type in patch.coords.items():
        if q_type in {"X-STAB", "X-STAB-BOUND-B"}:
            if q2i[qubit] in index_h:
                continue
            else:
                h_gates_rep.append(q2i[qubit])

        elif q_type == "Z-STAB":
            if q2i[qubit] in index_h:
                h_gates_rep.append(q2i[qubit])
            else:
                continue

        if q_type in {"Z-STAB-BOUND-U-H", "Z-STAB-BOUND-U"}:
            h_gates_rep.append(q2i[qubit])

    reversed_switch_circ.append("H", h_gates_rep)
    reversed_switch_circ.append("TICK")

    #########################################
    # Adding the XCY gates after the H switch
    #########################################

    reversed_switch_circ.append("SQRT_X_DAG", index_x_deg)
    reversed_switch_circ.append("TICK")

    #Adding half diagonal H
    reversed_switch_circ.append("H", index_h)
    reversed_switch_circ.append("TICK")

    #2) CX Operations

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "5TICK":

            if len(coord_pairs) == 2:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                reversed_switch_circ.append("CX", index_pairs)
            else:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                reversed_switch_circ.append("CX", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "4TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("CX", index_pairs)
   
    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3.5TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("XCY", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("CX", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "2TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("CX", index_pairs)

    reversed_switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch_xcy.items():
   
        #Parallel Implementation of CX
        if order == "1TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            reversed_switch_circ.append("XCY", index_pairs)

    reversed_switch_circ.append("TICK")

    #-------Continue-Circuit------------
    


    # Adding stabs h
    reversed_switch_circ.append("H", x_stab_index + r_h_stabs)
    reversed_switch_circ.append("TICK")

    #Adding Resets
    reversed_switch_circ.append("MZ", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    reversed_switch_circ.append("SHIFT_COORDS", arg = (0,0,1))

    ########################
    # Adding final Detectors
    ########################

    ###############################
    # Define Observable measurement
    ###############################

    logical_x_string = []
    logical_z_string = []

    # Finding logical Strings for x and z
    for imag in range(3, distance * 2, 2):
        logical_x_string.append(q2i[1 + 1j * imag])

    for real in range(3, distance * 2, 2):
        logical_z_string.append(q2i[real + 1j])

    # Adding logical z string
    targets = []
    for j, idz in enumerate(logical_z_string):
        targets.append(stim.target_z(idz))
        if j < len(logical_z_string) - 1:
            targets.append(stim.target_combiner())

    # Adding single Y index
    targets.append(stim.target_combiner())
    targets.append(stim.target_y(y_index))
    targets.append(stim.target_combiner())

    # Adding logical x string
    for j, idx in enumerate(logical_x_string):
        targets.append(stim.target_x(idx))
        if j < len(logical_x_string) - 1:
            targets.append(stim.target_combiner())

    reversed_switch_circ.append("TICK")

    reversed_switch_circ.append("MPP", targets)

    ###########################
    # Adding logical Observable
    ###########################

    logical_x_string = []
    logical_z_string = []

    # Finding logical Strings for x and z
    for imag in range(3, distance * 2, 2):
        logical_x_string.append(q2i[1 + 1j * imag])

    for real in range(3, distance * 2, 2):
        logical_z_string.append(q2i[real + 1j])

    rec_meas = [- i for i, q in enumerate(logical_z_string + [y_index] + logical_x_string)]

    reversed_switch_circ.append("OBSERVABLE_INCLUDE", [stim.target_rec(i - 1) for i in rec_meas], 0)

    return CircuitResult(circuit=reversed_switch_circ)

    """for q_coord, index in patch.coords.items():
        
        # Filtering out most right x stabs
        if q_coord.real in [i for i in range(2, distance * 2, 4)]:

            #First find the paris for the x_stabs
            if index == "X-STAB":

                new_cord1 = (q_coord.real + 1) + (q_coord.imag - 1) * 1j
                new_cord2 = (q_coord.real - 1) + (q_coord.imag - 1) * 1j
                new_cord3 = (q_coord.real + 1) + (q_coord.imag + 1) * 1j
                new_cord4 = (q_coord.real - 1) + (q_coord.imag + 1) * 1j

                x_pauli_index = [q2i[coord] for coord in [new_cord1, new_cord2, new_cord3, new_cord4]]

                # Build targets list with combiners
                targets = []
                for j, idx in enumerate(x_pauli_index):
                    targets.append(stim.target_x(idx))
                    if j < len(x_pauli_index) - 1:
                        targets.append(stim.target_combiner())

                final_measurement.append("MPP", targets)

                # Adding corresponding detecotrs
                for i, q in enumerate(x_pauli_index):
                    current_tar = - len(x_pauli_index) + i

                    for i_2,q_2 in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
                        
                        if i_2 == i:
                            previous_tar = - len(x_pauli_index) - len(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs) + i_2

                            #final_measurement.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_tar)])

            elif index == "Z-STAB-BOUND-U-H":
                
                new_cord3 = (q_coord.real + 1) + (q_coord.imag + 1) * 1j
                new_cord4 = (q_coord.real - 1) + (q_coord.imag + 1) * 1j

                x_pauli_index = [q2i[coord] for coord in [new_cord3, new_cord4]]

                # Build targets list with combiners
                targets = []
                for j, idx in enumerate(x_pauli_index):
                    targets.append(stim.target_x(idx))
                    if j < len(x_pauli_index) - 1:
                        targets.append(stim.target_combiner())

                final_measurement.append("MPP", targets)

    final_measurement.append("TICK")

    ######################################
    # Define Second Stab X MPP measurement
    ######################################

    for q_coord, index in patch.coords.items():
        
        # Filtering out most right x stabs
        if q_coord.real in [i for i in range(4, distance * 2, 4)]:

            #First find the paris for the x_stabs
            if index == "X-STAB":

                new_cord1 = (q_coord.real + 1) + (q_coord.imag - 1) * 1j
                new_cord2 = (q_coord.real - 1) + (q_coord.imag - 1) * 1j
                new_cord3 = (q_coord.real + 1) + (q_coord.imag + 1) * 1j
                new_cord4 = (q_coord.real - 1) + (q_coord.imag + 1) * 1j

                x_pauli_index = [q2i[coord] for coord in [new_cord1, new_cord2, new_cord3, new_cord4]]

                # Build targets list with combiners
                targets = []
                for j, idx in enumerate(x_pauli_index):
                    targets.append(stim.target_x(idx))
                    if j < len(x_pauli_index) - 1:
                        targets.append(stim.target_combiner())

                final_measurement.append("MPP", targets)

            elif index == "X-STAB-BOUND-B":
                
                new_cord1 = (q_coord.real + 1) + (q_coord.imag - 1) * 1j
                new_cord2 = (q_coord.real - 1) + (q_coord.imag - 1) * 1j

                x_pauli_index = [q2i[coord] for coord in [new_cord1, new_cord2]]

                # Build targets list with combiners
                targets = []
                for j, idx in enumerate(x_pauli_index):
                    targets.append(stim.target_x(idx))
                    if j < len(x_pauli_index) - 1:
                        targets.append(stim.target_combiner())

                final_measurement.append("MPP", targets)

    final_measurement.append("TICK")

    #####################################
    # Define First Stab Z MPP measurement
    #####################################

    for q_coord, index in patch.coords.items():
        
        # Filtering out most right x stabs
        if q_coord.imag in [i for i in range(4, distance * 2, 4)]:

            #First find the paris for the x_stabs
            if index == "Z-STAB":

                new_cord1 = (q_coord.real + 1) + (q_coord.imag - 1) * 1j
                new_cord2 = (q_coord.real - 1) + (q_coord.imag - 1) * 1j
                new_cord3 = (q_coord.real + 1) + (q_coord.imag + 1) * 1j
                new_cord4 = (q_coord.real - 1) + (q_coord.imag + 1) * 1j

                z_pauli_index = [q2i[coord] for coord in [new_cord1, new_cord2, new_cord3, new_cord4]]

                # Build targets list with combiners
                targets = []
                for j, idz in enumerate(z_pauli_index):
                    targets.append(stim.target_z(idz))
                    if j < len(z_pauli_index) - 1:
                        targets.append(stim.target_combiner())

                final_measurement.append("MPP", targets)

            elif index == "Z-STAB-BOUND-L":
                
                new_cord1 = (q_coord.real + 1) + (q_coord.imag - 1) * 1j
                new_cord3 = (q_coord.real + 1) + (q_coord.imag + 1) * 1j

                z_pauli_index = [q2i[coord] for coord in [new_cord1, new_cord3]]

                # Build targets list with combiners
                targets = []
                for j, idz in enumerate(z_pauli_index):
                    targets.append(stim.target_z(idz))
                    if j < len(z_pauli_index) - 1:
                        targets.append(stim.target_combiner())

                final_measurement.append("MPP", targets)

    final_measurement.append("TICK")

    ######################################
    # Define Second Stab z MPP measurement
    ######################################

    for q_coord, index in patch.coords.items():
        
        # Filtering out most right x stabs
        if q_coord.imag in [i for i in range(2, distance * 2, 4)]:

            #First find the paris for the x_stabs
            if index == "Z-STAB":

                new_cord1 = (q_coord.real + 1) + (q_coord.imag - 1) * 1j
                new_cord2 = (q_coord.real - 1) + (q_coord.imag - 1) * 1j
                new_cord3 = (q_coord.real + 1) + (q_coord.imag + 1) * 1j
                new_cord4 = (q_coord.real - 1) + (q_coord.imag + 1) * 1j

                z_pauli_index = [q2i[coord] for coord in [new_cord1, new_cord2, new_cord3, new_cord4]]

                # Build targets list with combiners
                targets = []
                for j, idz in enumerate(z_pauli_index):
                    targets.append(stim.target_z(idz))
                    if j < len(z_pauli_index) - 1:
                        targets.append(stim.target_combiner())

                final_measurement.append("MPP", targets)

            elif index == "X-STAB-BOUND-R-H":
                
                new_cord2 = (q_coord.real - 1) + (q_coord.imag - 1) * 1j
                new_cord4 = (q_coord.real - 1) + (q_coord.imag + 1) * 1j

                z_pauli_index = [q2i[coord] for coord in [new_cord2, new_cord4]]

                # Build targets list with combiners
                targets = []
                for j, idz in enumerate(z_pauli_index):
                    targets.append(stim.target_z(idz))
                    if j < len(z_pauli_index) - 1:
                        targets.append(stim.target_combiner())

                final_measurement.append("MPP", targets)

    final_measurement.append("TICK")


    ###############################
    # Define Observable measurement
    ###############################

    logical_x_string = []
    logical_z_string = []

    # Finding logical Strings for x and z
    for imag in range(3, distance * 2, 2):
        logical_x_string.append(q2i[1 + 1j * imag])

    for real in range(3, distance * 2, 2):
        logical_z_string.append(q2i[real + 1j])

    # Adding logical z string
    targets = []
    for j, idz in enumerate(logical_z_string):
        targets.append(stim.target_z(idz))
        if j < len(logical_z_string) - 1:
            targets.append(stim.target_combiner())

    # Adding single Y index
    targets.append(stim.target_combiner())
    targets.append(stim.target_y(y_index))
    targets.append(stim.target_combiner())

    # Adding logical x string
    for j, idx in enumerate(logical_x_string):
        targets.append(stim.target_x(idx))
        if j < len(logical_x_string) - 1:
            targets.append(stim.target_combiner())

    final_measurement.append("MPP", targets)

    ###################################
    # Getting Observable with stim.flow
    ###################################

    logical_x_string = []
    logical_z_string = []

    # Finding logical Strings for x and z
    for imag in range(3, distance * 2, 2):
        logical_x_string.append(q2i[1 + 1j * imag])

    for real in range(3, distance * 2, 2):
        logical_z_string.append(q2i[real + 1j])

    rec_meas = [- i for i, q in enumerate(logical_z_string + [y_index] + logical_x_string)]

    #final_measurement.append("OBSERVABLE_INCLUDE", [stim.target_rec(i - 1) for i in rec_meas], 0)"""

    return final_measurement