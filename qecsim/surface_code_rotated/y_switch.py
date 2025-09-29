import stim
from .data_models import Config, Patch, Context, NoiseModel, CircuitResult

Coord = complex

__all__ = ["y_switch_circ"]

def y_switch_circ(*, 
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

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab
    switch_stab_apply_h = patch.stab_switch_apply_h

    # Finding Upper right qubit index -> need to look in 2-CX
    y_coords = 1 + 1j
    y_index = q2i[y_coords]

    ###########################
    # Define Repetition Circuit
    ###########################

    #-----BUILDING-REPETITION-CIRC------

    pre_switch_circ = stim.Circuit()

    pre_switch_circ.append("TICK")
    pre_switch_circ.append("R", x_stab_index + z_stab_index)

    #-------Adding-Before-Round-Depol.-Data------------

    if noise.before_round_depol > 0:
        pre_switch_circ.append("DEPOLARIZE1", data, noise.before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    pre_switch_circ.append("TICK")
    pre_switch_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        pre_switch_circ.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    pre_switch_circ.append("TICK")

    #2) CX Operations

    #2) CX Operations

    def _pairs_for(order: str) -> list[list[int]]:
        # Return the Pairs needed at the current order
        return [[q2i[cp[1]], q2i[cp[0]]] for cp, o in stab_to_data.items() if o == order]

    def _append_by_order(op: str, order: str, noise: float = 0.0) -> None:
        # Getting pair info
        for pair in _pairs_for(order):
            # Checking for upper corner CX and leave it out!
            if y_index not in pair:
                #Adding Pair on Operation
                if op == "CX":
                    pre_switch_circ.append(op, pair)
                elif op == "DEPOLARIZE2":
                    pre_switch_circ.append(op, pair, noise)

    # Adding all the CX gates
    for order in ("1-CX", "2-CX", "3-CX", "4-CX"):
        _append_by_order("CX", order)
        if noise.after_c_depol_prob > 0:
            _append_by_order("DEPOLARIZE2", order, noise.after_c_depol_prob)
        pre_switch_circ.append("TICK")

    #-------Continue-Circuit------------

    #3) Basis/ Measurement
    pre_switch_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        pre_switch_circ.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    pre_switch_circ.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.-------

    if noise.before_m_flip_prob > 0:
        pre_switch_circ.append("X_ERROR", x_stab_index + z_stab_index, noise.before_m_flip_prob)

    #-------Continue-Circuit----------

    pre_switch_circ.append("MZ", x_stab_index + z_stab_index)

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    pre_switch_circ.append("SHIFT_COORDS", arg = (0,0,1))

    #4) Detectors
    num_measurements_repeat = len(x_stab_index + z_stab_index)

    for index, q_index in enumerate(x_stab_index + z_stab_index):
        prev_tar = -2 * num_measurements_repeat + index
        current_tar = -1 * num_measurements_repeat + index
        pre_switch_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                             (i2q[q_index].real, i2q[q_index].imag, 0))
        
    pre_switch_circ.append("TICK")

    ############################################
    # Adding the Switch after the init of the RY
    ############################################

    switch_circ = stim.Circuit()

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

    #---------Adding-RY-Gate----------------

    switch_circ.append("RY", y_index)

    #---------Adding-Resets-Old&New-Stabs---

    switch_circ.append("RZ", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)

    #-------Adding-Before-Round-Depol.-Data------------

    if noise.before_round_depol > 0:
        switch_circ.append("DEPOLARIZE1", data, noise.before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    switch_circ.append("TICK")
    switch_circ.append("H", x_stab_index + r_h_stabs)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        switch_circ.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    #-----------Adding-H-&-X-DAG-Gates------

    switch_circ.append("TICK")
    switch_circ.append("H", index_h)
    switch_circ.append("SQRT_X_DAG", index_x_deg)
    switch_circ.append("TICK")

    #########################################
    # Adding the XCY gates after the H switch
    #########################################

    #2) CX Operations

    for coord_pairs, order in stab_to_data_switch_xcy.items():
   
        #Parallel Implementation of CX
        if order == "1TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("XCY", index_pairs)

    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "2TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("CX", index_pairs)

    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("CX", index_pairs)

    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "3.5TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("XCY", index_pairs)

    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "4TICK":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            switch_circ.append("CX", index_pairs)
   
    switch_circ.append("TICK")

    for coord_pairs, order in stab_to_data_switch.items():
   
        #Parallel Implementation of CX
        if order == "5TICK":

            if len(coord_pairs) == 2:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                switch_circ.append("CX", index_pairs)
            else:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                switch_circ.append("CX", index_pairs)

    #-------Continue-Circuit------------
    
    switch_circ.append("TICK")

    #3) Basis/ Measurement
    switch_circ.append("H", switch_stab_apply_h)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        switch_circ.append("DEPOLARIZE1", x_stab_index + r_h_stabs, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    switch_circ.append("TICK")
    switch_circ.append("M", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    switch_circ.append("SHIFT_COORDS", arg = (0,0,1))

    ###########################################################
    # Adding Detectors (Newly Gen Boundary and Old Stabilizers)
    ###########################################################

    #1) Deterministic Detectors which can be build up by the old stabilizers

    previous_round = pre_switch_circ.num_measurements
    current_round = switch_circ.num_measurements

    for index, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):

        # Calc important info
        coord = i2q[q_index]
        role = patch.coords[coord]

        print(coord)
        
        if role == "X-STAB":

            # skip stabs on the Y-cut
            if q_index in index_nh:
                print(q_index)
                # find this stabilizers position in the previous block
                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))


        elif role == "Z-STAB":

            #skip stabs on the Y-cut
            if q_index in index_nh:

                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                

        elif role in {"X-STAB-BOUND-B", "Z-STAB-BOUND-L"}:
            prev_tar = - current_round - previous_round + index
            current_tar = - current_round + index
            switch_circ.append("DETECTOR",
                [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                (coord.real, coord.imag, 0))

        elif role in {"X-STAB-BOUND-R"}:
    
            """
            This Detector gets from weight 2 to weight 3 and includes the measurement from the X-STAB-BOUND-R-H
            """

            # Important Information of above ancilla
            coord_above = i2q[q_index] - 2j
            index_above = q2i[coord_above]

            # Finding record index
            for index2, q_index2 in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
                if q_index2 == index_above:
                    current_tar2 = - current_round + index2

            previous_tar = - previous_round - current_round + index
            current_tar1 = - current_round + index
            
            # Adding Detector out of all 3 measurements:
            #switch_circ.append(
            #        "DETECTOR",
            #        [stim.target_rec(current_tar1), stim.target_rec(previous_tar), stim.target_rec(current_tar2)],
            #        (coord.real, coord.imag, 0)
            #)

        #for flows in switch_circ.flow_generators():
        #    print(flows)


        #     if q2i[(coord.real - 1 + 1j * coord.imag + 1j)] in index_h:
        #         #switch_circ.append(
        #         #    "DETECTOR",
        #         #    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
        #         #    (coord.real, coord.imag, 0)
        #         #)

        # elif role in {"X-STAB-BOUND-B", "Z-STAB-BOUND-L"}:
        #     #prev_tar = - current_round - previous_round + prev_idx
        #     #current_tar = - current_round + index
        #     #switch_circ.append(
        #     #    "DETECTOR",
        #     #    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
        #     #    (coord.real, coord.imag, 0)
        #     #)

        # e

    full_switch = pre_switch_circ + switch_circ

    return CircuitResult(circuit=full_switch)

















































    """
    We need to look at the switch circuit and what detectors we need, we do this the following way

    ####################
    # FLOW MEASUREMENTS:
    ####################

    *  S -> 1 Contarction of stabilizers (These checks stop exisiting as we have a bsis switch)
    *  1 -> S Creation ofs tabilizers (These get created as we have a basis switch)
    """

    # Add Detector for the Z stabilizers on the Diagonal     
    for coords, label in patch.coords.items():
        if label in {"X-STAB"}:

            if q2i[coords] in index_h:

                new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
                new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
                new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
                new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j

                if coords in [2 + i + 1j * i for i in range(2,distance + 2, 2)]:
        
                    old_stab = [new_cord1, new_cord2, new_cord3, new_cord4]
                    new_stab_z = [new_cord2, new_cord3]
                    new_stab_y = [new_cord4]
                    
                    x_new = '*'.join([f"Y{q2i[q]}" for q in new_stab_y] + [f"X{q2i[q]}" for q in new_stab_z])
                    x_old = '*'.join(f"X{i}" for i in [q2i[current] for current in old_stab])
                    prev_round = f"{1} -> {x_old}"
                    switch_round_conc = f"{x_old} -> {1}"
                    switch_round_crea = f"{1} -> {x_new}"

                    # -> Frist Contract (Old Circuit before siwtch)
                    (included_measurements_prev,) = pre_switch_circ.solve_flow_measurements([
                        stim.Flow(prev_round),
                    ])

                    (included_measurements_conc,) = switch_circ.solve_flow_measurements([
                    stim.Flow(switch_round_conc),
                    ])

                    (included_measurements_crea,) = switch_circ.solve_flow_measurements([
                    stim.Flow(switch_round_crea),
                    ])

                    # Buidling new x stab off diagonal stabilizers
                    #print(included_measurements_prev, included_measurements_conc, included_measurements_crea)
                
                # Adding the other weight 3 Stabilizers
                else:

                    old_stab = [new_cord1, new_cord2, new_cord3, new_cord4]
                    new_stab_z = [new_cord2, new_cord3, new_cord4]
                    
                    x_new = '*'.join([f"Z{q2i[q]}" for q in new_stab_z])
                    x_old = '*'.join(f"Z{i}" for i in [q2i[current] for current in old_stab])
                    prev_round = f"{1} -> {x_old}"
                    switch_round_conc = f"{x_old} -> {1}"
                    switch_round_crea = f"{x_new} -> {1}"

                    # -> Frist Contract (Old Circuit before siwtch)
                    (included_measurements_prev,) = pre_switch_circ.solve_flow_measurements([
                    stim.Flow(prev_round),
                    ])

                    (included_measurements_conc,) = switch_circ.solve_flow_measurements([
                    stim.Flow(switch_round_conc),
                    ])

                    (included_measurements_crea,) = switch_circ.solve_flow_measurements([
                    stim.Flow(switch_round_crea),
                    ])

                    #print(included_measurements_prev, included_measurements_conc, included_measurements_crea)

                    #included_measurements_prev = [- pre_switch_circ.num_measurements + i for i in included_measurements_prev]

                    #full_meas = included_measurements_conc + included_measurements_prev

                    #for i in full_meas:
                        #current_tar = - switch_circ.num_measurements + i
                        #switch_circ.append("DETECTOR", [stim.target_rec(current_tar)], 
                                    #(coords.real, coords.imag, 0))
                    
        elif label in {"Z-STAB"}:

            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
            new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
                
            if coords in [2 + i + 2j + 1j * i for i in range(0,distance, 2)]:

                old_stab = [new_cord1, new_cord2, new_cord3, new_cord4]
                x_old = '*'.join(f"Z{i}" for i in [q2i[current] for current in old_stab])

                switch_round_conc = f"{x_old} -> {1}"

                (included_measurements_conc,) = pre_switch_circ.solve_flow_measurements([
                stim.Flow(switch_round_conc),
                ])

                #print(included_measurements_conc)

                # Buidling new Diagonal detectors
                #for index in included_measurements_conc:
                    #prev_tar = -1 * current_round - previous_round + index - 1
                    #current_tar = -1 * current_round + index - 1
                    #switch_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                                    #(i2q[q_index].real, i2q[q_index].imag, 0))
                    
            if q2i[coords] in index_h:

                old_stab = [new_cord1, new_cord2, new_cord3, new_cord4]
                new_stab = [new_cord2, new_cord3, new_cord4]
                z_old = '*'.join(f"X{i}" for i in [q2i[current] for current in old_stab])
                z_new = '*'.join(f"X{i}" for i in [q2i[current] for current in new_stab])

                stabs_contraction = f"{z_old} -> {1}"
                stabs_creation = f"{1} -> {z_new}"

                (included_measurements_conc,) = switch_circ.solve_flow_measurements([
                                                stim.Flow(stabs_contraction),
                                                ])
                
                (included_measurements_creat,) = switch_circ.solve_flow_measurements([
                                                stim.Flow(stabs_creation),
                                                ])
                
                #print(included_measurements_conc, included_measurements_creat)

                # Buidling new Diagonal detectors
                # for index_old, index_new in zip(included_measurements_conc, included_measurements_creat):
                #     prev_tar = -1 * current_round - previous_round + index_old - 1
                #     current_tar = -1 * current_round + index_new - 1
                #     switch_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                #                    (i2q[q_index].real, i2q[q_index].imag, 0))
                    
        elif label in {"Z-STAB-BOUND-U-H"}:

            new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
            new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_stab = [new_cord3, new_cord4]
            new_stab_index = [q2i[new_cord3], q2i[new_cord4]]

            # Filter for Boundary who has X and Z checks
            if y_index not in new_stab_index:

                z_new = '*'.join(f"X{i}" for i in [q2i[current] for current in new_stab])
                switch_round_conc = f"{1} -> {z_new}"

                (included_measurements_conc,) = switch_circ.solve_flow_measurements([
                stim.Flow(switch_round_conc),
                ])

                # Buidling new Diagonal detectors
                #for index in included_measurements_conc:
                #    current_tar = -1 * current_round + index - 1
                #    switch_circ.append("DETECTOR", [stim.target_rec(current_tar)], 
                #                    (i2q[q_index].real, i2q[q_index].imag, 0))

            else:
                z_new = '*'.join([f"X{new_stab_index[0]}"] + [f"X{new_stab_index[1]}"])
                switch_round_conc = f"{1} -> {z_new}"

                (included_measurements_conc,) = switch_circ.solve_flow_measurements([
                stim.Flow(switch_round_conc),
                ])

                # Buidling new Diagonal detectors
                #for index in included_measurements_conc:
                    #current_tar = -1 * current_round + index - 1
                    #switch_circ.append("DETECTOR", [stim.target_rec(current_tar)], 
                                    #(i2q[q_index].real, i2q[q_index].imag, 0))
                    
        elif label in {"X-STAB-BOUND-R-H"}:

            new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_stab = [new_cord2, new_cord4]
            new_stab_index = [q2i[new_cord2], q2i[new_cord4]]
                
            x_new = '*'.join(f"Z{i}" for i in [q2i[current] for current in new_stab])

            switch_round_conc = f"{1} -> {x_new}"

            (included_measurements_conc,) = switch_circ.solve_flow_measurements([
            stim.Flow(switch_round_conc),
            ])

            # Buidling new Diagonal detectors
            #for index in included_measurements_conc:
            #    current_tar = -1 * current_round + index - 1
            #    switch_circ.append("DETECTOR", [stim.target_rec(current_tar)], 
            #                    (i2q[q_index].real, i2q[q_index].imag, 0))

