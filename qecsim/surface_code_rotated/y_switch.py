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
        
        if role == "X-STAB":

            # skip stabs on the Y-cut
            if q_index in index_nh:

                # find this stabilizers position in the previous block
                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                
                
            # Everything but the off diagonal
            elif i2q[q_index] not in [2 + i + 1j * i for i in range(2,distance * 2 + 1, 2)]:

                bottom_neighbour = coord.real + coord.imag*1j + 2j
                bottom_index = q2i[bottom_neighbour]

                # find this stabilizers position in the previous block (Z stab ancilla nex to it)
                for curr_idx, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
                    if q_index == bottom_index:
                        prev_tar = - current_round - previous_round + curr_idx
                
                current_tar = - current_round + index
                switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                
            # Off Diagonal


        elif role == "Z-STAB":

            diagonal_indices = [2 + i + 2j + 1j * i for i in range(0,distance * 2 - 2, 2)]

            #skip stabs on the Y-cut
            if q_index in index_nh:

                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                
                
            # Skip Diagonal
            elif i2q[q_index] not in diagonal_indices:

                left_neighbour = coord.real - 2 + coord.imag*1j
                left_index = q2i[left_neighbour]

                # find this stabilizers position in the previous block (Z stab ancilla nex to it)
                for curr_idx, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
                    if q_index == left_index:
                        prev_tar = - current_round - previous_round + curr_idx
                
                current_tar = - current_round + index
                switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                
            # Diagonal (Except most bottom entry and first)


        elif role in {"X-STAB-BOUND-B", "Z-STAB-BOUND-L"}:
            prev_tar = - current_round - previous_round + index
            current_tar = - current_round + index
            switch_circ.append("DETECTOR",
                [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                (coord.real, coord.imag, 0))
            
        elif role in {"X-STAB-BOUND-R-H"}:
    
            """
            This Detector goes to the left into the Z stabs
            """

            bottom_imag = [i for i in range(4, distance * 2, 4)][-1]
            bottom_real = distance * 2

            btm_idx = bottom_real + bottom_imag * 1j

            if i2q[q_index] != btm_idx:

                left_neighbour = coord.real - 2 + coord.imag*1j
                left_index = q2i[left_neighbour]

                # find this stabilizers position in the previous block (Z stab ancilla nex to it)
                for curr_idx, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
                    if q_index == left_index:
                        prev_tar = - current_round - previous_round + curr_idx
                    
                current_tar = - current_round + index
                switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                
            # Bottom boundary weight 3 check

        elif role in {"X-STAB-BOUND-R"}:
    
            """
            This Detector gets from weight 2 to weight 3 and includes the measurement from the X-STAB-BOUND-R-H
            """


        elif role in {"Z-STAB-BOUND-U-H"}:

            far_right_real = [i for i in range(2, distance * 2, 4)][0]
            far_right_imag = 0

            fr_index = far_right_real + far_right_imag * 1j

            # Every H detector compares current M measurement with M on the left (Normal Upper-Boundary)
            if i2q[q_index] != fr_index:

                bottom_neighbour = coord.real + coord.imag*1j + 2j
                bottom_index = q2i[bottom_neighbour]

                # find this stabilizers position in the previous block
                for curr_idx, q_index in enumerate(x_stab_index + z_stab_index + r_h_stabs + u_h_stabs):
                    if q_index == bottom_index:
                        prev_tar = - current_round - previous_round + curr_idx

                current_tar = - current_round + index
                switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))


        elif role in {"Z-STAB-BOUND-U"}:

            far_left_real = [i for i in range(4, distance * 2, 4)][-1]
            far_left_imag = 0

            fl_index = far_left_real + far_left_imag * 1j

            # Far right detector gets treated normally
            if i2q[q_index] == fl_index:

                # find this stabilizers position in the previous block
                prev_tar = - current_round - previous_round + index
                current_tar = - current_round + index
                switch_circ.append("DETECTOR",
                    [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
                    (coord.real, coord.imag, 0))
                

    full_switch = pre_switch_circ + switch_circ

    return CircuitResult(circuit=full_switch)
