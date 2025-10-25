import stim

from qecsim.core.cx_builder import cx_builder
from qecsim.core.data_models import (
    CircuitResult,
    ConfigSurface as Config,
    Context,
    NoiseModel,
    Patch,
)

Coord = complex

__all__ = ["y_switch_circ"]

def y_switch_circ(*,
                  lct: Context,
                  patch: dict[str, Patch],
                  offset: complex = 0 + 0j,
                  cfg: Config,
                  noise: NoiseModel) -> CircuitResult:

    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Retrieving Global Infomration
    q2i = lct.q2i
    log_obs = cfg.obs
    distance = cfg.distance
    stab_to_data = lct.stab_to_data
    stab_to_data_switch = lct.stab_to_data_modified
    stab_to_data_switch_xcy = lct.stab_to_data_modified2

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab
    switch_stab_apply_h = patch.stab_switch_apply_h

    # Finding Upper right qubit index -> need to look in 2-CX
    y_coords = 1 + 1j + offset
    y_index = q2i[y_coords]

    ###########################
    # Define Repetition Circuit
    ###########################

    #-----BUILDING-REPETITION-CIRC------

    pre_switch_circ = stim.Circuit()

    pre_switch_circ.append("TICK")
    pre_switch_circ.append("R", x_stab_index + z_stab_index)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    pre_switch_circ.append("TICK")
    pre_switch_circ.append("H", x_stab_index)

    #-------Continue-Circuit------------

    pre_switch_circ.append("TICK")

    #2) CX Operations
    # Adding all the CX gates
    cx_builder(q2i= q2i,
               stab_to_data= stab_to_data,
               circuit= pre_switch_circ,
               excluded_index= y_index,
               noise= noise,
               noise_overwrite= True)

    #-------Continue-Circuit------------

    #3) Basis/ Measurement
    pre_switch_circ.append("H", x_stab_index)

    pre_switch_circ.append("TICK")

    pre_switch_circ.append("MZ", x_stab_index + z_stab_index)

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    pre_switch_circ.append("SHIFT_COORDS", arg = (0,0,1))

    #4) Detectors
    _num_measurements_repeat = len(x_stab_index + z_stab_index)
    """
    for index, q_index in enumerate(x_stab_index + z_stab_index):
        prev_tar = -2 * num_measurements_repeat + index
        current_tar = -1 * num_measurements_repeat + index
        pre_switch_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)], 
                             (i2q[q_index].real, i2q[q_index].imag, 0))
    """
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
            if cords.real - offset.real > cords.imag - offset.imag:
                index_h.append(q2i[cords])

            # Filtering out the X_DAG -> Not on Data
            elif cords.real - offset.real == cords.imag - offset.imag:
                if qtype != "DATA":
                    index_x_deg.append(q2i[cords])

            else:
                index_nh.append(q2i[cords])

    #---------Adding-RY-Gate----------------

    switch_circ.append("RY", y_index)

    #---------Adding-Resets-Old&New-Stabs---

    switch_circ.append("RZ", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)

    #1) Reset/ Basis
    switch_circ.append("TICK")
    switch_circ.append("H", x_stab_index + r_h_stabs)

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

    #-------Continue-Circuit------------

    switch_circ.append("TICK")
    switch_circ.append("M", x_stab_index + z_stab_index + r_h_stabs + u_h_stabs)
    switch_circ.append("TICK")

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    switch_circ.append("SHIFT_COORDS", arg = (0,0,1))

    ############################################################
    # Adding Observable Includes if logical basis differs from Y
    ############################################################

    if log_obs == "Z":

        """
        We need to remove the added Pauli measurement from the end of the circuit
        -> Else the Z paulis tring would anticommute with the RZ reset of the data
        """

        # Getting corresponding logical string and rec
        log_z = []

        for real in range(1, (distance * 2), 2):
            log_z.append(q2i[real + (distance * 2) *1j - 1j])

        # XORing the observable away
        switch_circ.append("OBSERVABLE_INCLUDE", [f"Z{index}" for index in log_z], 0)

    if log_obs == "X":

        """
        We need to remove the added Pauli measurement from the end of the circuit
        -> Else the X paulis tring would anticommute with the RZ reset of the data
        """

        # Getting corresponding logical string and rec
        log_x = []

        for imag in range(1, (distance * 2), 2):
            log_x.append(q2i[distance * 2 - 1 + imag*1j])

        # XORing the observable away
        switch_circ.append("OBSERVABLE_INCLUDE", [f"X{index}" for index in log_x], 0)

    full_switch = pre_switch_circ + switch_circ

    return CircuitResult(circuit=full_switch)
