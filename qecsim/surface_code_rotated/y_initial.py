from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch, Context

Coord = complex

__all__ = ["y_initial"]

def y_initial(*, lct : Context, patches: dict[str, Patch], 
            cfg : Config, before_round_depol : float, before_m_flip_prob : float, 
            after_r_flip : float, after_c_depol_prob : float) -> stim.Circuit:
    
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    patch = patches["patch"]

    #-Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    distance = cfg.distance
    init_state = cfg.state_init
    stab_to_data = lct.stab_to_data

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab

    """
    For the Y basis initilization we exclude the upper right data qubit
    -> One weight 4 stabilizer reduced to weight 3 and boundary stabilizer missing
    """

    # Finding Upper right qubit index -> need to look in 2-CX
    y_index = 1 + 1j

    ########################
    # Define Initial Circuit
    ########################

    initial_circuit = stim.Circuit()

    initial_circuit.append("TICK")

    #-------Adding-Before-Round-Depol.-Data------------

    if before_round_depol > 0:
        initial_circuit.append("DEPOLARIZE1", data, before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    initial_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", x_stab_index, after_c_depol_prob)
        
    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #2) CX Operations

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":

            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "1-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():

        #Parallel Implementation of CX
        if order == "2-CX":

            # Removing corner CX
            if coord_pairs[0] != y_index and coord_pairs[1] != y_index:
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("CX", index_pairs)
        
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "2-CX":

                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)
    
    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "3-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            initial_circuit.append("CX", index_pairs)

    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
                
        for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
            if order == "4-CX":
                index_pairs = []
                index_pairs.append(q2i[coord_pairs[1]])
                index_pairs.append(q2i[coord_pairs[0]])
                initial_circuit.append("DEPOLARIZE2", index_pairs, after_c_depol_prob)

    #-------Continue-Circuit------------
    
    initial_circuit.append("TICK")

    #3) Basis/ Measurement
    initial_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", x_stab_index, after_c_depol_prob)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.------------

    if before_m_flip_prob > 0:
        initial_circuit.append("X_ERROR", x_stab_index + z_stab_index, before_m_flip_prob)

    #-------Continue-Circuit------------
    
    initial_circuit.append("MR", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if after_r_flip > 0:
        initial_circuit.append("X_ERROR", x_stab_index + z_stab_index, after_r_flip)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #4) DETECTORS -> Measure only deterministic-Stabilizers!

    """
    As we have split the circuit in equal part x and z basis 
    -> we need to add detectors only to their respective deterministic basis

    ATTENTION:
    -> The diagonal where RZ and RX init meet cannot have detectors
    """

    #-----------Filter-stabs-------------

    """
    INEFFICIENT -> ALREADY CALCULATED IN RESET_CIRC
    """

    data_rx : list[complex] = []
    data_rz : list[complex] = []

    xs = [i2q[i].real for i in data]
    ys = [i2q[i].imag for i in data]

    # Calc threshold for diagonal cut        
    s0 = (min(xs)+max(xs))/2 + (min(ys)+max(ys))/2

    skip_coord = 1 + 1j

    # Filtering out Data qubits
    for data_index in data:
        c = i2q[data_index]
        if c == skip_coord:
            continue
            
        # Diagonal Cut at s0 +1 -> filter out
        if (c.real + c.imag) >= s0 + 1:
            data_rx.append(c)
        else:
            data_rz.append(c)

    # Filter Stabs -> Every X stabilizer with data_rx valid and z with only data_z 
    # -> if stabilizer with data_rx & data_rz -> invalid

    stab_rx : list= []
    stab_rz : list = []

    # Determine the neighbouring Data qubits to each Stabilizer
    OFFSETS = {
        "Z-STAB":           [(-1, -1), (+1, -1), (-1, +1), (+1, +1)],
        "Z-STAB-BOUND-L":   [(+1, -1), (+1, +1)],
        "X-STAB-BOUND-R":   [(-1, -1), (-1, +1)],
        "X-STAB":           [(-1, -1), (+1, -1), (-1, +1), (+1, +1)],
        "Z-STAB-BOUND-U":   [(-1, +1), (+1, +1)],
        "X-STAB-BOUND-B":   [(-1, -1), (+1, -1)],
    }

    # Unified detector construction
    for q, qtype in patch.coords.items():

        if qtype in {"Z-STAB", "X-STAB-BOUND-L", "Z-STAB-BOUND-R", "X-STAB", "X-STAB-BOUND-U", "Z-STAB-BOUND-B"}:
            
            #Needed Data Qubits
            neighbor_coords = [(q.real + dx) + (q.imag + dy) * 1j for dx, dy in OFFSETS[qtype]]
            # keep only those neighbors that are actual data qubits
            nb = [nc for nc in neighbor_coords if (nc in data_rx or nc in data_rz)]

            if nb and all(n in data_rx for n in nb):
                stab_rx.append(q2i[q])
            elif nb and all(n in data_rz for n in nb):
                stab_rz.append(q2i[q])

    #-Determine-the-postion-in-the-current-measurement-record-

    num_measurements_initial = len(z_stab_index)
        
    for index, q_index in enumerate(z_stab_index):

        if q_index in stab_rz:

            current_tar = -1 * num_measurements_initial + index
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    num_measurements_initial = len(x_stab_index + z_stab_index)
        
    for index, q_index in enumerate(x_stab_index):

        if q_index in stab_rx:

            current_tar = -1 * num_measurements_initial + index
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))

    return(initial_circuit)
