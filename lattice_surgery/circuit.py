import stim
from typing import Dict, Tuple

from .geometry import build_lattice
from .stabilizers import populate_stab_to_data

Coord = complex
Label = str

__all__ = ["surgery_circuit"]

# -------------------------
# Helper Functions
# -------------------------

def _add_boundary_labels(distance: int,
                         ancilla: Dict[Coord, Label],
                         target: Dict[Coord, Label],
                         control: Dict[Coord, Label],
                         surgery: Dict[Coord, Label]) -> None:
    
    """
    Adds the neseccary Boundary and Surgery Stabilizers needed
    """

    max_coord = 2 * distance

    # Z-boundary stabilizers
    for y in range(2, max_coord, 4):
        coord_ancilla = complex(0, y)
        coord_target = complex(distance * 2, y + 2)
        coord_control = complex(0, y + 2 + (distance * 2))
        coord_lattice = complex(y + 2, max_coord)
        ancilla[coord_ancilla] = "Z-STAB-BOUND-L-A"
        target[coord_target] = "Z-STAB-BOUND-L-T"
        control[coord_control] = "Z-STAB-BOUND-L-C"
        surgery[coord_lattice] = "Z-STAB-SURGERY-M"

    for y in range(4, max_coord, 4):
        coord_ancilla = complex(max_coord, y)
        coord_target = complex(max_coord + (distance * 2), y - 2)
        coord_control = complex(max_coord, y - 2 + (distance * 2)) 
        ancilla[coord_ancilla] = "Z-STAB-BOUND-R-A"
        target[coord_target] = "Z-STAB-BOUND-R-T"
        control[coord_control] = "Z-STAB-BOUND-R-C"

    # X-boundary stabilizers
    for y in range(4, max_coord, 4):
        coord_ancilla = complex(y, 0)
        coord_target = complex(y - 2 + (distance * 2), 0)
        coord_control = complex(y - 2, distance * 2)
        coord_lattice = complex(max_coord, y-2)
        ancilla[coord_ancilla] = "X-STAB-BOUND-A-A"
        target[coord_target] = "X-STAB-BOUND-A-T"
        control[coord_control] = "X-STAB-BOUND-A-C"
        surgery[coord_lattice] = "X-STAB-SURGERY-M"

    for y in range(2, max_coord, 4):
        coord_ancilla = complex(y, max_coord)
        coord_target = complex(y + 2 + (distance * 2), max_coord)
        coord_control = complex(y + 2, max_coord + (distance * 2))
        ancilla[coord_ancilla] = "X-STAB-BOUND-B-A"
        target[coord_target] = "X-STAB-BOUND-B-T"
        control[coord_control] = "X-STAB-BOUND-B-C"

    #Adding X Surgery Stabilizer Between Ancilla & Target
    surgery[complex(max_coord, max_coord)] = "X-STAB-SURGERY-B"

    #Adding Z Surgery Stabilizer Between Ancilla & Control
    surgery[complex(0, max_coord)] = "Z-STAB-SURGERY-L"

# -------------------------
# Public function
# -------------------------

def surgery_circuit(distance: int, *, target_state_init: str, control_state_init: str) -> stim.Circuit:
    """
    Returns the full lattice surgery circuit

    Arguments:
                -> target_state_init: In which basis should the target lattice be initlized?
                -> control_state_init: In which basis should the control lattice be initilized?

    Returns:
                -> Fully implemented CX-Gate in stim.Circuit format
    """

    # ------------------------------------------------------------------
    # 1. Build independent square patches (using geometry function)
    # ------------------------------------------------------------------
    qubit_coords_ancilla: Dict[Coord, Label] = build_lattice(distance, offset=0+0j, starting_stabilizer_x=True)
    qubit_coords_target: Dict[Coord, Label] = build_lattice(distance, offset=(distance*2)+0j, starting_stabilizer_x=False)
    qubit_coords_control: Dict[Coord, Label] = build_lattice(distance, offset=0+distance*2j, starting_stabilizer_x=False)
    qubit_coords_surgery: Dict[Coord, Label] = {}

    # ------------------------------------------------------------------
    # 2. Insert boundary & surgery labels (Only get activated in splitting/merging process)
    # ------------------------------------------------------------------
    _add_boundary_labels(distance, qubit_coords_ancilla, qubit_coords_target, qubit_coords_control, qubit_coords_surgery)

    # Merge into one big map
    merged: Dict[Coord, Label] = {
        **qubit_coords_ancilla,
        **qubit_coords_control,
        **qubit_coords_target,
        **qubit_coords_surgery
    }

    # ------------------------------------------------------------------
    # 3. Adding the Mapping from Stabilizer to Data for later CX gate implementation
    # ------------------------------------------------------------------
    stab_to_data: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(merged)

    ancilla_set = set(qubit_coords_ancilla)
    target_set  = set(qubit_coords_target)
    control_set = set(qubit_coords_control)

    # helper to filter the big dict
    get_view = lambda region: {
        pair: order
        for pair, order in stab_to_data.items()
        if pair[0] in region or pair[1] in region
    }

    stab_to_data_ancilla = get_view(ancilla_set)
    stab_to_data_target  = get_view(target_set)
    stab_to_data_control = get_view(control_set)

    # ------------------------------------------------------------------
    # 4. Indexing All Qubits From given Coordinates
    # ------------------------------------------------------------------

    #Indexing Qubits
    q2i: dict[complex, int] = {q: i for i, q in enumerate(
    sorted(merged, key=lambda v: (v.real, v.imag))
    )}

    #Reverse Indexing
    i2q: dict[int, complex] = {i: q for q, i in q2i.items()}

    #Indexing Z and X Stabilizers
    x_stab_index_ancilla = [q2i[q] for q, qtype in qubit_coords_ancilla.items() if qtype in {"X-STAB","X-STAB-BOUND-A-A", "X-STAB-BOUND-B-A"}]
    z_stab_index_ancilla = [q2i[q] for q, qtype in qubit_coords_ancilla.items() if qtype in {"Z-STAB","Z-STAB-BOUND-L-A", "Z-STAB-BOUND-R-A"}]
    x_stab_boundary_b_index_ancilla = [q2i[q] for q, qtype in qubit_coords_ancilla.items() if qtype == "X-STAB-BOUND-B-A"]
    x_stab_index_target = [q2i[q] for q, qtype in qubit_coords_target.items() if qtype in {"X-STAB","X-STAB-BOUND-A-T", "X-STAB-BOUND-B-T"}]
    z_stab_index_target = [q2i[q] for q, qtype in qubit_coords_target.items() if qtype in {"Z-STAB","Z-STAB-BOUND-L-T", "Z-STAB-BOUND-R-T"}]
    x_stab_index_control = [q2i[q] for q, qtype in qubit_coords_control.items() if qtype in {"X-STAB","X-STAB-BOUND-A-C", "X-STAB-BOUND-B-C"}]
    z_stab_index_control = [q2i[q] for q, qtype in qubit_coords_control.items() if qtype in {"Z-STAB","Z-STAB-BOUND-L-C", "Z-STAB-BOUND-R-C"}]

    #Indexing Data-Qubits
    data_ancilla = [q2i[q] for q, qtype in qubit_coords_ancilla.items() if qtype == "DATA"]
    data_target = [q2i[q] for q, qtype in qubit_coords_target.items() if qtype == "DATA"]
    data_control = [q2i[q] for q, qtype in qubit_coords_control.items() if qtype == "DATA"]
    
    #Building Circuit
    inital_circuit = stim.Circuit()
    merge_init_circuit = stim.Circuit()
    merge_round_circuit = stim.Circuit()

    # ------------------------------------------------------------------
    # 5. Building Initilization Circuit
    # ------------------------------------------------------------------

    #Appending Coords
    for q, i in q2i.items():
        inital_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

    """
    Looking at every state preperation seperatly seems to be inefficient
    ->  If not all Operators only once used one gets an incorrect formatting in the 
        timeslice view because of the Operations being in different timeslices in each TICK!
    """

    #Inilizing Ancilla in Plus (Reset) and Control/ Target in desired State
    init_patterns = {
    ("Z0", "Z0"): [("RX", data_ancilla), ("R", data_control + data_target), ("Z", data_ancilla)],
    ("Z0", "Z1"): [("RX", data_ancilla), ("R", data_control + data_target), ("X", data_target), ("Z", data_ancilla)],
    ("Z0", "X+"): [("RX", data_ancilla + data_target), ("R", data_control), ("Z", data_ancilla + data_target)],
    ("Z0", "X-"): [("RX", data_ancilla + data_target), ("R", data_control), ("Z", data_ancilla)],
    ("Z1", "Z0"): [("RX", data_ancilla), ("R", data_control + data_target), ("X", data_control), ("Z", data_ancilla)],
    ("Z1", "Z1"): [("RX", data_ancilla), ("R", data_control + data_target), ("X", data_control + data_target), ("Z", data_ancilla)],
    ("Z1", "X+"): [("RX", data_ancilla + data_target), ("R", data_control), ("X", data_control), ("Z", data_ancilla + data_target)],
    ("Z1", "X-"): [("RX", data_ancilla + data_target), ("R", data_control), ("X", data_control), ("Z", data_ancilla)],
    ("X+", "Z0"): [("RX", data_ancilla + data_control), ("R", data_target), ("Z", data_ancilla + data_control)],
    ("X+", "Z1"): [("RX", data_ancilla + data_control), ("R", data_target), ("X", data_target), ("Z", data_ancilla + data_control)],
    ("X+", "X+"): [("RX", data_ancilla + data_control + data_target), ("Z", data_ancilla + data_control + data_target)],
    ("X+", "X-"): [("RX", data_ancilla + data_control + data_target), ("Z", data_ancilla + data_target)],
    ("X-", "Z0"): [("RX", data_ancilla + data_control), ("R", data_target), ("Z", data_ancilla)],
    ("X-", "Z1"): [("RX", data_ancilla + data_control), ("R", data_target), ("X", data_target), ("Z", data_ancilla)],
    ("X-", "X+"): [("RX", data_ancilla + data_control + data_target), ("Z", data_ancilla + data_control)],
    ("X-", "X-"): [("RX", data_ancilla + data_control + data_target), ("Z", data_ancilla + data_control + data_target)],
    }

    # Apply the initialization pattern
    key = (control_state_init, target_state_init)

    if key not in init_patterns:
        raise ValueError(f"Invalid basis combination: {key}")

    for gate, qubits in init_patterns[key]:
        inital_circuit.append(gate, qubits)


    inital_circuit.append("TICK")

    #Adding h gate for X stabilizers -> Filtering out double coords
    combined_x_stab : list = []
    for coords in (x_stab_index_ancilla + x_stab_index_control + x_stab_index_target):
        if coords not in combined_x_stab:
            combined_x_stab.append(coords)

    inital_circuit.append("H", combined_x_stab)
    inital_circuit.append("TICK")

    #2) CX Operations

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)

    inital_circuit.append("TICK")
            
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)

    inital_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)
    
    inital_circuit.append("TICK")
        
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!):
    inital_circuit.append("TICK")
    inital_circuit.append("H", x_stab_index_ancilla)
    inital_circuit.append("TICK")
    inital_circuit.append("MR", x_stab_index_ancilla + z_stab_index_ancilla)
    inital_circuit.append("TICK")
    inital_circuit.append("H", x_stab_boundary_b_index_ancilla)
    inital_circuit.append("TICK")

    #-----Implementing Detectors for Ancilla (+ State -> X Basis is deterministic)-----

    #Determining Position in the measurement Run of only the Ancilla
    pos_to_index_ancilla_x : list = []
    pos_to_index_ancilla_z : list = []

    for pos, index in enumerate(x_stab_index_ancilla + z_stab_index_ancilla):
        if index in x_stab_index_ancilla:
            pos_to_index_ancilla_x.append([pos, index])

        elif index in z_stab_index_ancilla:
            pos_to_index_ancilla_z.append([pos, index])

    #Adding the needed Detectors
    for index_pos in pos_to_index_ancilla_x:
        current_tar = index_pos[0] - len(x_stab_index_ancilla + z_stab_index_ancilla)
        q_index = index_pos[1]
        inital_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))


    #Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "5-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)

    inital_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
   
        #Parallel Implementation of CX
        if order == "6-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            inital_circuit.append("CX", index_pairs)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    inital_circuit.append("TICK")
    inital_circuit.append("H", x_stab_index_control +  x_stab_index_target)
    inital_circuit.append("TICK")
    inital_circuit.append("MR", control_target_stabs)

    #4) DETECTORS -> Only record Targets from deterministic results. i.e stabilizers from the current basis!
    
    #Determining Postion in the measurement Run of Target & Control
    pos_to_index_control_x : list = []
    pos_to_index_control_z : list = []
    pos_to_index_target_x : list = []
    pos_to_index_target_z : list = []

    for pos, index in enumerate(control_target_stabs):
        if index in x_stab_index_control:
            pos_to_index_control_x.append([pos, index])

        elif index in z_stab_index_control:
            pos_to_index_control_z.append([pos, index])

        elif index in x_stab_index_target:
            pos_to_index_target_x.append([pos, index])

        elif index in z_stab_index_target:
            pos_to_index_target_z.append([pos, index])
    
    #-----Implementing Detectors for Control-----

    #Z-Basis (0/1 - state)
    if control_state_init in {"Z0", "Z1"}:

        for index_pos in pos_to_index_control_z:
            current_tar = index_pos[0] - len(control_target_stabs)
            q_index = index_pos[1]
            inital_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    
    #X-Basis (+/- - state)
    elif control_state_init in {"X-", "X+"}:

        for index_pos in pos_to_index_control_x:
            current_tar = index_pos[0] - len(control_target_stabs)
            q_index = index_pos[1]
            inital_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")

    #-----Implementing Detectors for Target-----
    
    #Z-Basis (0/1 - state)
    if target_state_init in {"Z0", "Z1"}:
  
        for index_pos in pos_to_index_target_z:
            current_tar = index_pos[0] - len(control_target_stabs)
            q_index = index_pos[1]
            inital_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
    
    #X-Basis (+/- - state)
    elif target_state_init in {"X-", "X+"}:

        for index_pos in pos_to_index_target_x:
            current_tar = index_pos[0] - len(control_target_stabs)
            q_index = index_pos[1]
            inital_circuit.append("DETECTOR", [stim.target_rec(current_tar)], (i2q[q_index].real, i2q[q_index].imag, 0))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Target Lattice")

    # ------------------------------------------------------------------
    # 4. Building Merging Ancilla Control Circuit
    # ------------------------------------------------------------------

    #Indexing of the additional Stabilizers included in the merging/splitting process
    z_stab_index_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "Z-STAB-SURGERY-M"]
    z_stab_boundary_l_surgery = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "Z-STAB-SURGERY-L"]
    
    stab_to_data_ancilla_control: dict[list[complex], str] = {}

    """
    We now redefine the CX implementation which now does Control + Ancilla as one lattice!
    """

    #------CX-GATES-ANCILLA-&-CONTROL------
    for coords,string in (qubit_coords_ancilla | qubit_coords_surgery | qubit_coords_control).items():
        #Already in Correct Orientation for Measurement of CX
        if string in {"X-STAB", "X-STAB-BOUND-A-C"}:
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord3 = (coords.real + 1) + (coords.imag + 1) * 1j
            new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
            stab_to_data_ancilla_control[new_cord1, coords] = "1-CX"
            stab_to_data_ancilla_control[new_cord2, coords] = "2-CX"
            stab_to_data_ancilla_control[new_cord3, coords] = "3-CX"
            stab_to_data_ancilla_control[new_cord4, coords] = "4-CX"

        elif string in {"Z-STAB", "Z-STAB-SURGERY-M"}:
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1) + (coords.imag + 1) * 1j
            new_cord3 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord4 = (coords.real - 1) + (coords.imag + 1) * 1j
            stab_to_data_ancilla_control[coords, new_cord1] = "1-CX"
            stab_to_data_ancilla_control[coords, new_cord2] = "2-CX"
            stab_to_data_ancilla_control[coords, new_cord3] = "3-CX"
            stab_to_data_ancilla_control[coords, new_cord4] = "4-CX"

        elif string in {"Z-STAB-BOUND-L-A", "Z-STAB-BOUND-L-C", "Z-STAB-SURGERY-L"}:
            new_cord1 = (coords.real + 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data_ancilla_control[coords, new_cord1] = "1-CX"
            stab_to_data_ancilla_control[coords, new_cord2] = "2-CX"

        elif string in {"Z-STAB-BOUND-R-A", "Z-STAB-BOUND-R-C"}:
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real - 1 ) + (coords.imag + 1) * 1j
            stab_to_data_ancilla_control[coords, new_cord1] = "3-CX"
            stab_to_data_ancilla_control[coords, new_cord2] = "4-CX"
        
        elif string == "X-STAB-BOUND-A-A":
            new_cord1 = (coords.real - 1) + (coords.imag + 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag + 1) * 1j
            stab_to_data_ancilla_control[new_cord1, coords] = "4-CX"
            stab_to_data_ancilla_control[new_cord2, coords] = "3-CX"

        elif string == "X-STAB-BOUND-B-C":
            new_cord1 = (coords.real - 1) + (coords.imag - 1) * 1j
            new_cord2 = (coords.real + 1 ) + (coords.imag - 1) * 1j
            stab_to_data_ancilla_control[new_cord1, coords] = "2-CX"
            stab_to_data_ancilla_control[new_cord2, coords] = "1-CX"

    #Adding h gate for X stabilizers -> Filtering out double coords
    combined_x_stab : list = []
    for coords in (x_stab_index_ancilla + x_stab_index_control + x_stab_index_target):
        if coords not in combined_x_stab:
            combined_x_stab.append(coords)

    merge_init_circuit.append("TICK")
    merge_init_circuit.append("H", combined_x_stab)
    merge_init_circuit.append("TICK")

    #2) CX Operations

    for coord_pairs, order in (stab_to_data_ancilla_control |stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    merge_init_circuit.append("TICK")
            
    for coord_pairs, order in (stab_to_data_ancilla_control |stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    merge_init_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_ancilla_control |stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)
    
    merge_init_circuit.append("TICK")
        
    for coord_pairs, order in (stab_to_data_ancilla_control |stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    #Adding h gate for X stabilizers -> Filtering out double coords in big lattice
    combined_x_stab_ac : list = []
    for coords in (x_stab_index_ancilla + x_stab_index_control):
        if coords not in combined_x_stab_ac:
            combined_x_stab_ac.append(coords)

    combined_z_stab_ac : list = []
    for coords in (z_stab_index_ancilla + z_stab_index_control + z_stab_boundary_l_surgery + z_stab_index_surgery):
        if coords not in combined_z_stab_ac:
            combined_z_stab_ac.append(coords)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement):
    merge_init_circuit.append("TICK")
    merge_init_circuit.append("H", combined_x_stab_ac)
    merge_init_circuit.append("TICK")
    merge_init_circuit.append("MR", combined_z_stab_ac + combined_x_stab_ac)
    merge_init_circuit.append("TICK")

    #Adding Detectors -> Firstly Stabilizers which measurement is already known i.e. outside of merging region

    #Determining Position in the measurement Run of only the Ancilla (Shared Stabilizers excluded)
    pos_to_index_ancilla_x_ac : list = []
    pos_to_index_ancilla_z_ac : list = []

    for pos, index in enumerate(combined_z_stab_ac + combined_x_stab_ac):
        if index in x_stab_index_ancilla:
            if index not in x_stab_boundary_b_index_ancilla:
                pos_to_index_ancilla_x_ac.append([pos, index])

        elif index in z_stab_index_ancilla:
            pos_to_index_ancilla_z_ac.append([pos, index])

    #Adding the needed Detectors
    for index_pos_ac in pos_to_index_ancilla_x_ac:
        current_tar = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac)
        q_index = index_pos_ac[1]
        for index_pos in pos_to_index_ancilla_x:
            if q_index == index_pos[1]:
                previous_target = index_pos[0] - len(combined_z_stab_ac + combined_x_stab_ac) - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
                merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #Determining Position in the measurement Run of only the Control (Shared Stabilizers excluded)
    pos_to_index_control_x_ac : list = []
    pos_to_index_control_z_ac : list = []

    for pos, index in enumerate(combined_z_stab_ac + combined_x_stab_ac):
        if index in x_stab_index_control:
            #Exclude Shared Ancillas i.e. below stabilizers from Ancilla
            if index not in x_stab_boundary_b_index_ancilla:
                pos_to_index_control_x_ac.append([pos, index])

        elif index in z_stab_index_control:
            pos_to_index_control_z_ac.append([pos, index])

    #Adding the needed Detectors dependent of logical state of the lattice

    #Z-Basis (0/1 - state)
    if control_state_init in {"Z0", "Z1"}:

        for index_pos_ac in pos_to_index_control_z_ac:
            current_tar = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac)
            q_index = index_pos_ac[1]
            for index_pos in pos_to_index_control_z:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(combined_z_stab_ac + combined_x_stab_ac) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #X-Basis (+/- - state)
    elif control_state_init in {"X-", "X+"}:

        for index_pos_ac in pos_to_index_control_x_ac:
            current_tar = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac)
            q_index = index_pos_ac[1]
            for index_pos in pos_to_index_control_x:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(combined_z_stab_ac + combined_x_stab_ac) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")
    
    """
    We now initlize the stabilizers newly introduce through the merging; one has to distinguish two types
    1) Boundary stabilizers -> Newly formed weight 4 x stabilizers take measurement record of both old boundary weight 2 stabilizers of control and ancilla
    2) No record history -> Newly formed weight 4 z stabilizers have no old record history and therefore have an non deterministic outcome!
    """
    
    #Adding shared Stabilizers
    pos_to_index_shared_x_stabs : list = []

    for pos, index in enumerate(combined_z_stab_ac + combined_x_stab_ac):
        if index in x_stab_index_control:
            #Include Shared Ancillas i.e. below stabilizers from Ancilla
            if index in x_stab_boundary_b_index_ancilla:
                pos_to_index_shared_x_stabs.append([pos, index])

    #Adding the needed Detectors
    for index_pos_ac in pos_to_index_shared_x_stabs:
        current_tar = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac)
        q_index = index_pos_ac[1]
        for index_pos_a in pos_to_index_ancilla_x:
            if q_index == index_pos_a[1]:
                previous_target_ancilla = index_pos_a[0] - len(combined_z_stab_ac + combined_x_stab_ac) - len(control_target_stabs) - len(x_stab_index_ancilla + z_stab_index_ancilla)
                for index_pos_c in pos_to_index_control_x:
                    if q_index == index_pos_c[1]:
                        previous_target_control = index_pos_c[0] - len(combined_z_stab_ac + combined_x_stab_ac) - len(control_target_stabs) 
                        merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target_ancilla), stim.target_rec(previous_target_control)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #Continue CX-Implementation for Target (As AC-Lattice already has a full run)
    for coord_pairs, order in (stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "5-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    merge_init_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "6-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_init_circuit.append("CX", index_pairs)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    merge_init_circuit.append("TICK")
    merge_init_circuit.append("H", x_stab_index_target)
    merge_init_circuit.append("TICK")
    merge_init_circuit.append("MR", x_stab_index_target + z_stab_index_target)

    #Determining Position in the measurement Run of only the Target (Excluded from merge -> Normal stabilizer measurement)
    pos_to_index_target_x_ac : list = []
    pos_to_index_target_z_ac : list = []

    for pos, index in enumerate(x_stab_index_target + z_stab_index_target):
        if index in x_stab_index_target:
            pos_to_index_target_x_ac.append([pos, index])

        elif index in z_stab_index_target:
            pos_to_index_target_z_ac.append([pos, index])

    #Adding the needed Detectors dependent of logical state of the lattice

    #Z-Basis (0/1 - state)
    if target_state_init in {"Z0", "Z1"}:

        for index_pos_ac in pos_to_index_target_z_ac:
            current_tar = index_pos_ac[0] - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            for index_pos in pos_to_index_target_z:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(x_stab_index_target + z_stab_index_target) - len(combined_z_stab_ac + combined_x_stab_ac) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))

    #X-Basis (+/- - state)
    elif target_state_init in {"X-", "X+"}:

        for index_pos_ac in pos_to_index_target_x_ac:
            current_tar = index_pos_ac[0] - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            for index_pos in pos_to_index_target_x:
                if q_index == index_pos[1]:
                    previous_target = index_pos[0] - len(x_stab_index_target + z_stab_index_target) - len(combined_z_stab_ac + combined_x_stab_ac) - len(control_target_stabs)
                    merge_init_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 1))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")

    #----Merging-Repeat-Circuit---
    #Adding repeat circuit with all stabilizers defined

    merge_round_circuit.append("TICK")
    merge_round_circuit.append("H", combined_x_stab)
    merge_round_circuit.append("TICK")

    #2) CX Operations

    for coord_pairs, order in (stab_to_data_ancilla_control |stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    merge_round_circuit.append("TICK")
            
    for coord_pairs, order in (stab_to_data_ancilla_control |stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "2-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    merge_round_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_ancilla_control |stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "3-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)
    
    merge_round_circuit.append("TICK")
        
    for coord_pairs, order in (stab_to_data_ancilla_control |stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    #Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement):
    merge_round_circuit.append("TICK")
    merge_round_circuit.append("H", combined_x_stab_ac)
    merge_round_circuit.append("TICK")
    merge_round_circuit.append("MR", combined_z_stab_ac + combined_x_stab_ac)
    merge_round_circuit.append("TICK")

    #Adding Detectors -> Firstly Stabilizers which measurement is already known i.e. outside of merging region

    #Adding the needed Detectors for Ancilla
    for index_pos_ac in pos_to_index_ancilla_x_ac:
        current_tar = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac)
        previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_ac + combined_x_stab_ac) - len(x_stab_index_target + z_stab_index_target)
        q_index = index_pos_ac[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 2))

    #Adding the needed Detectors dependent of logical state of the lattice (Control)

    #Z-Basis (0/1 - state)
    if control_state_init in {"Z0", "Z1"}:

        for index_pos_ac in pos_to_index_control_z_ac:
            current_tar = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac)
            previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_ac + combined_x_stab_ac) - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 2))

    #X-Basis (+/- - state)
    elif control_state_init in {"X-", "X+"}:

        for index_pos_ac in pos_to_index_control_x_ac:
            current_tar = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac)
            previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_ac + combined_x_stab_ac) - len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 2))
        
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")
    
    #Adding the Ancialla Control shared Detectors
    for index_pos_ac in pos_to_index_shared_x_stabs:
        current_tar = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac)
        previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_ac + combined_x_stab_ac) - len(x_stab_index_target + z_stab_index_target)
        q_index = index_pos_ac[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 2))
                        
    #Adding newly generated Stabilizers
    pos_to_index_new_z_stabs : list = []

    for pos, index in enumerate(combined_z_stab_ac + combined_x_stab_ac):
        if index in (z_stab_index_surgery + z_stab_boundary_l_surgery):
            pos_to_index_new_z_stabs.append([pos, index])

    #Adding the needed Detectors (Newly Z generated Stabs)
    for index_pos_ac in pos_to_index_new_z_stabs:
        current_tar = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac)
        previous_target = index_pos_ac[0] - 2 * len(combined_z_stab_ac + combined_x_stab_ac) - len(x_stab_index_target + z_stab_index_target)
        q_index = index_pos_ac[1]
        merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 2))

    #Continue CX-Implementation for Target (As AC-Lattice already has a full run)
    for coord_pairs, order in (stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "5-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    merge_round_circuit.append("TICK")

    for coord_pairs, order in (stab_to_data_target).items():
   
        #Parallel Implementation of CX
        if order == "6-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            merge_round_circuit.append("CX", index_pairs)

    #All Stabilizers from the Target and Control Lattice
    control_target_stabs = x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target

    #Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    merge_round_circuit.append("TICK")
    merge_round_circuit.append("H", x_stab_index_target)
    merge_round_circuit.append("TICK")
    merge_round_circuit.append("MR", x_stab_index_target + z_stab_index_target)

    #Adding the needed Detectors dependent of logical state of the lattice

    #Z-Basis (0/1 - state)
    if target_state_init in {"Z0", "Z1"}:

        for index_pos_ac in pos_to_index_target_z_ac:
            current_tar = index_pos_ac[0] - len(x_stab_index_target + z_stab_index_target)
            previous_target = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac) - 2 * len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 2))

    #X-Basis (+/- - state)
    elif target_state_init in {"X-", "X+"}:

        for index_pos_ac in pos_to_index_target_x_ac:
            current_tar = index_pos_ac[0] - len(x_stab_index_target + z_stab_index_target)
            previous_target = index_pos_ac[0] - len(combined_z_stab_ac + combined_x_stab_ac) - 2 * len(x_stab_index_target + z_stab_index_target)
            q_index = index_pos_ac[1]
            merge_round_circuit.append("DETECTOR", [stim.target_rec(current_tar), stim.target_rec(previous_target)], (i2q[q_index].real, i2q[q_index].imag, 2))
    
    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")
    
    merge_round_circuit.append("SHIFT_COORDS", arg=(0,0,1))
    
    #Implementing Final Measurement Round (All Data qubits)
    merge_final_circuit = stim.Circuit()

    #Detectors for Target

    # -> Defining Data to measurement indexing
    Index_to_rec_data : dict[int, int] = {q: i for i, q in enumerate(reversed(data_target))}
    Index_to_rec_target: dict[int, int] = {q: i for i, q in enumerate(reversed(x_stab_index_target + z_stab_index_target))}

    #Z-Basis (0/1 - state)
    if target_state_init in {"Z0", "Z1"}:

        merge_final_circuit.append("TICK")
        merge_final_circuit.append("M", data_target)

        for q, qtype in qubit_coords_target.items():

            if qtype == "Z-STAB":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = (q.real + 1) + (q.imag - 1) * 1j
                lower_left = (q.real - 1) + (q.imag + 1) * 1j
                lower_right = (q.real + 1) + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_upper_right] - 1,
                                -Index_to_rec_data[index_lower_left] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                stab_index = q2i[q]
                last_record = [- Index_to_rec_target[stab_index] - 1 - len(data_target)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                merge_final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 2))

            elif qtype == "Z-STAB-BOUND-L-T":
                #Needed Data Qubits
                upper_right = (q.real + 1) + (q.imag - 1) * 1j
                lower_right = (q.real + 1) + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_right = q2i[upper_right]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_right] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                stab_index = q2i[q]
                last_record = [- Index_to_rec_target[stab_index] - 1 - len(data_target)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                merge_final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 2))

            elif qtype == "Z-STAB-BOUND-R-T":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                lower_left = (q.real - 1) + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_lower_left = q2i[lower_left]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_lower_left] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                stab_index = q2i[q]
                last_record = [- Index_to_rec_target[stab_index] - 1 - len(data_target)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                merge_final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 2))

    #X-Basis (+/- - state)
    elif target_state_init in {"X-", "X+"}:

        #Basis Change
        merge_final_circuit.append("TICK")
        merge_final_circuit.append("H", data_target)
        merge_final_circuit.append("TICK")
        merge_final_circuit.append("M", data_target)

        for q, qtype in qubit_coords_target.items():

            if qtype == "X-STAB":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = (q.real + 1) + (q.imag - 1) * 1j
                lower_left = (q.real - 1) + (q.imag + 1) * 1j
                lower_right = (q.real + 1) + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_upper_right] - 1,
                                -Index_to_rec_data[index_lower_left] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                stab_index = q2i[q]
                last_record = [- Index_to_rec_target[stab_index] - 1 - len(data_target)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                merge_final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 2))

            elif qtype == "X-STAB-BOUND-A-T":
                #Needed Data Qubits
                lower_left = (q.real - 1) + (q.imag + 1) * 1j
                lower_right = (q.real + 1) + (q.imag + 1) * 1j

                #Finding the Correct Data Index
                index_lower_left = q2i[lower_left]
                index_lower_right = q2i[lower_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_lower_left] - 1, -Index_to_rec_data[index_lower_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                stab_index = q2i[q]
                last_record = [- Index_to_rec_target[stab_index] - 1 - len(data_target)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                merge_final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 2))

            elif qtype == "X-STAB-BOUND-B-T":
                #Needed Data Qubits
                upper_left = (q.real - 1) + (q.imag - 1) * 1j
                upper_right = (q.real + 1) + (q.imag - 1) * 1j

                #Finding the Correct Data Index
                index_upper_left = q2i[upper_left]
                index_upper_right = q2i[upper_right]

                #Defining the current record targets
                current_record = [-Index_to_rec_data[index_upper_left] - 1, -Index_to_rec_data[index_upper_right] - 1]

                #Defining the last record targets (Normal Detectors from last round)
                stab_index = q2i[q]
                last_record = [- Index_to_rec_target[stab_index] - 1 - len(data_target)]

                #Combining the record targets
                final_record = current_record + last_record
                
                #Appending Detector
                merge_final_circuit.append("DETECTOR", [stim.target_rec(i) for i in final_record], arg = (q.real, q.imag, 2))

    #Detectors for Ancilla (Non Merged)

    #Detectors for Control (Non Merged)

    #Detectors for merged Section
        
    #Adding conditional Z gate on target, if X_L is uneven

    

    # ------------------------------------------------------------------
    # 5. Building Splitting Ancilla Control Circuit
    # ------------------------------------------------------------------

    """
    -> Redefine Logical Operators
    -> After Split d rounds of Stabilizer Measurements for fault tolerance
    """

    # ------------------------------------------------------------------
    # 6. Building Merging Ancilla Target Circuit
    # ------------------------------------------------------------------

    #Indexing of the additional Stabilizers included in the merging process
    x_stab_index_lattice = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-M"]
    x_stab_boundary_b_lattice = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-B"]

    """
    -> After merge logical Operators have to be redefined
    -> After Merge d rounds of Stabilizer Measurements for fault tolerance
    """

    # ------------------------------------------------------------------
    # 7. Building splitting Ancilla Target Circuit
    # ------------------------------------------------------------------

    """
    -> Redefine Logical Operators
    -> After Split d rounds of Stabilizer Measurements for fault tolerance
    """

    # ------------------------------------------------------------------
    # 8. Appending all Circuits
    # ------------------------------------------------------------------

    inital_circuit += merge_init_circuit
    inital_circuit += merge_round_circuit * (distance - 2)
    inital_circuit += merge_final_circuit


    # ------------------------------------------------------------------
    # 9. Retrieving final Circuit with inlined feedback
    # ------------------------------------------------------------------
    """
    We use inlined feedback to track the necessary flips after measurement outcomes clasically, instead of phsically flipping the qubits
    """

    #return_circuit = inital_circuit.with_inlined_feedback()

    return inital_circuit
