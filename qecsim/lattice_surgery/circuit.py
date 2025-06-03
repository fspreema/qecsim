import stim
from typing import Dict, Tuple

from .geometry import build_lattice
from .stabilizers import populate_stab_to_data
from .inital import initial
from .merging import merge

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
        coord_surgery = complex(y + 2, max_coord)
        ancilla[coord_ancilla] = "Z-STAB-BOUND-L-A"
        target[coord_target] = "Z-STAB-BOUND-L-T"
        control[coord_control] = "Z-STAB-BOUND-L-C"
        surgery[coord_surgery] = "Z-STAB-SURGERY-M"

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
        coord_surgery = complex(y - 2, max_coord)
        ancilla[coord_ancilla] = "X-STAB-BOUND-A-A"
        target[coord_target] = "X-STAB-BOUND-A-T"
        control[coord_control] = "X-STAB-BOUND-A-C"
        surgery[coord_surgery] = "X-STAB-SURGERY-M"

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

    ###############################################################
    # 1. Build independent square patches (using geometry function)
    ###############################################################
    qubit_coords_ancilla: Dict[Coord, Label] = build_lattice(distance, offset=0+0j, starting_stabilizer_x=True)
    qubit_coords_target: Dict[Coord, Label] = build_lattice(distance, offset=(distance*2)+0j, starting_stabilizer_x=False)
    qubit_coords_control: Dict[Coord, Label] = build_lattice(distance, offset=0+distance*2j, starting_stabilizer_x=False)
    qubit_coords_surgery: Dict[Coord, Label] = {}

    #######################################################################################
    # 2. Insert boundary & surgery labels (Only get activated in splitting/merging process)
    #######################################################################################
    _add_boundary_labels(distance, qubit_coords_ancilla, qubit_coords_target, qubit_coords_control, qubit_coords_surgery)

    # Merge into one big map
    merged: Dict[Coord, Label] = {
        **qubit_coords_ancilla,
        **qubit_coords_control,
        **qubit_coords_target,
        **qubit_coords_surgery
    }

    ################################################################################
    # 3. Adding the Mapping from Stabilizer to Data for later CX gate implementation
    ################################################################################
    stab_to_data: Dict[Tuple[Coord, Coord], str] = populate_stab_to_data(merged, merging = False)

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

    ###############################################
    # 4. Indexing All Qubits From given Coordinates
    ###############################################

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

    ###################################
    # 5. Building Initilization Circuit
    ###################################

    initial_circuit = initial(data_ancilla = data_ancilla, data_control = data_control, 
            data_target = data_target, q2i = q2i, i2q = i2q,
            control_state_init = control_state_init, target_state_init = target_state_init,
            stab_to_data = stab_to_data,
            x_stab_index_ancilla = x_stab_index_ancilla, z_stab_index_ancilla = z_stab_index_ancilla, 
            x_stab_boundary_b_index_ancilla = x_stab_boundary_b_index_ancilla,
            x_stab_index_control = x_stab_index_control, z_stab_index_control = z_stab_index_control,
            x_stab_index_target = x_stab_index_target, z_stab_index_target = z_stab_index_target)
   
    #############################################
    # 6. Building Merging Ancilla Control Circuit
    #############################################

    merged_circuit_AC = merge(distance = distance, data_ancilla = data_ancilla,  data_target = data_target, data_control = data_control, 
                              q2i = q2i, i2q = i2q, control_state_init = control_state_init, target_state_init = target_state_init,
                              x_stab_index_ancilla = x_stab_index_ancilla, z_stab_index_ancilla = z_stab_index_ancilla, 
                              x_stab_boundary_b_index_ancilla = x_stab_boundary_b_index_ancilla, x_stab_index_control = x_stab_index_control, 
                              z_stab_index_control = z_stab_index_control, x_stab_index_target = x_stab_index_target, 
                              z_stab_index_target = z_stab_index_target, qubit_coords_surgery = qubit_coords_surgery, 
                              qubit_coords_ancilla = qubit_coords_ancilla, qubit_coords_control = qubit_coords_control,
                              stab_to_data_target = stab_to_data_target, qubit_coords_target = qubit_coords_target)

    ###############################################
    # 5. Building Splitting Ancilla Control Circuit
    ###############################################

    """
    -> Redefine Logical Operators
    -> After Split d rounds of Stabilizer Measurements for fault tolerance
    """

    ############################################
    # 6. Building Merging Ancilla Target Circuit
    ############################################

    #Indexing of the additional Stabilizers included in the merging process
    x_stab_index_lattice = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-M"]
    x_stab_boundary_b_lattice = [q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-B"]

    """
    -> After merge logical Operators have to be redefined
    -> After Merge d rounds of Stabilizer Measurements for fault tolerance
    """

    ##############################################
    # 7. Building splitting Ancilla Target Circuit
    ##############################################

    """
    -> Redefine Logical Operators
    -> After Split d rounds of Stabilizer Measurements for fault tolerance
    """

    ###########################
    # 8. Appending all Circuits
    ###########################

    initial_circuit += merged_circuit_AC

    ###################################################
    # 9. Retrieving final Circuit with inlined feedback
    ###################################################

    """
    We use inlined feedback to track the necessary flips after measurement outcomes clasically, instead of phsically flipping the qubits
    """

    #return_circuit = inital_circuit.with_inlined_feedback()

    return initial_circuit
