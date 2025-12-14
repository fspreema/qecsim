from collections.abc import Mapping

import stim

from src.core.cx_builder import cx_builder
from src.core.data_models import (
    ConfigLatticeSurgery as Config,
    LatticeContext,
    PatchAncilla,
    PatchControl,
    PatchSurgery,
    PatchTarget,
)

Coord = complex


def merge(
    *,
    lct: LatticeContext,
    patches: Mapping[str, PatchAncilla | PatchControl | PatchTarget | PatchSurgery],
    cfg: Config,
    merging_type: str,
) -> stim.Circuit:
    """
    modified_measurement: str
    -> Indicates if XX/ZZ measurements get modified to YY YX ZY etc. joint Measurements
    -> Dependent on the flow selected
    """

    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # -Loading in Patches
    ancilla_patch = patches["ancilla"]
    target_patch = patches["target"]
    control_patch = patches["control"]
    surgery_patch = patches["surgery"]

    # -Retrieving Global Infomration
    distance = cfg.distance
    q2i = lct.q2i

    rounds = distance

    # -Retrieving CX Gate Orders
    stab_to_data = lct.stab_to_data
    stab_to_data_surgery_ac = lct.stab_to_data_surgery_ac
    stab_to_data_surgery_at = lct.stab_to_data_surgery_at

    # -Retrieving Lattice Coords
    qubit_coords_control = control_patch.coords
    qubit_coords_target = target_patch.coords
    qubit_coords_surgery = surgery_patch.coords

    target_set = set(qubit_coords_target)
    control_set = set(qubit_coords_control)

    # helper to filter the big dict
    def get_view(region):
        return {
            pair: order
            for pair, order in stab_to_data.items()
            if pair[0] in region or pair[1] in region
        }

    stab_to_data_target = get_view(target_set)
    stab_to_data_control = get_view(control_set)

    # -Retrieving Index from Stabilizers of the Lattices
    x_stab_index_ancilla = ancilla_patch.x_stab
    z_stab_index_ancilla = ancilla_patch.z_stab
    x_stab_boundary_b_index_ancilla = ancilla_patch.x_bdy_b
    x_stab_index_control = control_patch.x_stab
    z_stab_index_control = control_patch.z_stab
    x_stab_index_target = target_patch.x_stab
    z_stab_index_target = target_patch.z_stab

    # All Stabilizers from the Target and Control Lattice
    control_target_stabs = (
        x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target
    )

    ################################
    # Define initial merging Circuit
    ################################

    merge_init_circuit = stim.Circuit()

    # Indexing of the additional Stabilizers included in the merging/splitting process
    z_stab_index_surgery = [
        q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "Z-STAB-SURGERY-M"
    ]
    z_stab_boundary_l_surgery = [
        q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "Z-STAB-SURGERY-L"
    ]
    x_stab_index_surgery = [
        q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-M"
    ]
    x_stab_boundary_b_surgery = [
        q2i[q] for q, qtype in qubit_coords_surgery.items() if qtype == "X-STAB-SURGERY-B"
    ]

    """
    We now redefine the CX implementation which now does Control + Ancilla or Target + Ancilla 
    as one lattice!
    """

    ###################################################
    # Set general global settings for each merging type
    ###################################################

    if merging_type == "AC":
        stab_to_data_curr_merg = stab_to_data_surgery_ac
        stab_to_data_untouched_circ = stab_to_data_target
        x_stab_index_untouched_circ = x_stab_index_target
        z_stab_index_untouched_circ = z_stab_index_target

    elif merging_type == "AT":
        stab_to_data_curr_merg = stab_to_data_surgery_at
        stab_to_data_untouched_circ = stab_to_data_control
        x_stab_index_untouched_circ = x_stab_index_control
        z_stab_index_untouched_circ = z_stab_index_control

    else:
        raise ValueError("No valid merging Type in Function selected!")

    ############################
    # CX-GATES-ANCILLA-&-CONTROL
    ############################

    combined_x_stab: list = []

    if merging_type == "AC":
        # Adding h gate for X stabilizers on all lattices -> Filtering out double coords
        for coords in x_stab_index_ancilla + x_stab_index_control + x_stab_index_target:
            if coords not in combined_x_stab:
                combined_x_stab.append(coords)

    elif merging_type == "AT":
        # Adding h gate for X stabilizers on all lattices -> Filtering out double coords
        for coords in (
            x_stab_index_ancilla
            + x_stab_index_control
            + x_stab_index_target
            + x_stab_boundary_b_surgery
            + x_stab_index_surgery
        ):
            if coords not in combined_x_stab:
                combined_x_stab.append(coords)

    # Adding reset from initial round and from AC split round
    merge_init_circuit.append("TICK")
    merge_init_circuit.append("R", control_target_stabs)

    merge_init_circuit.append("TICK")
    merge_init_circuit.append("H", combined_x_stab)

    merge_init_circuit.append("TICK")

    combined_x_stab_merging_lattices: list = []
    combined_z_stab_merging_lattices: list = []

    if merging_type == "AC":
        # Adding h gate for X stabilizers only on merging lattices
        # -> Filtering out double coords in big lattice
        for coords in x_stab_index_ancilla + x_stab_index_control:
            if coords not in combined_x_stab_merging_lattices:
                combined_x_stab_merging_lattices.append(coords)

        for coords in (
            z_stab_index_ancilla
            + z_stab_index_control
            + z_stab_boundary_l_surgery
            + z_stab_index_surgery
        ):
            if coords not in combined_z_stab_merging_lattices:
                combined_z_stab_merging_lattices.append(coords)

    elif merging_type == "AT":
        # Adding h gate for X stabilizers only on merging lattices
        # -> Filtering out double coords in big lattice
        for coords in (
            x_stab_index_ancilla
            + x_stab_index_target
            + x_stab_boundary_b_surgery
            + x_stab_index_surgery
        ):
            if coords not in combined_x_stab_merging_lattices:
                combined_x_stab_merging_lattices.append(coords)

        for coords in z_stab_index_ancilla + z_stab_index_target:
            if coords not in combined_z_stab_merging_lattices:
                combined_z_stab_merging_lattices.append(coords)

    # CX Operations
    joined_dict = stab_to_data_curr_merg | stab_to_data_untouched_circ

    cx_builder(
        q2i=q2i,
        stab_to_data=joined_dict,
        circuit=merge_init_circuit,
    )

    # Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement):
    # I.e. return to Z-Basis where needed and Measure
    merge_init_circuit.append("H", combined_x_stab_merging_lattices)
    merge_init_circuit.append("TICK")

    merge_init_circuit.append(
        "M",
        combined_z_stab_merging_lattices + combined_x_stab_merging_lattices,
    )
    merge_init_circuit.append("TICK")

    merge_init_circuit.append(
        "R",
        combined_z_stab_merging_lattices + combined_x_stab_merging_lattices,
    )
    merge_init_circuit.append("TICK")

    #########################################
    # Continue CX-Implementation for untouched lattice (As AC/AT-Lattice already has a full run)
    #########################################

    """
    Adding needed H-Gates for X-Stabs which are shared between merged lattice and untouched lattice
    """

    if merging_type == "AT":
        merge_init_circuit.append("H", x_stab_boundary_b_index_ancilla)

        merge_init_circuit.append("TICK")

    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data_untouched_circ,
        circuit=merge_init_circuit,
        orders=("5-CX", "6-CX"),
    )

    # Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    merge_init_circuit.append("H", x_stab_index_untouched_circ)
    merge_init_circuit.append("TICK")

    merge_init_circuit.append("M", x_stab_index_untouched_circ + z_stab_index_untouched_circ)

    ########################
    # Merging-Repeat-Circuit
    ########################

    # Defining Repeat Circuit
    merge_round_circuit = stim.Circuit()

    # Reinitializing Stabilizers and add basis change where needed
    merge_round_circuit.append("TICK")
    merge_round_circuit.append("R", x_stab_index_untouched_circ + z_stab_index_untouched_circ)

    merge_round_circuit.append("TICK")
    merge_round_circuit.append("H", combined_x_stab)

    merge_round_circuit.append("TICK")

    # CX Operations for Ancilla qubits
    cx_builder(
        q2i=q2i,
        stab_to_data=joined_dict,
        circuit=merge_round_circuit,
    )

    # Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement):
    # I.e. return to Z-Basis where needed and Measure
    merge_round_circuit.append("H", combined_x_stab_merging_lattices)
    merge_round_circuit.append("TICK")

    merge_round_circuit.append(
        "M",
        combined_z_stab_merging_lattices + combined_x_stab_merging_lattices,
    )
    merge_round_circuit.append("TICK")

    merge_round_circuit.append(
        "R",
        combined_z_stab_merging_lattices + combined_x_stab_merging_lattices,
    )
    merge_round_circuit.append("TICK")

    ############################################################################################
    # Continue CX-Implementation for untouched lattice (As AC/AT-Lattice already has a full run)
    ############################################################################################

    """
    Adding needed H-Gates for X-Stabs which are shared between merged lattice and untouched lattice
    """

    if merging_type == "AT":
        merge_round_circuit.append("H", x_stab_boundary_b_index_ancilla)
        merge_round_circuit.append("TICK")

    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data_untouched_circ,
        circuit=merge_round_circuit,
        orders=("5-CX", "6-CX"),
    )

    # Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    merge_round_circuit.append("H", x_stab_index_untouched_circ)
    merge_round_circuit.append("TICK")

    merge_round_circuit.append("M", x_stab_index_untouched_circ + z_stab_index_untouched_circ)

    #####################################################
    # Adding Circuits & Receving the MXX/MZZ Measurements
    #####################################################

    merge_init_circuit += merge_round_circuit * (rounds - 1)

    return merge_init_circuit
