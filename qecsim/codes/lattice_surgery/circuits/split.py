from collections.abc import Mapping

import stim

from qecsim.codes.lattice_surgery.logical_strings import get_logical_strings
from qecsim.core.cx_builder import cx_builder
from qecsim.core.data_models import (
    ConfigLatticeSurgery as Config,
    LatticeContext,
    PatchAncilla,
    PatchControl,
    PatchSurgery,
    PatchTarget,
)

Coord = complex


def split(
    *,
    lct: LatticeContext,
    patches: Mapping[str, PatchAncilla | PatchControl | PatchTarget | PatchSurgery],
    cfg: Config,
    split_type: str,
) -> stim.Circuit:
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

    # -Retrieving Lattice Coords
    qubit_coords_surgery = surgery_patch.coords

    # -Retrieving Data Coords
    data_ancilla = ancilla_patch.data

    # -Retrieving Index from Stabilizers of the Lattices
    x_stab_index_ancilla = ancilla_patch.x_stab
    z_stab_index_ancilla = ancilla_patch.z_stab
    x_stab_boundary_b_index_ancilla = ancilla_patch.x_bdy_b
    x_stab_index_control = control_patch.x_stab
    z_stab_index_control = control_patch.z_stab
    x_stab_index_target = target_patch.x_stab
    z_stab_index_target = target_patch.z_stab

    ###################################################
    # Set general global settings for each merging type
    ###################################################

    if split_type == "AC":
        x_stab_index_untouched_circ = x_stab_index_target
        z_stab_index_untouched_circ = z_stab_index_target

    elif split_type == "AT":
        x_stab_index_untouched_circ = x_stab_index_control
        z_stab_index_untouched_circ = z_stab_index_control

    else:
        raise ValueError("No valid merging Type in Function selected!")

    ################################
    # Define initial split Circuit
    ################################

    split_init_circuit = stim.Circuit()

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

    ##############################
    # Pre Round resets
    ##############################

    # Adding resets from merge ac & at
    split_init_circuit.append("TICK")
    split_init_circuit.append("R", x_stab_index_untouched_circ + z_stab_index_untouched_circ)

    # Adding h gate for X stabilizers -> Filtering out double coords
    combined_x_stab: list = []
    for coords in x_stab_index_ancilla + x_stab_index_control + x_stab_index_target:
        if coords not in combined_x_stab:
            combined_x_stab.append(coords)

    split_init_circuit.append("TICK")

    split_init_circuit.append("H", combined_x_stab)

    ####################################################
    # CX Operations
    ####################################################

    split_init_circuit.append("TICK")

    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=split_init_circuit,
    )

    # Retreive Boundary + Normal Stabilizers Ancilla
    # (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!)
    split_init_circuit.append("H", x_stab_index_ancilla)

    split_init_circuit.append("TICK")

    split_init_circuit.append("M", x_stab_index_ancilla + z_stab_index_ancilla)
    split_init_circuit.append("TICK")
    split_init_circuit.append("R", x_stab_index_ancilla + z_stab_index_ancilla)

    split_init_circuit.append("TICK")
    split_init_circuit.append("H", x_stab_boundary_b_index_ancilla)

    split_init_circuit.append("TICK")

    # Determining Position in the measurement Run of only the Ancilla
    combined_x_stab_merging_lattices: list = []
    combined_z_stab_merging_lattices: list = []

    """
    Only Comparing the inner stabilizers and the boundarys who are not involed in the merging with
    one another
    """

    if split_type == "AT":
        # Adding h gate for X stabilizers only on merging lattices -> Filtering out double coords in
        # big lattice
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

    # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=split_init_circuit,
        orders=("5-CX", "6-CX"),
    )

    # All Stabilizers from the Target and Control Lattice
    control_target_stabs = (
        x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target
    )

    # Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    split_init_circuit.append("H", x_stab_index_control + x_stab_index_target)

    split_init_circuit.append("TICK")

    split_init_circuit.append("M", control_target_stabs)

    ###########################
    # Implementing Repeat Block
    ###########################

    split_repeat_circuit = stim.Circuit()

    # Adding Reset operations from previous round
    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("R", control_target_stabs)

    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("H", combined_x_stab)

    split_repeat_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=split_repeat_circuit,
    )

    # Retreive Boundary + Normal Stabilizers Ancilla:
    split_repeat_circuit.append("H", x_stab_index_ancilla)

    split_repeat_circuit.append("TICK")

    split_repeat_circuit.append("M", x_stab_index_ancilla + z_stab_index_ancilla)
    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("R", x_stab_index_ancilla + z_stab_index_ancilla)

    split_repeat_circuit.append("TICK")
    split_repeat_circuit.append("H", x_stab_boundary_b_index_ancilla)

    split_repeat_circuit.append("TICK")

    ####################################################
    # Implementing Detectors for Ancilla
    ####################################################

    # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=split_repeat_circuit,
        orders=("5-CX", "6-CX"),
    )

    # All Stabilizers from the Target and Control Lattice
    control_target_stabs = (
        x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target
    )

    # Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    split_repeat_circuit.append("H", x_stab_index_control + x_stab_index_target)

    split_repeat_circuit.append("TICK")

    split_repeat_circuit.append("M", control_target_stabs)

    ###########################
    # Adding Final Circ
    ###########################
    """
    In this Section we add the Conditional X_L and Z_L depending on the XX and ZZ Measurements
    """

    split_final_circuit = stim.Circuit()

    # Adding Reset operations from previous round
    split_final_circuit.append("TICK")
    split_final_circuit.append("R", control_target_stabs)

    split_final_circuit.append("TICK")
    split_final_circuit.append("H", combined_x_stab)

    split_final_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=split_final_circuit,
    )

    # Retreive Boundary + Normal Stabilizers Ancilla:
    split_final_circuit.append("H", x_stab_index_ancilla)

    split_final_circuit.append("TICK")

    split_final_circuit.append("M", x_stab_index_ancilla + z_stab_index_ancilla)
    split_final_circuit.append("TICK")
    split_final_circuit.append("R", x_stab_index_ancilla + z_stab_index_ancilla)

    split_final_circuit.append("TICK")
    split_final_circuit.append("H", x_stab_boundary_b_index_ancilla)

    split_final_circuit.append("TICK")

    # -------------------------------------------------------------------------------
    """
    As these are dependent on the measurement outcome and also on the type of merge
    -> Seperated into two Parts conditioned on the Merge/Split type

    The logical ZZ Observable is already defined by the newly implemented Z stabilizers on the merge
    -> Product of the stabilizers give measurement result
    """
    ################################################################################
    # Finding newly generated Stabilizer postion in the Measurement-Rec of the Merge
    ################################################################################

    pos_to_index_newly_gen_stabs: list = []

    if split_type == "AC":
        for pos, index in enumerate(
            combined_z_stab_merging_lattices + combined_x_stab_merging_lattices,
        ):
            if index in (z_stab_index_surgery + z_stab_boundary_l_surgery):
                pos_to_index_newly_gen_stabs.append([pos, index])

    elif split_type == "AT":
        for pos, index in enumerate(
            combined_z_stab_merging_lattices + combined_x_stab_merging_lattices,
        ):
            if index in (x_stab_index_surgery + x_stab_boundary_b_surgery):
                pos_to_index_newly_gen_stabs.append([pos, index])

    logical_obs_rec_tar = []

    inner_record = len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
    skipped_records = len(
        x_stab_index_untouched_circ
        + z_stab_index_untouched_circ
        + rounds * (control_target_stabs + x_stab_index_ancilla + z_stab_index_ancilla),
    )

    for index_pos_merge in pos_to_index_newly_gen_stabs:
        logical_obs_rec_tar.append(index_pos_merge[0] - inner_record - skipped_records)

    #####################################################
    # Defining Logical Data Qubit string for all lattices
    #####################################################

    log_strings = get_logical_strings(q2i, distance)
    a_log_obs_z_index = log_strings["a_z"]
    t_log_obs_z_index = log_strings["t_z"]
    c_log_obs_x_index = log_strings["c_x"]

    ####################################
    # Adding Conditional CZ/CX-Operators
    ####################################

    if split_type == "AC":
        for records in logical_obs_rec_tar:
            for data in c_log_obs_x_index:
                split_final_circuit.append("CX", [stim.target_rec(records), data])

    elif split_type == "AT":
        for records in logical_obs_rec_tar:
            for data in t_log_obs_z_index:
                split_final_circuit.append("CZ", [stim.target_rec(records), data])

        ###################################
        # Meassuring Ancilla in the Z Basis
        ###################################

        # Meassuring Data
        split_final_circuit.append("TICK")
        split_final_circuit.append("MZ", data_ancilla)
        split_final_circuit.append("TICK")

        # Adding the conditional Gate on Control
        for rec_tar, index in enumerate(data_ancilla):
            if index in a_log_obs_z_index:
                for data in c_log_obs_x_index:
                    split_final_circuit.append(
                        "CX",
                        [stim.target_rec(-len(data_ancilla) + rec_tar), data],
                    )

    # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    split_final_circuit.append("TICK")

    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=split_final_circuit,
        orders=("5-CX", "6-CX"),
    )

    # All Stabilizers from the Target and Control Lattice
    control_target_stabs = (
        x_stab_index_control + x_stab_index_target + z_stab_index_control + z_stab_index_target
    )

    # Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    split_final_circuit.append("H", x_stab_index_control + x_stab_index_target)

    split_final_circuit.append("TICK")

    split_final_circuit.append("M", control_target_stabs)

    ##########################################
    # Adding Repeat Circ and returning circuit
    ##########################################

    split_init_circuit += split_repeat_circuit * (rounds - 2)
    split_init_circuit += split_final_circuit

    return split_init_circuit
