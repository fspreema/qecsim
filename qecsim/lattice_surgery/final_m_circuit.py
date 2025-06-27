from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery, LatticeContext
from .stabilizers import populate_stab_to_data

Coord = complex

def final_m(*, lct : LatticeContext, patches: dict[str, Patch_Ancilla, Patch_Control, Patch_Target, Patch_Surgery,], cfg : Config, before_m_flip_prob : float) -> stim.Circuit:

    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    ancilla_patch = patches["ancilla"]
    target_patch = patches["target"]
    control_patch = patches["control"]
    surgery_patch = patches["surgery"]

    #-Retrieving Global Infomration
    distance = cfg.distance
    q2i = lct.q2i
    i2q = lct.i2q
    control_state_init = cfg.control_state_init
    target_state_init = cfg.target_state_init

    #-Retrieving CX Gate Orders
    stab_to_data = lct.stab_to_data
    stab_to_data_surgery_ac = lct.stab_to_data_surgery_ac
    stab_to_data_surgery_at = lct.stab_to_data_surgery_at

    #-Retrieving Lattice Coords
    qubit_coords_ancilla = ancilla_patch.coords
    qubit_coords_control = control_patch.coords
    qubit_coords_target = target_patch.coords
    qubit_coords_surgery = surgery_patch.coords

    ancilla_set = set(qubit_coords_ancilla)
    target_set  = set(qubit_coords_target)
    control_set = set(qubit_coords_control)

    #-Retrieving Data Coords
    data_ancilla = ancilla_patch.data
    data_control = control_patch.data
    data_target = target_patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index_ancilla = ancilla_patch.x_stab
    z_stab_index_ancilla = ancilla_patch.z_stab
    x_stab_boundary_b_index_ancilla = ancilla_patch.x_bdyB
    z_stab_boundary_r_index_ancilla = ancilla_patch.z_bdyR
    x_stab_index_control = control_patch.x_stab
    z_stab_index_control = control_patch.z_stab
    x_stab_index_target = target_patch.x_stab
    z_stab_index_target = target_patch.z_stab

    ##############################
    # Initlize Measurement Circuit
    ##############################

    measure_circuit = stim.Circuit()

    #############################
    # Meassuring all Data Qubits:
    #############################

    measure_circuit.append("TICK")

    if control_state_init in {"X+", "X-"}:
        measure_circuit.append("MX", data_control)

    elif control_state_init in {"Z0", "Z1"}:
        measure_circuit.append("MZ", data_control)

    if target_state_init in {"X+", "X-"}:
        measure_circuit.append("MX", data_target)

    elif target_state_init in {"Z0", "Z1"}:
        measure_circuit.append("MZ", data_target)

    ##############################
    # Defining Logical Observables
    ##############################

    if control_state_init in {"Z0", "Z1"}:

        # Control stabilized by x logical
        log_z_c = []

        for real in range(1, (distance * 2), 2):
            log_z_c.append(q2i[real + ((distance * 2) + 1) * 1j])

        tar_rec = []

        for rec_pos, index in enumerate(data_control + data_target):
            if index in log_z_c:
                tar_rec.append(rec_pos)

        measure_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(- len(data_control + data_target) + k) for k in tar_rec], 0)

    elif control_state_init in {"X+", "X-"}:

        # Control stabilized by x logical
        log_x_ct = []

        for imag in range(((distance * 2) + 1), (distance * 4), 2):
            log_x_ct.append(q2i[1 + imag*1j])

        for imag in range(1, (distance * 2), 2):
            log_x_ct.append(q2i[((distance * 2) + 1) + imag*1j])

        tar_rec = []

        for rec_pos, index in enumerate(data_control + data_target):
            if index in log_x_ct:
                tar_rec.append(rec_pos)

        measure_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(-len(data_control + data_target) + k) for k in tar_rec], 0)

    """
    if target_state_init in {"X+", "X-"}:

        # Control stabilized by x logical
        log_x_c = []

        for imag in range(1, (distance * 2), 2):
            log_x_c.append(q2i[1 + imag*1j])

        #Rewriting in correct form i.e. X1 X2 etc...
        targets_c = [f"X{i}" for i in log_x_c]

        measure_circuit.append("OBSERVABLE_INCLUDE", targets_c, 0)

    elif target_state_init in {"Z0", "Z1"}:

        # Control stabilized by x logical
        log_x_c = []

        for imag in range(1, (distance * 2), 2):
            log_x_c.append(q2i[1 + imag*1j])

        #Rewriting in correct form i.e. X1 X2 etc...
        targets_c = [f"X{i}" for i in log_x_c]

        measure_circuit.append("OBSERVABLE_INCLUDE", targets_c, 0)

    """

    return measure_circuit