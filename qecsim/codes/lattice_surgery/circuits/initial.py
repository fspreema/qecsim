from collections.abc import Mapping

import stim

from qecsim.core.cx_builder import cx_builder
from qecsim.core.data_models import (
    ConfigLatticeSurgery as Config,
    LatticeContext,
    NoiseModel,
    PatchAncilla,
    PatchControl,
    PatchSurgery,
    PatchTarget,
)

Coord = complex

__all__ = ["initial"]


def initial(
    *,
    lct: LatticeContext,
    patches: Mapping[str, PatchAncilla | PatchTarget | PatchControl | PatchSurgery],
    cfg: Config,
    noise: NoiseModel,
) -> stim.Circuit:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # -Loading in Patches
    ancilla_patch = patches["ancilla"]
    target_patch = patches["target"]
    control_patch = patches["control"]

    # -Retrieving Global Infomration
    q2i = lct.q2i
    distance = cfg.distance
    control_state_init = cfg.control_state_init
    target_state_init = cfg.target_state_init
    stab_to_data = lct.stab_to_data

    rounds = distance

    # -Retrieving Data Coords
    data_ancilla = ancilla_patch.data
    data_control = control_patch.data
    data_target = target_patch.data

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

    ########################
    # Define Initial Circuit
    ########################

    initial_circuit = stim.Circuit()

    """
    Looking at every state preperation seperatly seems to be inefficient
    ->  If not all Operators only once used one gets an incorrect formatting in the 
        timeslice view because of the Operations being in different timeslices in each TICK!
    """
    ########################################################################
    # Inilizing Ancilla in Plus (Reset) and Control/ Target in desired State
    ########################################################################

    init_patterns = {
        (a, b): [("RX", data_ancilla)]
        for a in ["Z0", "Z1", "X+", "X-", "Y+", "Y-"]
        for b in ["Z0", "Z1", "X+", "X-", "Y+", "Y-"]
    }

    # Apply the initialization pattern
    key = (control_state_init, target_state_init)

    if key not in init_patterns:
        raise ValueError(f"Invalid basis combination: {key}")

    for gate, qubits in init_patterns[key]:
        initial_circuit.append(gate, qubits)

    # -------Adding Before Round Data Depol.------------
    if noise.before_round_depol > 0:
        initial_circuit.append(
            "DEPOLARIZE1",
            data_ancilla + data_control + data_target,
            noise.before_round_depol,
        )
    # --------------------------------------------------

    initial_circuit.append("TICK")

    # Adding h gate for X stabilizers -> Filtering out double coords
    combined_x_stab: list = []
    for coords in x_stab_index_ancilla + x_stab_index_control + x_stab_index_target:
        if coords not in combined_x_stab:
            combined_x_stab.append(coords)

    initial_circuit.append("H", combined_x_stab)

    initial_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    cx_builder(q2i=q2i, stab_to_data=stab_to_data, circuit=initial_circuit, noise=noise)

    # Retreive Boundary + Normal Stabilizers Ancilla
    # (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!)
    initial_circuit.append("H", x_stab_index_ancilla)

    initial_circuit.append("TICK")

    initial_circuit.append("M", x_stab_index_ancilla + z_stab_index_ancilla)
    initial_circuit.append("TICK")
    initial_circuit.append("R", x_stab_index_ancilla + z_stab_index_ancilla)

    initial_circuit.append("TICK")
    initial_circuit.append("H", x_stab_boundary_b_index_ancilla)

    initial_circuit.append("TICK")

    # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=initial_circuit,
        orders=("5-CX", "6-CX"),
        noise=noise,
    )

    # Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    initial_circuit.append("H", x_stab_index_control + x_stab_index_target)

    initial_circuit.append("TICK")

    initial_circuit.append("M", control_target_stabs)

    ###########################################
    # Adding Repeat Block
    ###########################################

    initial_repeat_circuit = stim.Circuit()

    # Adding reset from initial round
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("R", control_target_stabs)

    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", combined_x_stab)

    initial_repeat_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    cx_builder(q2i=q2i, stab_to_data=stab_to_data, circuit=initial_repeat_circuit, noise=noise)

    # Retreive Boundary + Normal Stabilizers Ancilla
    # (Basis change and Measurement -> Measurement only in the x Basis UPDATE!!!!!)
    initial_repeat_circuit.append("H", x_stab_index_ancilla)

    initial_repeat_circuit.append("TICK")

    initial_repeat_circuit.append("M", x_stab_index_ancilla + z_stab_index_ancilla)
    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("R", x_stab_index_ancilla + z_stab_index_ancilla)

    initial_repeat_circuit.append("TICK")
    initial_repeat_circuit.append("H", x_stab_boundary_b_index_ancilla)

    initial_repeat_circuit.append("TICK")

    # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=initial_repeat_circuit,
        orders=("5-CX", "6-CX"),
        noise=noise,
    )

    # Retreive Boundary + Normal Stabilizers from Target and Control (Basis Change + Measurement):
    initial_repeat_circuit.append("H", x_stab_index_control + x_stab_index_target)

    initial_repeat_circuit.append("TICK")

    initial_repeat_circuit.append("M", control_target_stabs)

    initial_circuit += initial_repeat_circuit * (rounds - 1)

    return initial_circuit
