import stim

from src.core.cx_builder import cx_builder
from src.core.data_models import ConfigXZZX, PatchXZZX, XZZXContext

__all__ = ["initial"]


def initial(
    *,
    lct: XZZXContext,
    cfg: ConfigXZZX,
    patch: PatchXZZX,
) -> stim.Circuit:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    state_init = cfg.state_init

    # Retrieving Stabilizer to Data CX/CZ order
    stab_to_data = lct.stab_to_data

    # Retrieving Stabilizer Indices
    stab_index = patch.stab_index
    stab_index_ver = patch.stab_index_ver
    stab_index_hor = patch.stab_index_hor

    ##########################
    # Defining Initial Circuit
    ##########################

    initial_circuit = stim.Circuit()

    initial_circuit.append("TICK")
    initial_circuit.append("R", stab_index)

    initial_circuit.append("TICK")
    initial_circuit.append("H", stab_index)

    initial_circuit.append("TICK")

    ###############
    # CX Operations
    ###############

    cx_builder(
        circuit=initial_circuit,
        q2i=q2i,
        stab_to_data=stab_to_data,
        orders=("1-CX", "2-CZ", "3-CZ", "4-CX"),
    )

    # Retreive Boundary + Normal Stabilizers Ancilla:
    initial_circuit.append("H", stab_index)
    initial_circuit.append("TICK")

    initial_circuit.append("M", stab_index)
    initial_circuit.append("TICK")

    ##########################################################################
    # Implementing Detectors for Ancilla (inital Basis is deterministic)
    ##########################################################################

    # Determining Position in the measurement Run of only the Ancilla
    pos_to_index_ver: list = []
    pos_to_index_hor: list = []

    for pos, index in enumerate(stab_index):
        if index in stab_index_ver:
            pos_to_index_ver.append([pos, index])
        elif index in stab_index_hor:
            pos_to_index_hor.append([pos, index])

    if state_init in {"XZZX-VER"}:
        # Adding the needed Detectors
        for index_pos in pos_to_index_ver:
            current_tar = index_pos[0] - len(stab_index)
            q_index = index_pos[1]
            initial_circuit.append(
                "DETECTOR",
                [stim.target_rec(current_tar)],
                (i2q[q_index].real, i2q[q_index].imag, 0),
            )

    elif state_init in {"XZZX-HOR"}:
        # Adding the needed Detectors
        for index_pos in pos_to_index_hor:
            current_tar = index_pos[0] - len(stab_index)
            q_index = index_pos[1]
            initial_circuit.append(
                "DETECTOR",
                [stim.target_rec(current_tar)],
                (i2q[q_index].real, i2q[q_index].imag, 0),
            )

    else:
        raise ValueError("Not a valid Basis for initlization in the Control Lattice")

    return initial_circuit
