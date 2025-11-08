import stim

from qecsim.core.data_models import ConfigXZZX, PatchXZZX, XZZXContext

__all__ = ["repetition"]


def repetition(
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
    rounds = cfg.rounds

    # Retrieving Stabilizer to Data CX/CZ order
    stab_to_data = lct.stab_to_data

    # Retrieving Stabilizer Indices
    stab_index = patch.stab_index
    stab_index_ver = patch.stab_index_ver
    stab_index_hor = patch.stab_index_hor

    ###########################################
    # Adding Repeat Block
    ###########################################

    repeat_circuit = stim.Circuit()

    repeat_circuit.append("SHIFT_COORDS", arg=(0, 0, 1))
    repeat_circuit.append("R", stab_index)

    repeat_circuit.append("TICK")

    repeat_circuit.append("H", stab_index)

    repeat_circuit.append("TICK")

    ####################################################
    # CX Operations
    ####################################################

    for coord_pairs, order in stab_to_data.items():
        # Parallel Implementation of CX
        if order == "1-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            repeat_circuit.append("CX", index_pairs)

    repeat_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
        # Parallel Implementation of CX
        if order == "2-CZ":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            repeat_circuit.append("CZ", index_pairs)

    repeat_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
        # Parallel Implementation of CX
        if order == "3-CZ":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            repeat_circuit.append("CZ", index_pairs)

    repeat_circuit.append("TICK")

    for coord_pairs, order in stab_to_data.items():
        # Parallel Implementation of CX
        if order == "4-CX":
            index_pairs = []
            index_pairs.append(q2i[coord_pairs[1]])
            index_pairs.append(q2i[coord_pairs[0]])
            repeat_circuit.append("CX", index_pairs)

    # Retreive Boundary + Normal Stabilizers Ancilla:
    repeat_circuit.append("TICK")
    repeat_circuit.append("H", stab_index)

    repeat_circuit.append("TICK")

    repeat_circuit.append("M", stab_index)

    repeat_circuit.append("TICK")

    ######################################################
    # Implementing Detectors (All Basis are deterministic)
    ######################################################

    # Determining Position in the measurement Run of only the Ancilla
    pos_to_index_ver: list = []
    pos_to_index_hor: list = []

    for pos, index in enumerate(stab_index):
        if index in stab_index_ver:
            pos_to_index_ver.append([pos, index])
        elif index in stab_index_hor:
            pos_to_index_hor.append([pos, index])

    """
    As intial circuit run is completed, now we define all stabilizers in every basis
    -> Detectors on both basis are now deterministic
    """

    # Adding the needed Detectors (Z-Basis)
    for index_pos in pos_to_index_ver:
        current_tar = index_pos[0] - len(stab_index)
        previous_tar = index_pos[0] - 2 * len(stab_index)
        q_index = index_pos[1]
        repeat_circuit.append(
            "DETECTOR",
            [stim.target_rec(current_tar), stim.target_rec(previous_tar)],
            (i2q[q_index].real, i2q[q_index].imag, 0),
        )

    # Adding the needed Detectors (X-Basis)
    for index_pos in pos_to_index_hor:
        current_tar = index_pos[0] - len(stab_index)
        previous_tar = index_pos[0] - 2 * len(stab_index)
        q_index = index_pos[1]
        repeat_circuit.append(
            "DETECTOR",
            [stim.target_rec(current_tar), stim.target_rec(previous_tar)],
            (i2q[q_index].real, i2q[q_index].imag, 0),
        )

    return repeat_circuit * (rounds - 1)
