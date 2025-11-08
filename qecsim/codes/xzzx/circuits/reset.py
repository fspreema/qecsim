import numpy as np
import stim

from qecsim.core.data_models import ConfigXZZX, PatchXZZX, XZZXContext, XZZXNoise

__all__ = ["reset"]


def reset(
    *,
    lct: XZZXContext,
    cfg: ConfigXZZX,
    patch: PatchXZZX,
    noise: XZZXNoise,
) -> stim.Circuit:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # Retrieving Global Infomration
    q2i = lct.q2i
    state_init = cfg.state_init

    # Retrievin Data Coords from Patch
    data_x = patch.data_x
    data_z = patch.data_z
    data = patch.data

    ########################
    # Define Initial Circuit
    ########################

    initial_circuit = stim.Circuit()

    """
    Looking at every state preperation seperatly seems to be inefficient
    ->  If not all Operators only once used one gets an incorrect formatting in the 
        timeslice view because of the Operations being in different timeslices in each TICK!
    """

    ##################
    # Appending Coords
    ##################

    for q, i in q2i.items():
        initial_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

    ########################################################################
    # Inilizing Ancilla in Plus (Reset) and Control/ Target in desired State
    ########################################################################

    init_patterns = {
        ("Ver"): [("RZ", data_z), ("RX", data_x)],
        ("Hor"): [("RX", data_x), ("RZ", data_z)],
    }

    # Apply the initialization pattern
    key = state_init

    if key not in init_patterns:
        raise ValueError(f"Invalid basis combination: {key}")

    for gate, qubits in init_patterns[key]:
        initial_circuit.append(gate, qubits)

    # -------Adding Before Round Data Depol.------------
    if noise.before_round_depol > 0:
        initial_circuit.append("DEPOLARIZE1", data, noise.before_round_depol)
    # --------------------------------------------------

    # -------Adding Before Round Data Depol.------------
    if np.any(noise.before_round_p_xyz):
        initial_circuit.append("PAULI_CHANNEL_1", data, noise.before_round_p_xyz)
    # --------------------------------------------------

    return initial_circuit
