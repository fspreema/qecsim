import stim

from qecsim.core.add_noise import BiasNoise, CircuitNoise
from qecsim.core.data_models import ConfigXZZX, PatchXZZX, XZZXContext, XZZXNoise
from qecsim.core.geometry import build_lattice
from qecsim.core.stabilizers import populate_stab_to_data

from .circuits.final_measure import final_measure
from .circuits.initial import initial
from .circuits.repetition import repetition
from .circuits.reset import reset

Coord = complex
Label = str

__all__ = ["xzzx_code"]

# -----------------
# Helper Fuctions
# -----------------


def _add_boundary_labels(distance: int, ancilla: dict[Coord, Label]) -> None:
    """
    Adds the neseccary Boundary and Surgery Stabilizers needed
    """

    max_coord = 2 * distance

    # Z-boundary stabilizers
    for y in range(2, max_coord, 4):
        coord_ancilla = complex(0, y)
        ancilla[coord_ancilla] = "STAB-BOUND-L-Hor"

    for y in range(4, max_coord, 4):
        coord_ancilla = complex(max_coord, y)
        ancilla[coord_ancilla] = "STAB-BOUND-R-Hor"

    # X-boundary stabilizers
    for y in range(4, max_coord, 4):
        coord_ancilla = complex(y, 0)
        ancilla[coord_ancilla] = "STAB-BOUND-A-Ver"

    for y in range(2, max_coord, 4):
        coord_ancilla = complex(y, max_coord)
        ancilla[coord_ancilla] = "STAB-BOUND-B-Ver"


# -----------------------------------------
# Public function -> Building final circuit
# -----------------------------------------


def xzzx_code(
    distance: int,
    rounds: int,
    *,
    state_init: str,
    before_round_depol: float = 0.0,
    before_round_p_xyz: list | None = None,
    before_m_flip_prob: float = 0.0,
    after_r_flip: float = 0.0,
    after_c_depol_prob: float = 0.0,
    noise_bias: list | None = None,
    after_c_pauli_channel_prob: float = 0.0,
) -> stim.Circuit:
    """
    Returns XZZX-Code circuit

    Arguments:
                -> state_init: In which basis should the lattice be initlized?

    Returns:
                -> Fully implemented XZZX-Code in stim.Circuit format
    """

    ###############################################################
    # 2. Build independent square patches (using geometry function)
    ###############################################################
    # Use the canonical core geometry builder. Pass state_init positionally
    # to match the existing xzzx API (state_init in {'Ver','Hor'}).
    qubit_coords: dict[Coord, Label] = build_lattice(distance, state_init)

    #######################################################################################
    # 3. Insert boundary & surgery labels (Only get activated in splitting/merging process)
    #######################################################################################
    _add_boundary_labels(distance, qubit_coords)

    ################################################################################
    # 4. Adding the Mapping from Stabilizer to Data for later CX gate implementation
    ################################################################################
    stab_to_data: dict[tuple[Coord, Coord], str] = populate_stab_to_data(qubit_coords)

    ###############################################
    # 5. Indexing All Qubits From given Coordinates
    ###############################################

    # Indexing Qubits
    q2i: dict[complex, int] = {
        q: i
        for i, q in enumerate(
            sorted(qubit_coords, key=lambda v: (v.real, v.imag)),
        )
    }

    # Reverse Indexing
    i2q: dict[int, complex] = {i: q for q, i in q2i.items()}

    # Filling Dataclasses
    noise = XZZXNoise(
        before_round_p_xyz=before_round_p_xyz,
        before_round_depol=before_round_depol,
        before_m_flip_prob=before_m_flip_prob,
        after_r_flip=after_r_flip,
        after_c_depol_prob=after_c_depol_prob,
    )

    lct = XZZXContext(
        q2i=q2i,
        i2q=i2q,
        stab_to_data=stab_to_data,
        coords=qubit_coords,
    )

    cfg = ConfigXZZX(
        distance=distance,
        state_init=state_init,
        rounds=rounds,
    )

    patches: dict[str, PatchXZZX] = {"patch": PatchXZZX.from_coords(qubit_coords, q2i)}

    ###############################
    # 6. Adding State Reset Circuit
    ###############################

    return_circuit = stim.Circuit()

    reset_circuit = reset(
        lct=lct,
        cfg=cfg,
        patch=patches["patch"],
        noise=noise,
    )

    return_circuit += reset_circuit

    ###########################
    # 7. Adding Initial Circuit
    ###########################

    initial_circuit = initial(
        lct=lct,
        cfg=cfg,
        patch=patches["patch"],
    )

    return_circuit += initial_circuit

    ##############################
    # 8. Adding Repetition Circuit
    ##############################

    repetition_circuit = repetition(
        lct=lct,
        cfg=cfg,
        patch=patches["patch"],
    )

    return_circuit += repetition_circuit

    #####################################
    # 9. Adding Final Measurement Circuit
    #####################################

    final_measure_circuit = final_measure(
        lct=lct,
        cfg=cfg,
        patch=patches["patch"],
    )

    return_circuit += final_measure_circuit

    #############################
    # 10. Returning Final Circuit
    #############################

    # Check what Noise needs to be added

    # 1) Circuit Noise Model
    if (
        before_m_flip_prob > 0.0
        or after_r_flip > 0.0
        or after_c_depol_prob > 0.0
        or before_round_depol > 0.0
    ):
        noise_dict = {
            "before_round_depol": before_round_depol,
            "before_m_flip_prob": before_m_flip_prob,
            "after_r_flip": after_r_flip,
            "after_c_depol_prob": after_c_depol_prob,
        }

        circuit_noise_builder = CircuitNoise(circuit=return_circuit, noise=noise_dict)

        return_circuit = circuit_noise_builder.apply()

    # 2) Biased Noise Model
    if after_c_pauli_channel_prob not in (0.0, None) or noise_bias not in (None, []):
        noise_dict = {
            "after_c_custom_noise": after_c_pauli_channel_prob,
            "bias": noise_bias,
        }

        bias_noise_builder = BiasNoise(circuit=return_circuit, noise=noise_dict)

        return_circuit = bias_noise_builder.apply()

    return return_circuit
