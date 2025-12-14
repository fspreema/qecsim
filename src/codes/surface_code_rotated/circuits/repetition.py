import stim

from src.core.cx_builder import cx_builder
from src.core.data_models import (
    CircuitResult,
    ConfigSurface as Config,
    Context,
    NoiseModel,
    Patch,
)

Coord = complex

__all__ = ["repetition_circ"]


def repetition_circ(
    *,
    lct: Context,
    patches: dict[str, Patch],
    cfg: Config,
    noise: NoiseModel,
) -> CircuitResult:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # -Loading in Patches
    patch = patches["patch"]

    # -Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    rounds = cfg.rounds
    stab_to_data = lct.stab_to_data

    # -Retrieving Data Coords
    data = patch.data

    # -Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab

    ###########################
    # Define Repetition Circuit
    ###########################

    # -----BUILDING-REPETITION-CIRC------

    round_circuit = stim.Circuit()

    round_circuit.append("R", x_stab_index + z_stab_index)
    round_circuit.append("TICK")

    # -------Adding-Before-Round-Depol.-Data------------

    if noise.before_round_depol > 0:
        round_circuit.append("DEPOLARIZE1", data, noise.before_round_depol)

    # -------Continue-Circuit------------

    # 1) Reset/ Basis
    round_circuit.append("H", x_stab_index)

    round_circuit.append("TICK")

    # 2) CX Operations

    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=round_circuit,
    )

    # -------Continue-Circuit------------

    # 3) Basis/ Measurement
    round_circuit.append("H", x_stab_index)

    round_circuit.append("TICK")

    round_circuit.append("M", x_stab_index + z_stab_index)

    # -> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    round_circuit.append("SHIFT_COORDS", arg=(0, 0, 1))

    # 4) Detectors
    num_measurements_repeat = len(x_stab_index + z_stab_index)

    for index, q_index in enumerate(x_stab_index + z_stab_index):
        prev_tar = -2 * num_measurements_repeat + index
        current_tar = -1 * num_measurements_repeat + index
        # round_circuit.append(
        #     "DETECTOR",
        #     [stim.target_rec(current_tar), stim.target_rec(prev_tar)],
        #     (i2q[q_index].real, i2q[q_index].imag, 0),
        # )

    round_circuit.append("TICK")

    rep_circ = round_circuit * (rounds - 1)

    return CircuitResult(circuit=rep_circ)
