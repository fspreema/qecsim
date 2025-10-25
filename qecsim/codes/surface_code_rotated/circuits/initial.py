import stim

from qecsim.core.cx_builder import cx_builder
from qecsim.core.data_models import (
    CircuitResult,
    ConfigSurface as Config,
    Context,
    NoiseModel,
    Patch,
)

Coord = complex

__all__ = ["initial"]

def initial(*,
            lct: Context,
            patches: dict[str, Patch],
            cfg: Config,
            noise: NoiseModel) -> CircuitResult:

    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    patch = patches["patch"]

    #-Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    init_state = cfg.state_init
    stab_to_data = lct.stab_to_data

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab

    ########################
    # Define Initial Circuit
    ########################

    initial_circuit = stim.Circuit()

    initial_circuit.append("TICK")

    #-------Adding-Before-Round-Depol.-Data------------

    if noise.before_round_depol > 0:
        initial_circuit.append("DEPOLARIZE1", data, noise.before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    initial_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #2) CX Operations

    cx_builder(q2i= q2i,
               stab_to_data= stab_to_data,
               circuit= initial_circuit,
               noise= noise)

    #-------Continue-Circuit------------

    #3) Basis/ Measurement
    initial_circuit.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        initial_circuit.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.------------

    if noise.before_m_flip_prob > 0:
        initial_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.before_m_flip_prob)

    #-------Continue-Circuit------------

    initial_circuit.append("MR", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if noise.after_r_flip > 0:
        initial_circuit.append("X_ERROR", x_stab_index + z_stab_index, noise.after_r_flip)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #4) DETECTORS -> Measure only deterministic-Stabilizers!

    if init_state in {"0", "1"}:
        num_measurements_initial = len(z_stab_index)

        for index, q_index in enumerate(z_stab_index):
            current_tar = -1 * num_measurements_initial + index
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], 
                                   (i2q[q_index].real, i2q[q_index].imag, 0))

    elif init_state in {"+", "-"}:
        num_measurements_initial = len(x_stab_index + z_stab_index)

        for index, q_index in enumerate(x_stab_index):
            current_tar = -1 * num_measurements_initial + index
            initial_circuit.append("DETECTOR", [stim.target_rec(current_tar)], 
                                   (i2q[q_index].real, i2q[q_index].imag, 0))

    return CircuitResult(circuit=initial_circuit)
