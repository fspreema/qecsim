import stim

from qecsim.core.cx_builder import cx_builder
from qecsim.core.data_models import (
    CircuitResult,
    Context,
    NoiseModel,
    Patch,
)

Coord = complex

__all__ = ["h_switched_circ_init"]

def h_switched_circ_init(*,
                        lct: Context,
                        patches: dict[str, Patch],
                        noise: NoiseModel) -> CircuitResult:

    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    patch = patches["patch"]

    #-Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    stab_to_data = lct.stab_to_data_modified

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.z_stab
    z_stab_index = patch.x_stab

    ###########################
    # Define Repetition Circuit
    ###########################

    #-----BUILDING-REPETITION-CIRC------

    switched_init_circ = stim.Circuit()

    # Adding H gate on data

    switched_init_circ.append("H", data)
    switched_init_circ.append("TICK")

    #-------Adding-Before-Round-Depol.-Data------------

    if noise.before_round_depol > 0:
        switched_init_circ.append("DEPOLARIZE1", data, noise.before_round_depol)

    #-------Continue-Circuit------------

    #1) Reset/ Basis
    switched_init_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        switched_init_circ.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    switched_init_circ.append("TICK")

    #2) CX Operations

    cx_builder(q2i= q2i,
               stab_to_data= stab_to_data,
               circuit= switched_init_circ,
               noise= noise)

    #3) Basis/ Measurement
    switched_init_circ.append("H", x_stab_index)

    #-------Adding-After-Clifford-Depol.------------

    if noise.after_c_depol_prob > 0:
        switched_init_circ.append("DEPOLARIZE1", x_stab_index, noise.after_c_depol_prob)

    #-------Continue-Circuit------------

    switched_init_circ.append("TICK")

    #-------Adding-Before-Measurement-Flip-Prob.-------

    if noise.before_m_flip_prob > 0:
        switched_init_circ.append("X_ERROR", x_stab_index + z_stab_index, noise.before_m_flip_prob)

    #-------Continue-Circuit----------

    """
    Due to the fact that x and z indicies are flipped, measurement order is changed!
    """

    switched_init_circ.append("MR", x_stab_index + z_stab_index)

    #-------Adding-After-Reset-Flip-Prob.------------

    if noise.after_r_flip > 0:
        switched_init_circ.append("X_ERROR", x_stab_index + z_stab_index, noise.after_r_flip)

    #-------Continue-Circuit------------

    #-> Shifting Coords in Time-Dimension to have 3D timelike Detector graph (Needed for decoding)
    switched_init_circ.append("SHIFT_COORDS", arg= (0, 0, 1))

    """
    For the detectors we now compare the newly formed stabilizers with the old ones on the same position!
    -> Z stabilizers (Here after switch) lays on the old x stabilizers
    -> We compare this Z stabilizer to the old X stabilizer
    -> We compare current qubit index with the same qubit index of old emasurement round
    """

    #4) Detectors
    num_measurements_repeat = len(x_stab_index + z_stab_index)

    # X -> Z Stabilizers & Z -> X Stabilizers
    for index, q_index in enumerate(x_stab_index + z_stab_index):

        # Currently Z Stab
        if q_index in z_stab_index:
            prev_tar = int(-2.5 * num_measurements_repeat + index)
            current_tar = -1 * num_measurements_repeat + index
            switched_init_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)],
                                (i2q[q_index].real, i2q[q_index].imag, 0))

        # Currently Z Stab
        if q_index in x_stab_index:
            prev_tar = int(-1.5 * num_measurements_repeat + index)
            current_tar = -1 * num_measurements_repeat + index
            switched_init_circ.append("DETECTOR", [stim.target_rec(current_tar),stim.target_rec(prev_tar)],
                                (i2q[q_index].real, i2q[q_index].imag, 0))

    switched_init_circ.append("TICK")

    return CircuitResult(circuit=switched_init_circ)
