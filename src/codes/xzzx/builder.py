import stim

from src.codes.xzzx.get_stab_pairings import XZZXPairings
from src.codes.xzzx.xzzx_geom import XZZXGeometry
from src.core.data_models import ConfigXZZX, PatchXZZX, XZZXContext, XZZXNoise
from src.core.noise_models import BiasNoise, CircuitNoise

from .circuits.final_measure import final_measure
from .circuits.initial import initial
from .circuits.repetition import repetition
from .circuits.reset import reset

Coord = complex
Label = str

__all__ = ["xzzx_code"]

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

    ########################
    # 2. Build XZZX Geometry
    ########################
    class_geometry = XZZXGeometry(
        distance=distance,
        state_init=state_init,
    )

    ################################################################################
    # 4. Adding the Mapping from Stabilizer to Data for later CX gate implementation
    ################################################################################

    pairings = XZZXPairings(patch=class_geometry.coords)

    # Filling Dataclasses
    noise = XZZXNoise(
        before_round_p_xyz=before_round_p_xyz,
        before_round_depol=before_round_depol,
        before_m_flip_prob=before_m_flip_prob,
        after_r_flip=after_r_flip,
        after_c_depol_prob=after_c_depol_prob,
    )

    lct = XZZXContext(
        q2i=class_geometry._get_q2i(),
        i2q=class_geometry._get_i2q(),
        stab_to_data=pairings.get_schedule(),
        coords=class_geometry.coords,
    )

    cfg = ConfigXZZX(
        distance=distance,
        state_init=state_init,
        rounds=rounds,
    )

    patches: dict[str, PatchXZZX] = {
        "patch": PatchXZZX.from_coords(
            class_geometry.coords,
            class_geometry._get_q2i(),
        ),
    }

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
