from dataclasses import dataclass

# --------------------------
# Shared Helper Dataclasses
# --------------------------

__all__ = ["NoiseParameters"]


@dataclass
class NoiseParameters:
    """
    Extended noise model to support XZZX-specific pauli-channel parameters.

    Notes:
    - before_round_p_xyz: per-qubit 1-body Pauli channel before a round (pX, pY, pZ)
    - after_c_p_xyz: per-qubit 1-body Pauli channel after Clifford (pX, pY, pZ)
    - after_c_p_xyz_multi: 2-qubit 15-probability Pauli channel after 2-qubit Clifford
    """

    # For Bias Noise Model
    before_round_p_xyz: list = None

    # For regular circuit Noise Model
    before_round_depol: float = 0.0
    before_m_flip_prob: float = 0.0
    after_r_flip: float = 0.0
    after_c_depol_prob: float = 0.0
    after_c_pauli_channel_prob: float = 0.0
    noise_bias: list[float] = None
