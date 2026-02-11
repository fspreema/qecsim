from dataclasses import dataclass

import stim

__all__ = ["NoiseParameters"]

# Dataclass for saving Noise Parameters which will be given to the individual Circuit Builders


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


# Dataclasses for construction of the PTM


@dataclass
class PTMCircuitsSurface:
    """
    Data model to store all circuits required for PTM calculation.

    Attributes:
        ideal_circuit: Ideal circuit without noise
        -> These are only needed for the diagonal entries
        noisy_circuit: Noisy circuit with noise model applied
    """

    # Noisy Circuits
    circ_xx_noisy: stim.Circuit
    circ_xy_noisy: stim.Circuit
    circ_xz_noisy: stim.Circuit
    circ_yx_noisy: stim.Circuit
    circ_yy_noisy: stim.Circuit
    circ_yz_noisy: stim.Circuit
    circ_zx_noisy: stim.Circuit
    circ_zy_noisy: stim.Circuit
    circ_zz_noisy: stim.Circuit

    # Ideal Circuits
    circ_xx_ideal: stim.Circuit
    circ_yy_ideal: stim.Circuit
    circ_zz_ideal: stim.Circuit


@dataclass
class PTMCircuitsSurgery:
    """
    Data model to store all circuits required for PTM calculation in lattice surgery.

    Attributes:
        ideal_circuit: Ideal circuit without noise
        noisy_circuit: Noisy circuit with noise model applied

    These are essentially buidling up the full PTM by seperatly calculating the different
    logical input to logical output combinations.

      ,II,IX,IY,IZ,XI,XX,XY,XZ,YI,YX,YY,YZ,ZI,ZX,ZY,ZZ
    II,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
    IX,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
    IY,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1
    IZ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1
    XI,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
    XX,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
    XY,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0
    XZ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,-1 ,0 ,0 ,0 ,0 ,0
    YI,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0
    YX,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0
    YY,0 ,0 ,0 ,0 ,0 ,0 ,0 ,-1,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
    YZ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
    ZI,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0
    ZX,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0
    ZY,0 ,0 ,-1,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0
    ZZ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0


    It should be noted that all II can't be produced as a measurement needs to be made
    or rather some obseravble needs to be declared as well as we need to reset the qubits

    """

    # Noisy Circuits
    # Keys should follow the pattern: f"{input_pauli}_{output_pauli}" (e.g., "XI_XX")
    circuits_noisy: dict[str, stim.Circuit]

    # Ideal Circuits
    circuits_ideal: dict[str, stim.Circuit]


@dataclass
class PTMExpectedResultsSurface:
    """
    Data model to store expected results of the deterministic diagonal
    entries of the circuit for PTM calculations.
    """


@dataclass
class PTMExpectedResultsSurgery:
    """
    Data model to store expected results of the deterministic diagonal
    entries of the lattice surgery circuit for PTM calculations.
    """
