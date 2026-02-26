from dataclasses import dataclass, field
from enum import Enum, auto

import stim

__all__ = ["NoiseParameters",
           "PTMCircuits"]


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
    noise_bias: list[float] | None = None


### FOR FUTURE IMPLEMENTATION ###


class MeasurementBasis(Enum):
    X = auto()
    Y = auto()
    Z = auto()


class InitilizationState(Enum):
    PLUS = auto()
    MINUS = auto()
    PLUS_I = auto()
    MINUS_I = auto()
    ZERO = auto()
    ONE = auto()


class SurgeryStabilizerType(Enum):
    X_STABILIZER = auto()
    Z_STABILIZER = auto()


class SurgeryOperationType(Enum):
    ANCILLA_CONTROL = auto()
    ANCILLA_TARGET = auto()


class SurgeryPatchTypes(Enum):
    CONTROL = auto()
    TARGET = auto()
    ANCILLA = auto()


##################################


# Dataclasses for construction of the PTM
@dataclass
class PTMCircuits:
    """
    Data model to store all circuits required for PTM calculation.

    Functionality:
    - Stores both deterministic and non-deterministic circuits for PTM calculations
    - Str: Organizes Circuits by input-output Pauli Combinations
    - For each
    - Returns a tuple with Circuit and Measurement Records
        of logical operator to alter xor with observable flipped value given by dem

    ### Lattice Surgery ###

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

    ### Surface Code ###
    This is just a 4 x 4 Matrix

     ,I,X,Y,Z
    I,1,0,0,0
    X,0,1,0,0
    Y,0,0,1,0
    Z,0,0,0,1

    """

    # Circuits
    # Structure: {basis_combination: ({init_state_label:
    # stim_circuit}, [meas_rec] )}
    circuits: dict[str, tuple[dict[str, stim.Circuit], list[int]]] = field(
        default_factory=dict,
    )


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
