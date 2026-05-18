from typing import Any

import numpy as np
import pytest

from src.codes.lattice_surgery.builder import SurgeryBuilder
from src.tools.qem_estimator.logical_level.calc_ptm import PTMCalculator

PAULI_STRING_TO_INDEX_CASES = [
    # Single Pauli case
    ("I", 0),
    ("X", 1),
    ("Y", 2),
    ("Z", 3),
    # Two-Pauli case
    ("II", 0),
    ("IX", 1),
    ("IY", 2),
    ("IZ", 3),
    ("XI", 4),
    ("XX", 5),
    ("XY", 6),
    ("XZ", 7),
    ("YI", 8),
    ("YX", 9),
    ("YY", 10),
    ("YZ", 11),
    ("ZI", 12),
    ("ZX", 13),
    ("ZY", 14),
    ("ZZ", 15),
]

@pytest.mark.parametrize(("pauli_str", "expected_index"), PAULI_STRING_TO_INDEX_CASES)
def test_map_pauli_string_to_indices(pauli_str: str, expected_index: int) -> None:
    map_pauli_string_to_index = PTMCalculator._map_pauli_string_to_indices
    assert map_pauli_string_to_index(pauli_str) == expected_index

SGN_CASES_TWO_QUBIT = {
    "XI->XI": {"control": -1, "target": 1},
    "IX->IX": {"control": 1, "target": -1},
    "XI->XX": {"control": -1, "target": 1},
    "IX->ZX": {"control": 1, "target": -1},
}

@pytest.mark.parametrize(("curr_flow", "expected_sgn"), SGN_CASES_TWO_QUBIT.items())
def test_get_sgn_two_qubit(curr_flow: str, expected_sgn: dict[str, int]) -> None:
    sgn = PTMCalculator._get_sgn_two_qubit(curr_basis_combination= curr_flow)
    assert sgn == expected_sgn

SGN_CASES_ONE_QUBIT = {
    "I->X": 1,
    "X->X": -1,
    "Y->Z": -1,
    "Z->Z": -1,
    "I->Y": 1
}

@pytest.mark.parametrize(("curr_flow", "expected_sgn"), SGN_CASES_ONE_QUBIT.items())
def test_get_sgn_one_qubit(curr_flow: str, expected_sgn: dict[str, int]) -> None:
    sgn = PTMCalculator._get_sgn_one_qubit(curr_basis_combination= curr_flow)
    assert sgn == expected_sgn

LABEL_CASES = {
    "X+,Z0": ("X+", "Z0"),
    "Y-,Z1": ("Y-", "Z1"),
    "Z+,X0": ("Z+", "X0"),
}

@pytest.mark.parametrize(("joint_label", "split_label"), LABEL_CASES.items())
def test_split_label(joint_label: str, split_label: tuple[str, str]) -> None:
    label_splitted = PTMCalculator._split_label(init_state_label= joint_label)
    assert label_splitted == split_label

XOR_TEST_CASES = [
    # (decoder_prediction, logical_state_meas, expected)
    (
        np.array([[False], [False], [False]]), 
        np.array([0, 1, 0]), 
        np.array([1, -1, 1])
    ),
    (
        np.array([[False], [True], [False]]), 
        np.array([0, 1, 0]), 
        np.array([1, 1, 1])
    )
]

@pytest.mark.parametrize(("decoder_prediction", "logical_state_meas", "expected"), XOR_TEST_CASES)
def test_xor_meas_and_decoder_no_flip(decoder_prediction, logical_state_meas, expected) -> None:
    result = PTMCalculator._xor_meas_and_decoder(decoder_prediction, logical_state_meas)
    np.testing.assert_array_equal(result, expected)
