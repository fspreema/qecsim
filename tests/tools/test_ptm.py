from typing import Any

import numpy as np
import pytest

from src.codes.lattice_surgery.builder import SurgeryBuilder
from src.tools.qem_estimator.logical_level.calc_ptm import PTMCalculator

PAULI_STRING_TO_INDEX_CASES = [
    # Single Pauli case
    ("X", 0),
    ("Y", 1),
    ("Z", 2),
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

SGN_CASES = {
    "XI->XI": {"control": -1, "target": 1},
    "IX->IX": {"control": 1, "target": -1},
    "XI->XX": {"control": -1, "target": 1},
    "IX->ZX": {"control": 1, "target": -1},
}

@pytest.mark.parametrize(("curr_flow", "expected_sgn"), SGN_CASES.items())
def test_get_sgn(curr_flow: str, expected_sgn: dict[str, int]):
    sgn = PTMCalculator._get_sgn(curr_basis_combination= curr_flow)
    assert sgn == expected_sgn


LABEL_CASES = {
    "X+,Z0": ("X+", "Z0"),
    "Y-,Z1": ("Y-", "Z1"),
    "Z+,X0": ("Z+", "X0"),
}

@pytest.mark.parametrize(("joint_label", "split_label"), LABEL_CASES.items())
def test_split_label(joint_label: str, split_label: tuple[str, str]):
    label_splitted = PTMCalculator._split_label(init_state_label= joint_label)
    assert label_splitted == split_label

def test_get_entires_for_surgery():
    """
    Test if the get entries for surgery is working as expted
    -> Using flows which should have expectation value of 1
    -> Test Flow IX -> IX
    """

    # Create Circuits for testing
    SurgeryBuilder(
        distance=3,
        control_state_init="Z0",
        target_state_init="X+",
        control_measure_basis="Z",
        target_measure_basis="X",
    )

    # Adding all Circuit to the Dict

    pass
