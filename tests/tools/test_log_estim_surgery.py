import numpy as np
import pytest

from src.tools.qem_estimator.logical_level.old.logical_estimator_surgery import LogicalEstimatorSurgery

COMMUTE_CASES = [
    ("XZ", "IZ", True),
    ("ZX", "XI", False),
]

@pytest.mark.parametrize(("basis_1", "basis_2", "expected_commute"), COMMUTE_CASES)
def test_does_commute(basis_1, basis_2, expected_commute):
    does_commute = LogicalEstimatorSurgery._does_commute(basis_1, basis_2)
    assert does_commute == expected_commute


def test_walsh_hadamard_transform():
    # Building Mtx
    mtx = LogicalEstimatorSurgery._walsh_hadamard_16()

    assert mtx.shape == (16, 16)
    assert np.all((mtx == 1) | (mtx == -1))

TEST_MTX = np.diag([0.3, 0.5, 0.2, 0.3, 0.8, 0.5, 0.2, 0.9, 0.7, 0.2, 0.1, 0.4, 0.7, 0.2, 0.3, 0.3])

def test_get_pauli_fidelities():
    fidelities = LogicalEstimatorSurgery._get_pauli_fidelities(TEST_MTX)

    assert fidelities.shape == (16,)
    assert np.all(fidelities >= 1)
