import numpy as np
import pytest

from src.tools.qem_estimator.logical_level.logical_estimator import GeneralLogicalEstimator

COMMUTE_CASES = [
    ("X", "X", True),
    ("X", "Y", False),
    ("I", "Z", True),
    ("XZ", "YY", True),
    ("II", "II", True),
    ("XZ", "IZ", True),
    ("ZX", "XI", False),
    # Although not used, testing is still beneficial!
    ("XZI", "YYI", True),
    ("XZI", "YII", False),
]

@pytest.mark.parametrize(("basis_1", "basis_2", "expected_commute"), COMMUTE_CASES)
def test_does_commute(basis_1, basis_2, expected_commute) -> None:
    does_commute = GeneralLogicalEstimator._does_commute(basis_1, basis_2)
    assert does_commute == expected_commute

VECTOR_STACKING_CASES= [
    (np.array([[1, 0], [0, 1]]), np.array([1, 0, 0, 1])),
    (np.array([[1, 0, 1], [0, 1, 1], [1, 0, 0]]), np.array([1, 0, 1, 0, 1, 0, 1, 1, 0]))
]

@pytest.mark.parametrize(("input_ptm", "expected_vec"), VECTOR_STACKING_CASES)
def test_stack_vectors_of_mtrx(input_ptm, expected_vec) -> None:
    estimator = GeneralLogicalEstimator.__new__(GeneralLogicalEstimator)
    estimator.n_rows, _ = input_ptm.shape
    result = estimator._stack_vectors_of_mtrx(input_ptm)
    # Check if stacking is as expected and dim is correct
    np.testing.assert_array_equal(result, expected_vec)

TEST_NUM_QBT_CASES= [
    (np.zeros((4,4)), 1),
    (np.zeros((16,16)), 2)
]

@pytest.mark.parametrize(("ptm_dummy", "expected_qbt_num"), TEST_NUM_QBT_CASES)
def test_num_qubits(ptm_dummy, expected_qbt_num) -> None:
    estimator = GeneralLogicalEstimator(ptm_dummy, ptm_dummy)
    assert estimator.num_qubits == expected_qbt_num

def test_get_inv_noise_mtx_no_noise() -> None:
    # When noisy == ideal, noise matrix should be identity
    # -> Therefore I^-1 = I is expected!
    ptm = np.diag([1.0, 0.9, 0.8, 0.95])
    estimator = GeneralLogicalEstimator(ptm_noisy=ptm, ptm_ideal=ptm)
    result = estimator._get_inv_noise_mtx()
    np.testing.assert_allclose(result, np.eye(4), atol=1e-10)