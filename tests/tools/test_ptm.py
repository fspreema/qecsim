import numpy as np

from src.core.data_models import PTMCircuits
from src.tools.qem_estimator.logical_level.calc_ptm import PTMCalculator


def test_pauli_to_string():
    # Test for Single Pauli Case
    assert PTMCalculator._map_pauli_string_to_indices("X") == 0
    assert PTMCalculator._map_pauli_string_to_indices("Y") == 1
    assert PTMCalculator._map_pauli_string_to_indices("Z") == 2

    # Test for Two Pauli Case
    assert PTMCalculator._map_pauli_string_to_indices("II") == 0
    assert PTMCalculator._map_pauli_string_to_indices("IX") == 1
    assert PTMCalculator._map_pauli_string_to_indices("IY") == 2
    assert PTMCalculator._map_pauli_string_to_indices("IZ") == 3
    assert PTMCalculator._map_pauli_string_to_indices("XI") == 4
    assert PTMCalculator._map_pauli_string_to_indices("XX") == 5
    assert PTMCalculator._map_pauli_string_to_indices("XY") == 6
    assert PTMCalculator._map_pauli_string_to_indices("XZ") == 7
    assert PTMCalculator._map_pauli_string_to_indices("YI") == 8
    assert PTMCalculator._map_pauli_string_to_indices("YX") == 9
    assert PTMCalculator._map_pauli_string_to_indices("YY") == 10
    assert PTMCalculator._map_pauli_string_to_indices("YZ") == 11
    assert PTMCalculator._map_pauli_string_to_indices("ZI") == 12
    assert PTMCalculator._map_pauli_string_to_indices("ZX") == 13
    assert PTMCalculator._map_pauli_string_to_indices("ZY") == 14
    assert PTMCalculator._map_pauli_string_to_indices("ZZ") == 15


def test_construct_ptm_matrix_surface():
    # Init Dummy Expectation Values for both Lattice Surgery and Surface Code
    exp_vals_surface = {
        "X->X": 1.0,
        "Y->Y": 1.0,
        "Z->Z": 1.0,
    }

    # Get PTM Matrices
    ptm_calculator = PTMCalculator(PTMCircuits(circuits={}), samples=1_000)
    ptm_matrix_surface = ptm_calculator._build_mtx_for_ptm(exp_vals_surface)

    # Assert that the PTM Matrices are correct
    expected_ptm_surface = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    assert np.array_equal(ptm_matrix_surface, expected_ptm_surface)
