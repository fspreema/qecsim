from src.codes.lattice_surgery.builder import SurgeryBuilder
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


def test_get_sgn():
    # Test for Identity Output
    assert PTMCalculator._get_sgn("XI->XI") == {"control": -1, "target": 1}
    assert PTMCalculator._get_sgn("IX->IX") == {"control": 1, "target": -1}

    # Test for Non-Identity Output
    assert PTMCalculator._get_sgn("XI->XX") == {"control": -1, "target": 1}
    assert PTMCalculator._get_sgn("IX->ZX") == {"control": 1, "target": -1}


def test_split_label():
    # Test Number of State Init Combinations
    assert PTMCalculator._split_label("X+,Z0") == ("X+", "Z0")
    assert PTMCalculator._split_label("Y-,Z1") == ("Y-", "Z1")
    assert PTMCalculator._split_label("Z+,X0") == ("Z+", "X0")


def test_get_init_pairing_value():
    # Create dummy logical estimate dictionary
    logical_estimate_dict = {"X+,Z0": 0.8, "X+,Z1": 0.2, "X-,Z0": 0.3, "X-,Z1": 0.7}

    # Get value for pairings
    assert (
        PTMCalculator._get_init_pairing_value(
            pair=(0, 0),
            average_logical_state=logical_estimate_dict,
        )
        == 0.8
    )
    assert (
        PTMCalculator._get_init_pairing_value(
            pair=(0, 1),
            average_logical_state=logical_estimate_dict,
        )
        == 0.2
    )
    assert (
        PTMCalculator._get_init_pairing_value(
            pair=(1, 0),
            average_logical_state=logical_estimate_dict,
        )
        == 0.3
    )
    assert (
        PTMCalculator._get_init_pairing_value(
            pair=(1, 1),
            average_logical_state=logical_estimate_dict,
        )
        == 0.7
    )


def test_get_entires_for_surgery():
    """
    Test if the get entries for surgery is working as expted
    -> Using flows which should have expectation value of 1
    -> Test Flow IX -> IX
    """

    # Create Circuits for testing
    circ_1 = SurgeryBuilder(
        distance=3,
        control_state_init="Z0",
        target_state_init="X+",
        control_measure_basis="Z",
        target_measure_basis="X",
    )

    # Adding all Circuit to the Dict
    circuits_surgery = PTMCircuits(circuits=[])

    pass
