import stim

from qecsim.core.cls_circuit import GenCircuit

"""
TEST OVERLAPPING MEASUREMENTS
"""


def test_gen_circuit_addition_with_overlapping_measurements():
    """
    Cretaete two GenCircuits with overlapping measurements
    """

    test_c_1 = GenCircuit(circuit=stim.Circuit())

    test_c_1.add_operator("H", [0])
    test_c_1.add_measurement([0])

    test_c_2 = GenCircuit(circuit=stim.Circuit())
    test_c_2.add_operator("H", [0])
    test_c_2.add_measurement([0])

    # Build new circuit with overlapping detectors added
    combined_circuit = test_c_1 + test_c_2

    # Check that the combined circuit has the correct number of measurements and detectors
    assert combined_circuit.circuit.num_measurements == 2
    assert combined_circuit.circuit.num_detectors == 1
