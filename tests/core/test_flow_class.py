import stim

from src.core.unused.flow_builder import CircuitChunk

"""
This tests needs to be updated!!
-> We do not have the i2q available in this test case which we need for the funtion
for the locality extraction.
"""


def test_surgery_chunk_flow_extraction():
    """Test flow extraction from a weight-4 surface code stabilizer."""

    # Create a circuit with a single weight-4 X stabilizer measurement
    circuit = stim.Circuit()

    # Initialize data qubits and ancilla
    circuit.append("R", [0, 1, 2, 3, 4])  # 0,1,2,3 are data qubits, 4 is ancilla

    # Weight-4 X stabilizer: X0*X1*X2*X3
    circuit.append("H", [4])  # Put ancilla in |+> state
    circuit.append("CX", [4, 0])
    circuit.append("CX", [4, 1])
    circuit.append("CX", [4, 2])
    circuit.append("CX", [4, 3])
    circuit.append("H", [4])  # Measure in X basis
    circuit.append("M", [4])

    # Create a dummy i2q mapping
    i2q = {0: 0j, 1: 1 + 0j, 2: 0 + 1j, 3: 1 + 1j, 4: 0.5 + 0.5j}

    chunk = CircuitChunk(circuit=circuit, i2q=i2q)

    assert chunk.flows is not None
    assert isinstance(chunk.flows, (list, dict))
