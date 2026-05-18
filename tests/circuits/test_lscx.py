import pytest
import stim

from src.codes.lattice_surgery.builder import SurgeryBuilder
from src.core.get_measurement_recs import get_measurement_recs
from src.core.data_models import NoiseParameters

# Example Circuits
CIRCUITS = {1:
(stim.Circuit("""
R 0
R 1
M 0
M 1
OBSERVABLE_INCLUDE(0) rec[-1]
"""),[-1]),
2: (stim.Circuit("""
R 0
R 1
M 0
M 1
OBSERVABLE_INCLUDE(0) rec[-1]
M 0
M 1
OBSERVABLE_INCLUDE(0) rec[-1]
"""), [-3, -1] ),
3: (stim.Circuit("""
R 0
R 1
M 0
M 1
OBSERVABLE_INCLUDE(0) rec[-1]
R 0
R 1
M 0
M 1
OBSERVABLE_INCLUDE(0) rec[-1]
M 1
"""), [-4, -2]),
}

@pytest.mark.parametrize("circuit_id,expected_recs", CIRCUITS.items())
def test_get_measurement_recs(circuit_id: int, expected_recs: tuple[stim.Circuit, list[int]]) -> None:
    circuit = CIRCUITS[circuit_id][0]
    recs = get_measurement_recs(circuit=circuit, observable_index=0)
    assert recs == expected_recs[1]
    assert isinstance(recs, list)
    assert all(r < 0 for r in recs)
    assert all(-circuit.num_measurements <= r <= -1 for r in recs)

# Construct all possible surgery types (256 circuits)
# Pauli alphabet for input and output states
PAULIS = ["I", "X", "Y", "Z"]

# Create Mapping for input states
input_to_init_state = {
    "X": ["X+", "X-"],
    "Y": ["Y+", "Y-"],
    "Z": ["Z0", "Z1"],
    "I": ["I0", "I1"],
}

# Define Builder
def _build(distance: int,
           p_in_c: str,
           p_in_t: str,
           p_out_c: str,
           p_out_t: str,
           sign_idx: int,
           noise = None) -> stim.Circuit:

    builder = SurgeryBuilder(
        distance=distance,
        control_state_init=input_to_init_state[p_in_c][sign_idx],
        target_state_init=input_to_init_state[p_in_t][sign_idx],
        control_measure_basis=p_out_c,
        target_measure_basis=p_out_t,
        noise= noise,
    )
    return builder.build_circuit()

@pytest.mark.parametrize("distance", [3,5], ids=["d3", "d5"])
@pytest.mark.parametrize(
    "p_in_c,p_in_t,p_out_c,p_out_t",
    [
        ("I", "X", "I", "X"),
        ("X", "I", "X", "X"),
        ("Z", "X", "Z", "X"),
        ("X", "X", "X", "I"),
        ("Z", "Y", "I", "Y"),
    ],
)
def test_surgery_builder(distance: int,
                         p_in_c: str,
                         p_in_t: str,
                         p_out_c: str,
                         p_out_t: str,
                         ) -> None:
    circuit_1 = _build(distance=distance,
                     p_in_c=p_in_c,
                     p_in_t=p_in_t,
                     p_out_c=p_out_c,
                     p_out_t=p_out_t,
                     sign_idx=0)
    circuit_2 = _build(distance=distance,
                       p_in_c=p_in_c,
                       p_in_t=p_in_t,
                       p_out_c=p_out_c,
                       p_out_t=p_out_t,
                       sign_idx=1)
    assert isinstance(circuit_1, stim.Circuit)
    assert isinstance(circuit_2, stim.Circuit)
    assert circuit_1.num_qubits > 0
    assert circuit_2.num_qubits > 0
    assert circuit_1.num_measurements > 0
    assert circuit_2.num_measurements > 0

    # Different logical initilizations should only add strings of pauli
    # operators to the circuit
    assert circuit_1.num_qubits == circuit_2.num_qubits
    assert circuit_1.num_measurements == circuit_2.num_measurements

@pytest.mark.parametrize("distance", [3,5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize(
    "p_in_c,p_in_t,p_out_c,p_out_t",
    [
        ("I", "X", "I", "X"),
        ("X", "I", "X", "X"),
        ("Z", "X", "Z", "X"),
        ("X", "X", "X", "I"),
        ("Z", "Y", "I", "Y"),
    ],
)
def test_noisy_surgery_builder(distance: int,
                         p_in_c: str,
                         p_in_t: str,
                         p_out_c: str,
                         p_out_t: str,
                         ) -> None:
    
    noise = 1e-5

    noise_class_circuit = NoiseParameters(before_m_flip_prob=noise,
                                after_r_flip=noise,
                                after_c_depol_prob=noise,
                                before_round_depol=noise)
    
    circuit_1 = _build(distance=distance,
                     p_in_c=p_in_c,
                     p_in_t=p_in_t,
                     p_out_c=p_out_c,
                     p_out_t=p_out_t,
                     sign_idx=0,
                     noise= noise_class_circuit)
    
    # Shortest logical error should be equal to the distance
    assert len(circuit_1.shortest_graphlike_error()) == distance