import pytest
import stim

from src.codes.xzzx.builder import XZZXBuilder
from src.core.data_models import NoiseParameters

# Define builder
def _build(
    distance: int,
    state_init: str,
    noise: NoiseParameters | None = None,
) -> stim.Circuit:
    builder = XZZXBuilder(
        distance=distance,
        state_init=state_init,
        noise=noise,
    )
    return builder.build_circuit()


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize("state_init",  ["XZZX-VER", "XZZX-HOR"], ids= ["XZZX-VER", "XZZX-HOR"])
def test_xzzx_builder(distance: int, state_init: str) -> None:
    return_circuit = _build(distance=distance, state_init=state_init)

    assert isinstance(return_circuit, stim.Circuit)
    assert return_circuit.num_qubits > 0
    assert return_circuit.num_measurements > 0

    # XZZX code: d² data qubits + (d²-1) ancilla = 2d²-1 total
    expected_qubits = 2 * distance**2 - 1
    assert return_circuit.num_qubits == expected_qubits

    # Circuit must define at least one logical observable
    assert "OBSERVABLE_INCLUDE" in str(return_circuit)

@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
def test_xzzx_builder_orientation_invariants(distance: int) -> None:
    """VER and HOR orientations must produce structurally identical circuits."""
    circuit_ver = _build(distance=distance, state_init="XZZX-VER")
    circuit_hor = _build(distance=distance, state_init="XZZX-HOR")

    # qubit and measurement counts must be identical!
    assert circuit_ver.num_qubits == circuit_hor.num_qubits
    assert circuit_ver.num_measurements == circuit_hor.num_measurements

@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize("state_init", ["XZZX-VER", "XZZX-HOR"], ids=["XZZX-VER", "XZZX-HOR"])
@pytest.mark.skip
def test_noisy_xzzx_builder(distance: int, state_init: str) -> None:
    noise = 1e-5
    noise_params = NoiseParameters(
        before_m_flip_prob=noise,
        after_r_flip=noise,
        after_c_depol_prob=noise,
        before_round_depol=noise,
    )

    circuit = _build(distance=distance, state_init=state_init, noise=noise_params)

    # Shortest graphlike error weight must equal the code distance
    assert len(circuit.shortest_graphlike_error()) == distance