import pytest
import stim

from src.codes.surface_code_rotated.builder import SurfaceBuilder
from src.core.data_models import NoiseParameters

# Define builder
def _build(
    distance: int,
    state_init: str,
    log_obs: str,
    noise: NoiseParameters | None = None,
) -> stim.Circuit:
    builder = SurfaceBuilder(
        distance=distance,
        state_init=state_init,
        log_obs=log_obs,
        noise=noise,
    )
    return builder.build_circuit()

@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize(
    "state_init,log_obs",
    [
        ("Z0", "Z"),
        ("Z1", "Z"),
        ("X+", "X"),
        ("X-", "X"),
        ("Y+", "Y"),
        ("Y-", "Y"),
    ],
)
def test_rotated_surface_code_qubits(distance: int, state_init: str, log_obs: str) -> None:

    return_circuit = _build(distance=distance, state_init=state_init, log_obs=log_obs)

    assert isinstance(return_circuit, stim.Circuit)
    assert return_circuit.num_qubits > 0
    assert return_circuit.num_measurements > 0

    # Rotated surface code: d² data qubits + (d²-1) ancilla = 2d²-1 total
    expected_qubits = 2 * distance**2 - 1
    assert return_circuit.num_qubits == expected_qubits

    # Circuit must define at least one logical observable
    assert "OBSERVABLE_INCLUDE" in str(return_circuit)


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize("basis", ["X", "Y", "Z"], ids= ["X", "Y", "Z"])
@pytest.mark.parametrize("log_obs", ["X", "Y", "Z"], ids= ["X", "Y", "Z"])
def test_surface_builder_invariants(distance, basis, log_obs) -> None:
    """
    Different eigenstates should return, other than the logical pauli string
    at the beginning, the exact same circuit! 
    """

    signs = {"Z": ("Z0", "Z1"), "X": ("X+", "X-"), "Y": ("Y+", "Y-")}
    state_a, state_b = signs[basis]

    circuit_a = _build(distance=distance, state_init=state_a, log_obs=log_obs)
    circuit_b = _build(distance=distance, state_init=state_b, log_obs=log_obs)

    # Different logical initialisations should only differ in Pauli frame —
    # qubit count and measurement count must be identical.
    assert circuit_a.num_qubits == circuit_b.num_qubits
    assert circuit_a.num_measurements == circuit_b.num_measurements

@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize(
    "state_init,log_obs",
    [
        ("Z0", "Z"),
        ("Z1", "Z"),
        ("X+", "X"),
        ("X-", "X"),
        ("Y+", "Y"),
        ("Y-", "Y"),
    ],
)
def test_noisy_surface_builder(distance: int, state_init: str, log_obs: str) -> None:
    noise = 1e-5
    noise_params = NoiseParameters(circuit_noise_prob= noise)

    circuit = _build(
        distance=distance,
        state_init=state_init,
        log_obs=log_obs,
        noise=noise_params,
    )

    # Shortest graphlike error weight must equal the code distance
    assert len(circuit.shortest_graphlike_error()) == distance