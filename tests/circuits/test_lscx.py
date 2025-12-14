import pytest
import stim

from src.codes.lattice_surgery.builder import surgery_circuit


def _has_observable(circuit: stim.Circuit) -> bool:
    return "OBSERVABLE_INCLUDE" in str(circuit)


@pytest.mark.parametrize("distance", [3, 5], ids=["d3", "d5"])
@pytest.mark.parametrize(
    "flow,control,target",
    [
        # All-X: X -> XX
        ("X -> XX", "X+", "X+"),
        ("X -> XX", "X+", "X-"),
        ("X -> XX", "X-", "X+"),
        ("X -> XX", "X-", "X-"),
        # All-X: XX -> X
        ("XX -> X", "X+", "X+"),
        ("XX -> X", "X+", "X-"),
        ("XX -> X", "X-", "X+"),
        ("XX -> X", "X-", "X-"),
        # All-X: X -> X
        ("X -> X", "X+", "X+"),
        ("X -> X", "X+", "X-"),
        ("X -> X", "X-", "X+"),
        ("X -> X", "X-", "X-"),
        ("X -> X", "Z0", "X+"),
        ("X -> X", "Z0", "X-"),
        ("X -> X", "Z1", "X+"),
        ("X -> X", "Z1", "X-"),
        # All-Z: Z -> ZZ
        ("Z -> ZZ", "Z0", "Z0"),
        ("Z -> ZZ", "Z0", "Z1"),
        ("Z -> ZZ", "Z1", "Z0"),
        ("Z -> ZZ", "Z1", "Z1"),
        # All-Z: ZZ -> Z
        ("ZZ -> Z", "Z0", "Z0"),
        ("ZZ -> Z", "Z0", "Z1"),
        ("ZZ -> Z", "Z1", "Z0"),
        ("ZZ -> Z", "Z1", "Z1"),
        # All-Z: Z -> Z
        ("Z -> Z", "Z0", "Z0"),
        ("Z -> Z", "Z0", "Z1"),
        ("Z -> Z", "Z1", "Z0"),
        ("Z -> Z", "Z1", "Z1"),
        ("Z -> Z", "Z0", "X+"),
        ("Z -> Z", "Z0", "X-"),
        ("Z -> Z", "Z1", "X+"),
        ("Z -> Z", "Z1", "X-"),
        # Mixed ZX: ZX -> ZX
        ("ZX -> ZX", "Z0", "X+"),
        ("ZX -> ZX", "Z0", "X-"),
        ("ZX -> ZX", "Z1", "X+"),
        ("ZX -> ZX", "Z1", "X-"),
    ],
)
def test_lscx_valid_flows_compile_and_have_observable(distance, flow, control, target):
    """Each supported flow compiles for valid state combos and includes observable."""
    circuit = surgery_circuit(
        distance=distance,
        target_state_init=target,
        control_state_init=control,
        flow_observable=flow,
    )

    assert isinstance(circuit, stim.Circuit)
    assert _has_observable(circuit)


@pytest.mark.parametrize("d", [3, 5, 7], ids=["d3", "d5", "d7"])
def test_lscx_qubit_count_lower_bound(d):
    """Qubit count should be at least 3 patches worth (ancilla, control, target)."""
    circuit = surgery_circuit(
        distance=d,
        target_state_init="X+",
        control_state_init="X+",
        flow_observable="X -> XX",
    )

    assert isinstance(circuit, stim.Circuit)

    # One rotated surface code patch uses 2*d^2 - 1 qubits.
    expected_min_qubits = 3 * (2 * d**2 - 1)

    assert circuit.num_qubits >= expected_min_qubits


@pytest.mark.parametrize("distance", [3, 5], ids=["d3", "d5"])
def test_lscx_invalid_flow_and_state_combinations_return_error(distance):
    """Invalid combinations should raise ValueError."""

    # Invalid flow name raises ValueError
    with pytest.raises(ValueError, match="Invalid Flow selected"):
        surgery_circuit(
            distance=distance,
            target_state_init="X+",
            control_state_init="X+",
            flow_observable="INVALID",
        )

    # Using Y-basis isn't supported - raises ValueError
    with pytest.raises(
        ValueError,
        match="Invalid control/target state initialization",
    ):
        surgery_circuit(
            distance=distance,
            target_state_init="Y+",
            control_state_init="X+",
            flow_observable="X -> X",
        )

    # Wrong basis per flow rules - raises ValueError
    # X -> XX requires both control and target in X basis
    with pytest.raises(ValueError, match="Wrong target basis for selected flow"):
        surgery_circuit(
            distance=distance,
            target_state_init="Z0",
            control_state_init="X+",
            flow_observable="X -> XX",
        )

    # Z -> ZZ requires both control and target in Z basis
    with pytest.raises(ValueError, match="Wrong target basis for selected flow"):
        surgery_circuit(
            distance=distance,
            target_state_init="X+",
            control_state_init="Z0",
            flow_observable="Z -> ZZ",
        )

    # ZX -> ZX requires control in Z basis and target in X basis
    with pytest.raises(ValueError, match="Wrong control basis for selected flow"):
        surgery_circuit(
            distance=distance,
            target_state_init="Z0",
            control_state_init="X+",
            flow_observable="ZX -> ZX",
        )


@pytest.mark.parametrize("distance", [2, 4, 10], ids=["d2", "d4", "d10"])
def test_lscx_even_distance_is_invalid(distance):
    """Even distances are invalid (geometry build requires odd d)."""
    with pytest.raises(ValueError):
        surgery_circuit(
            distance=distance,
            target_state_init="X+",
            control_state_init="X+",
            flow_observable="X -> XX",
        )
