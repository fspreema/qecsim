import pytest
import stim

from qecsim.codes.lattice_surgery.builder import surgery_circuit


def _has_observable(circuit: stim.Circuit) -> bool:
    return "OBSERVABLE_INCLUDE" in str(circuit)


def test_lscx_valid_flows_compile_and_have_observable():
    """Each supported flow compiles for valid state combos and includes observable."""
    distances = [3, 5]

    valid_cases = [
        # All-X
        ("X -> XX", {"control": ["X+", "X-"], "target": ["X+", "X-"]}),
        ("XX -> X", {"control": ["X+", "X-"], "target": ["X+", "X-"]}),
        ("X -> X", {"control": ["X+", "X-", "Z0", "Z1"], "target": ["X+", "X-"]}),
        # All-Z
        ("Z -> ZZ", {"control": ["Z0", "Z1"], "target": ["Z0", "Z1"]}),
        ("ZZ -> Z", {"control": ["Z0", "Z1"], "target": ["Z0", "Z1"]}),
        ("Z -> Z", {"control": ["Z0", "Z1"], "target": ["Z0", "Z1", "X+", "X-"]}),
        # Mixed ZX
        ("ZX -> ZX", {"control": ["Z0", "Z1"], "target": ["X+", "X-"]}),
    ]

    for d in distances:
        for flow, roles in valid_cases:
            # test a couple of specific combos for speed
            for control in roles["control"][:2]:
                for target in roles["target"][:2]:
                    circuit = surgery_circuit(
                        distance=d,
                        target_state_init=target,
                        control_state_init=control,
                        flow_observable=flow,
                    )

                    assert isinstance(circuit, stim.Circuit)
                    assert _has_observable(circuit)


def test_lscx_qubit_count_lower_bound():
    """Qubit count should be at least 3 patches worth (ancilla, control, target)."""
    for d in [3, 5, 7]:
        circuit = surgery_circuit(
            distance=d,
            target_state_init="X+",
            control_state_init="X+",
            flow_observable="X -> XX",
        )

        assert isinstance(circuit, stim.Circuit)
        # One rotated surface code patch uses 2*d^2 - 1 qubits.
        # We have 3 patches plus surgery extras.
        expected_min_qubits = 3 * (2 * d**2 - 1)
        assert circuit.num_qubits >= expected_min_qubits


def test_lscx_invalid_flow_and_state_combinations_return_error():
    """Invalid combinations should raise ValueError."""
    d = 3

    # Invalid flow name raises ValueError
    with pytest.raises(ValueError, match="Invalid Flow selected"):
        surgery_circuit(
            distance=d,
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
            distance=d,
            target_state_init="Y+",
            control_state_init="X+",
            flow_observable="X -> X",
        )

    # Wrong basis per flow rules - raises ValueError
    # X -> XX requires both control and target in X basis
    with pytest.raises(ValueError, match="Wrong target basis for selected flow"):
        surgery_circuit(
            distance=d,
            target_state_init="Z0",
            control_state_init="X+",
            flow_observable="X -> XX",
        )

    # Z -> ZZ requires both control and target in Z basis
    with pytest.raises(ValueError, match="Wrong target basis for selected flow"):
        surgery_circuit(
            distance=d,
            target_state_init="X+",
            control_state_init="Z0",
            flow_observable="Z -> ZZ",
        )

    # ZX -> ZX requires control in Z basis and target in X basis
    with pytest.raises(ValueError, match="Wrong control basis for selected flow"):
        surgery_circuit(
            distance=d,
            target_state_init="Z0",
            control_state_init="X+",
            flow_observable="ZX -> ZX",
        )


def test_lscx_even_distance_is_invalid():
    """Even distances are invalid (geometry build requires odd d)."""

    for distance in [2, 4, 6, 8]:
        with pytest.raises(ValueError):
            surgery_circuit(
                distance=distance,
                target_state_init="X+",
                control_state_init="X+",
                flow_observable="X -> XX",
            )
