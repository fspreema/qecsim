import pytest
import stim

from src.codes.xzzx.builder import XZZXBuilder


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
def test_xzzx_qubits(distance):
    """Test XZZX code basic compilation"""
    circuit = XZZXBuilder(
        distance=distance,
        state_init="XZZX-VER",
    ).build_circuit()

    assert isinstance(circuit, stim.Circuit)

    # XZZX code: d^2 data + (d-1)^2 interior ancilla
    # + 2*(d-1) boundary = 2d^2 - 1
    expected_qubits = 2 * distance**2 - 1

    assert circuit.num_qubits == expected_qubits


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize("state_init", ["XZZX-VER", "XZZX-HOR"])
def test_xzzx_different_states(distance, state_init):
    """Test XZZX code with different state initializations"""
    circuit = XZZXBuilder(
        distance=distance,
        state_init=state_init,
    ).build_circuit()

    assert isinstance(circuit, stim.Circuit)


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize("state_init", ["XZZX-VER", "XZZX-HOR"])
def test_xzzx_observable_present(distance, state_init):
    """Test that XZZX circuit includes observable"""
    circuit = XZZXBuilder(
        distance=distance,
        state_init=state_init,
    ).build_circuit()
    # Check for observable in circuit
    assert "OBSERVABLE_INCLUDE" in str(circuit)


@pytest.mark.parametrize("state_init", ["XZZX-VER", "XZZX-HOR"], ids=["XZZX-VER", "XZZX-HOR"])
@pytest.mark.parametrize("distance", [2, 4, 6, 8, 10], ids=["d2", "d4", "d6", "d8", "d10"])
def test_xzzx_invalid_params_even_distance(state_init, distance):
    """Test that XZZX code raises error for even distance"""
    # Even distance should fail
    with pytest.raises(ValueError):
        XZZXBuilder(
            distance=distance,
            state_init=state_init,
        ).build_circuit()


@pytest.mark.parametrize("state_init", ["XZZX-VER", "XZZX-HOR"], ids=["XZZX-VER", "XZZX-HOR"])
def test_xzzx_invalid_params_distance_too_small(state_init):
    """Test that XZZX code raises error for distance < 3"""
    # Distance < 3 should fail
    with pytest.raises(ValueError):
        XZZXBuilder(
            distance=1,
            state_init=state_init,
        ).build_circuit()


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
def test_xzzx_invalid_params_invalid_state(distance):
    """Test that XZZX code raises error for invalid state initialization"""
    # Invalid state initialization should fail
    with pytest.raises(ValueError):
        XZZXBuilder(
            distance=distance,
            state_init="Invalid",
        ).build_circuit()
