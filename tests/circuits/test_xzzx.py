import pytest
import stim

from qecsim.codes.xzzx.circuit import xzzx_code


def test_xzzx_qubits():
    """Test XZZX code basic compilation"""
    for distance in [3, 5, 7]:
        circuit = xzzx_code(
            distance=distance,
            rounds=distance,
            state_init="Ver",
        )
        assert isinstance(circuit, stim.Circuit)

        # XZZX code: d^2 data + (d-1)^2 interior ancilla
        # + 2*(d-1) boundary = 2d^2 - 1
        expected_qubits = 2 * distance**2 - 1

        assert circuit.num_qubits == expected_qubits


def test_xzzx_different_states():
    """Test XZZX code with different state initializations"""

    for distance in [3, 5, 7]:
        for state_init in ["Ver", "Hor"]:
            circuit = xzzx_code(
                distance=distance,
                rounds=distance,
                state_init=state_init,
            )

            assert isinstance(circuit, stim.Circuit)


def test_xzzx_observable_present():
    """Test that XZZX circuit includes observable"""
    for distance in [3, 5, 7]:
        for state_init in ["Ver", "Hor"]:
            circuit = xzzx_code(
                distance=distance,
                rounds=distance,
                state_init=state_init,
            )

            # Check for observable in circuit
            assert "OBSERVABLE_INCLUDE" in str(circuit)


def test_xzzx_invalid_params():
    """Test that XZZX code raises error for invalid parameters"""

    for state_init in ["Ver", "Hor"]:
        for distance in [2, 4, 6, 8, 10]:
            # Even distance should fail
            with pytest.raises(ValueError):
                xzzx_code(
                    distance=distance,
                    rounds=distance,
                    state_init=state_init,
                )

        # Distance < 3 should fail
        with pytest.raises(ValueError):
            xzzx_code(
                distance=1,
                rounds=1,
                state_init=state_init,
            )

        for distance in [3, 5, 7]:
            # Invalid state initialization should fail
            with pytest.raises(ValueError):
                xzzx_code(
                    distance=distance,
                    rounds=distance,
                    state_init="Invalid",
                )
