import pytest
import stim

from qecsim.codes.surface_code_rotated.builder import rotated_surface_code


def test_rotated_surface_code_qubits():
    """Test rotated surface code basic compilation"""
    for distance in [3, 5, 7]:
        result = rotated_surface_code(
            distance=distance,
            rounds=distance,
            state_init="0",
            log_obs="Z",
        )
        # When state_init and log_obs are in same basis, returns just circuit
        assert isinstance(result, stim.Circuit)

        # Surface code: d^2 + (d-1)^2 + 2(d-1) = 2d^2 - 1
        expected_qubits = 2 * distance**2 - 1

        assert result.num_qubits == expected_qubits


def test_rotated_surface_code_different_states():
    """Test rotated surface code with different state initializations"""

    for distance in [3, 5, 7]:
        # Test Z-basis states (0, 1) with Z observable
        for state_init in ["0", "1"]:
            result = rotated_surface_code(
                distance=distance,
                rounds=distance,
                state_init=state_init,
                log_obs="Z",
            )
            assert isinstance(result, stim.Circuit)

        # Test X-basis states (+, -) with X observable
        for state_init in ["+", "-"]:
            result = rotated_surface_code(
                distance=distance,
                rounds=distance,
                state_init=state_init,
                log_obs="X",
            )
            assert isinstance(result, stim.Circuit)

        # Test Y-basis states (+i, -i) with Y observable
        for state_init in ["+i", "-i"]:
            result = rotated_surface_code(
                distance=distance,
                rounds=distance,
                state_init=state_init,
                log_obs="Y",
            )
            # Y-basis returns just circuit when log_obs is Y
            assert isinstance(result, stim.Circuit)


def test_rotated_surface_code_mismatched_basis():
    """Test Y-basis states return tuple when observable is different"""
    for distance in [3, 5, 7]:
        # Z-basis states with X or Y observable should return
        # (circuit, obs_indices)
        for state_init in ["0", "1"]:
            for log_obs in ["X", "Y"]:
                result = rotated_surface_code(
                    distance=distance,
                    rounds=distance,
                    state_init=state_init,
                    log_obs=log_obs,
                )

                # Should return a tuple (circuit, obs_indices)
                assert isinstance(result, tuple)
                assert len(result) == 2
                circuit, obs_indices = result
                assert isinstance(circuit, stim.Circuit)
                assert isinstance(obs_indices, list)

        # X-basis states with Z or Y observable should return
        # (circuit, obs_indices)
        for state_init in ["+", "-"]:
            for log_obs in ["Z", "Y"]:
                result = rotated_surface_code(
                    distance=distance,
                    rounds=distance,
                    state_init=state_init,
                    log_obs=log_obs,
                )

                # Should return a tuple (circuit, obs_indices)
                assert isinstance(result, tuple)
                assert len(result) == 2
                circuit, obs_indices = result
                assert isinstance(circuit, stim.Circuit)
                assert isinstance(obs_indices, list)

        # Y-basis states with X or Z observable should return
        # (circuit, obs_indices)
        for state_init in ["+i", "-i"]:
            for log_obs in ["X", "Z"]:
                result = rotated_surface_code(
                    distance=distance,
                    rounds=distance,
                    state_init=state_init,
                    log_obs=log_obs,
                )

                # Should return a tuple (circuit, obs_indices)
                assert isinstance(result, tuple)
                assert len(result) == 2
                circuit, obs_indices = result
                assert isinstance(circuit, stim.Circuit)
                assert isinstance(obs_indices, list)


def test_rotated_surface_code_observable_present():
    """Test rotated surface code circuit includes observable when appropriate"""
    distance = 5

    # OBSERVABLE_INCLUDE should be present
    for state_init in ["0", "1", "+", "-", "+i", "-i"]:
        for log_obs in ["X", "Y", "Z"]:
            result = rotated_surface_code(
                distance=distance,
                rounds=distance,
                state_init=state_init,
                log_obs=log_obs,
            )

        if isinstance(result, tuple):
            circuit, _ = result
            assert "OBSERVABLE_INCLUDE" in str(circuit)


def test_rotated_surface_code_with_logical_h():
    """Test rotated surface code with logical Hadamard"""

    for distance in [3, 5, 7]:
        # Test with Z-basis states
        for state_init in ["0", "1"]:
            result = rotated_surface_code(
                distance=distance,
                rounds=distance,
                state_init=state_init,
                log_obs="Z",
                logical_h=True,
            )
            assert isinstance(result, stim.Circuit)

        # Test with X-basis states
        for state_init in ["+", "-"]:
            result = rotated_surface_code(
                distance=distance,
                rounds=distance,
                state_init=state_init,
                log_obs="X",
                logical_h=True,
            )
            assert isinstance(result, stim.Circuit)


def test_rotated_surface_code_invalid_params():
    """Test that rotated surface code raises error for invalid parameters"""

    for distance in [2, 4, 6, 8, 10]:
        # Even distance should fail
        with pytest.raises(ValueError):
            rotated_surface_code(
                distance=distance,
                rounds=distance,
                state_init="0",
                log_obs="Z",
            )

        # invalid state_init should fail
        with pytest.raises(ValueError):
            rotated_surface_code(
                distance=distance,
                rounds=distance,
                state_init="Invalid",
                log_obs="Invalid",
            )

    # Distance < 3 should fail
    with pytest.raises(ValueError):
        rotated_surface_code(
            distance=1,
            rounds=1,
            state_init="0",
            log_obs="Z",
        )
