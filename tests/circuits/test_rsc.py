import pytest
import stim

from src.codes.surface_code_rotated.builder import rotated_surface_code


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
def test_rotated_surface_code_qubits(distance):
    """Test rotated surface code basic compilation"""
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


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize(
    "state_init,log_obs",
    [
        # Z-basis states (0, 1) with Z observable
        ("0", "Z"),
        ("1", "Z"),
        # X-basis states (+, -) with X observable
        ("+", "X"),
        ("-", "X"),
        # Y-basis states (+i, -i) with Y observable
        ("+i", "Y"),
        ("-i", "Y"),
    ],
)
def test_rotated_surface_code_different_states(distance, state_init, log_obs):
    """Test rotated surface code with different state initializations"""
    result = rotated_surface_code(
        distance=distance,
        rounds=distance,
        state_init=state_init,
        log_obs=log_obs,
    )
    assert isinstance(result, stim.Circuit)


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize(
    "state_init,log_obs",
    [
        # Z-basis states with X or Y observable
        ("0", "X"),
        ("0", "Y"),
        ("1", "X"),
        ("1", "Y"),
        # X-basis states with Z or Y observable
        ("+", "Z"),
        ("+", "Y"),
        ("-", "Z"),
        ("-", "Y"),
        # Y-basis states with X or Z observable
        ("+i", "X"),
        ("+i", "Z"),
        ("-i", "X"),
        ("-i", "Z"),
    ],
)
def test_rotated_surface_code_mismatched_basis(distance, state_init, log_obs):
    """Test Y-basis states return tuple when observable is different"""
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


@pytest.mark.parametrize(
    "state_init,log_obs",
    [
        ("0", "X"),
        ("0", "Y"),
        ("0", "Z"),
        ("1", "X"),
        ("1", "Y"),
        ("1", "Z"),
        ("+", "X"),
        ("+", "Y"),
        ("+", "Z"),
        ("-", "X"),
        ("-", "Y"),
        ("-", "Z"),
        ("+i", "X"),
        ("+i", "Y"),
        ("+i", "Z"),
        ("-i", "X"),
        ("-i", "Y"),
        ("-i", "Z"),
    ],
)
def test_rotated_surface_code_observable_present(state_init, log_obs):
    """Test rotated surface code circuit includes observable when appropriate"""
    distance = 5

    result = rotated_surface_code(
        distance=distance,
        rounds=distance,
        state_init=state_init,
        log_obs=log_obs,
    )

    if isinstance(result, tuple):
        circuit, meas = result
        assert "OBSERVABLE_INCLUDE" in str(circuit)
        assert len(meas) >= 1


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize(
    "state_init,log_obs",
    [
        # Z-basis states
        ("0", "Z"),
        ("1", "Z"),
        # X-basis states
        ("+", "X"),
        ("-", "X"),
    ],
)
def test_rotated_surface_code_with_logical_h(distance, state_init, log_obs):
    """Test rotated surface code with logical Hadamard"""
    result = rotated_surface_code(
        distance=distance,
        rounds=distance,
        state_init=state_init,
        log_obs=log_obs,
        logical_h=True,
    )
    assert isinstance(result, stim.Circuit)


@pytest.mark.parametrize("distance", [2, 4, 10], ids=["d2", "d4", "d10"])
def test_rotated_surface_code_invalid_params_even_distance(distance):
    """Test that rotated surface code raises error for even distance"""
    # Even distance should fail
    with pytest.raises(ValueError):
        rotated_surface_code(
            distance=distance,
            rounds=distance,
            state_init="0",
            log_obs="Z",
        )


@pytest.mark.parametrize("distance", [2, 4, 10], ids=["d2", "d4", "d10"])
def test_rotated_surface_code_invalid_params_invalid_state(distance):
    """Test that rotated surface code raises error for invalid state_init"""
    # invalid state_init should fail
    with pytest.raises(ValueError):
        rotated_surface_code(
            distance=distance,
            rounds=distance,
            state_init="Invalid",
            log_obs="Invalid",
        )


def test_rotated_surface_code_invalid_params_distance_too_small():
    """Test that rotated surface code raises error for distance < 3"""
    # Distance < 3 should fail
    with pytest.raises(ValueError):
        rotated_surface_code(
            distance=1,
            rounds=1,
            state_init="0",
            log_obs="Z",
        )
