# tests/core/test_geometry.py

import pytest

from qecsim.core.geometry import build_lattice


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
def test_build_lattice_basic(distance):
    """Test basic lattice building returns correct structure."""

    lattice = build_lattice(distance=distance)

    assert isinstance(lattice, dict)
    # All keys should be complex coordinates
    assert all(isinstance(k, complex) for k in lattice.keys())
    # All values should be strings (labels)
    assert all(isinstance(v, str) for v in lattice.values())


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
def test_build_lattice_qubit_count(distance):
    """Test that lattice has correct number of qubits for various distances."""

    lattice = build_lattice(distance=distance)

    # Count data qubits and internal stabilizers
    data_qubits = sum(1 for v in lattice.values() if "DATA" in v)
    stabs = sum(1 for v in lattice.values() if "STAB" in v)

    # For distance d, we expect d^2 data qubits
    assert data_qubits == distance**2
    # Interior stabilizers should be (d-1)^2
    assert stabs == (distance - 1) ** 2


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
def test_build_lattice_state_norm(distance):
    """Test lattice with no state initialization."""
    lattice = build_lattice(distance=distance, state_init=None)

    # Should have generic DATA labels
    assert any(v == "DATA" for v in lattice.values())
    # Should have X-STAB and Z-STAB
    assert any(v == "X-STAB" for v in lattice.values())
    assert any(v == "Z-STAB" for v in lattice.values())


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize("state_init", ["Ver", "Hor"], ids=["Ver", "Hor"])
def test_build_lattice_state_init_xzzx(distance, state_init):
    """Test lattice with vertical state initialization."""
    lattice = build_lattice(distance=distance, state_init=state_init)

    # Should have DATA_X and DATA_Z labels
    assert any(v == "DATA_X" for v in lattice.values())
    assert any(v == "DATA_Z" for v in lattice.values())
    # Should have STAB-Ver and STAB-Hor
    assert any(v == "STAB-Ver" for v in lattice.values())
    assert any(v == "STAB-Hor" for v in lattice.values())


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize(
    "offset",
    [0 + 0j, 1 + 2j, -2 + 3j],
    ids=["offset_0", "offset_1_2", "offset_-2_3"],
)
def test_build_lattice_offset(distance, offset):
    """Offset should translate all coordinates without changing labels."""

    base = build_lattice(distance=distance, offset=0 + 0j)
    shifted = build_lattice(distance=distance, offset=offset)

    # The translated base must exactly equal the shifted lattice
    translated = {
        complex(
            int(k.real) + int(offset.real),
            int(k.imag) + int(offset.imag),
        ): v
        for k, v in base.items()
    }

    assert len(translated) == len(shifted)
    assert translated == shifted


@pytest.mark.parametrize("distance", [3, 5, 7], ids=["d3", "d5", "d7"])
@pytest.mark.parametrize("starting_stabilizer_x", [True, False], ids=["start_x", "start_z"])
@pytest.mark.parametrize("state_init", [None, "Ver"], ids=["none", "ver"])
def test_build_lattice_starting_stabilizer_general(distance, starting_stabilizer_x, state_init):
    """All interior stabilizers follow the alternating pattern for any distance.

    Check every even-even interior coordinate (excluding the 0-axis)
    -> verify label matches the expected row/column alternating scheme.
    """
    lattice = build_lattice(
        distance=distance,
        starting_stabilizer_x=starting_stabilizer_x,
        state_init=state_init,
        offset=0 + 0j,
    )

    def expected_label(real: int, imag: int) -> str:
        """Expected stabilizer label at (real, imag) for given start/initialization.

        Interior stabilizers form a checkerboard. If starting_stabilizer_x is True,
        the top-left interior stabilizer (2,2) is X; otherwise it's Z. The pattern
        alternates by column and flips each interior row.
        """
        # Zero-based interior row/col indices: 2,4,...,2d-2 -> 0,1,...,d-2
        row_idx = (real - 2) // 2
        col_idx = (imag - 2) // 2

        # Checkerboard parity across interior grid -> One type odd, the other even
        parity = (row_idx + col_idx) % 2
        use_x = (parity == 0) if starting_stabilizer_x else (parity == 1)

        if state_init is None:
            label_x, label_z = "X-STAB", "Z-STAB"
        else:
            label_x, label_z = "STAB-Ver", "STAB-Hor"
        return label_x if use_x else label_z

    # Run through all lattice coords and check validity
    for real in range(2, distance * 2, 2):
        if real == 0:
            continue
        for imag in range(2, distance * 2, 2):
            if imag == 0:
                continue
            coord = complex(real, imag)
            assert coord in lattice
            assert lattice[coord] == expected_label(real, imag)


def test_build_lattice_coordinate_positions():
    """Test that data qubits and stabilizers are at correct positions."""
    lattice = build_lattice(distance=3, offset=0 + 0j)

    # Data qubits should be at odd real, odd imag
    for coord, label in lattice.items():
        if "DATA" in label:
            assert int(coord.real) % 2 == 1
            assert int(coord.imag) % 2 == 1

    # Interior stabilizers should be at even real, even imag
    for coord, label in lattice.items():
        if "STAB" in label and "BOUND" not in label:
            assert int(coord.real) % 2 == 0
            assert int(coord.imag) % 2 == 0
            assert int(coord.real) > 0
            assert int(coord.imag) > 0


@pytest.mark.parametrize("distance_even", [2, 4, 10], ids=["d2", "d4", "d10"])
def test_build_lattice_invalid_distance_even(distance_even):
    """Test that even distances raise ValueError."""
    with pytest.raises(ValueError, match="distance must be odd"):
        build_lattice(distance=distance_even)


def test_build_lattice_invalid_distance_too_small():
    """Test that distance <= 2 raises ValueError."""
    for distance in [1, 2]:
        with pytest.raises(ValueError, match="distance must be odd"):
            build_lattice(distance=distance)


def test_build_lattice_invalid_state_init():
    """Test that invalid state_init raises ValueError."""
    invalid_states = ["X", "Z", "Invalid", "x", "ver", "hor", ""]

    for state in invalid_states:
        with pytest.raises(ValueError, match="state_init must be either"):
            build_lattice(distance=3, state_init=state)
