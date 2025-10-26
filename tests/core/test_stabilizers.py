from qecsim.core.stabilizers import populate_stab_to_data


def get_weight(mapping: dict, stab_coords: dict) -> dict[complex, int]:
    """Compute how many data neighbors each stabilizer has from the returned mapping.

    Function returns dict with coord tuple and implementation order as str.
    -> Count how many times each Stabilizer coord is appearing
    -> Boundary Stabs should have weight 2, interior weight 4
    """

    deg: dict[complex, int] = {}

    for k in mapping.keys():
        # Check both positions in the key tuple
        if k[1] in stab_coords:
            stab = k[1]
            deg[stab] = deg.get(stab, 0) + 1
        elif k[0] in stab_coords:
            stab = k[0]
            deg[stab] = deg.get(stab, 0) + 1

    return deg


def test_xzzx_basic():
    patch = {
        1 + 1j: "STAB-Ver",
        3 + 1j: "STAB-Hor",
        0 + 2j: "STAB-BOUND-L-Hor",
    }

    m = populate_stab_to_data(patch)

    # Expect keys to be present for neighbours of 1+1j and 3+1j
    assert any(k[1] == 1 + 1j for k in m.keys())
    assert any(k[1] == 3 + 1j for k in m.keys())

    # Weight checks via black-box counting
    deg = get_weight(m, patch)
    assert deg.get(1 + 1j, 0) == 4
    assert deg.get(3 + 1j, 0) == 4
    assert deg.get(0 + 2j, 0) == 2


def test_surface_basic():
    patch = {
        1 + 1j: "X-STAB",
        3 + 1j: "Z-STAB",
        2 + 2j: "X-STAB",
    }
    # default surface call
    m = populate_stab_to_data(
        patch,
        is_flipped=False,
        y_basis=False,
        y_switch=False,
        y_memory=False,
        distance=2,
        offset=0 + 0j,
    )
    # expect entries for these stabs
    assert any(v.startswith("1-") or v.startswith("2-") for v in m.values())
    assert any(k[1] == 1 + 1j or k[1] == 2 + 2j for k in m.keys())

    # Interior stabilizer should have weight 4
    deg = get_weight(m, patch)
    assert deg.get(2 + 2j, 0) == 4


def test_lattice_surgery_basic():
    patch = {
        1 + 1j: "X-STAB",
        3 + 1j: "Z-STAB",
        0 + 1j: "X-STAB-BOUND-B-A",
    }
    m = populate_stab_to_data(patch, merging=True, merging_type="AT")

    assert isinstance(m, dict)
    # expect at least one CX order string
    assert any(isinstance(v, str) and "CX" in v for v in m.values())

    # Degrees of any stabilizers that appear should be small and valid
    deg = get_weight(m, patch)

    assert deg.get(1 + 1j, 0) == 4
    assert deg.get(3 + 1j, 0) == 4
    assert deg.get(0 + 1j, 0) == 2
