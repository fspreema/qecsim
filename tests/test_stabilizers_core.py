from qecsim.core.stabilizers import populate_stab_to_data

# Simple patches to exercise the three styles

def test_xzzx_basic():
    patch = {
        1+1j: "STAB-Ver",
        3+1j: "STAB-Hor",
        0+2j: "STAB-BOUND-L-Hor",
    }
    m = populate_stab_to_data(patch)
    # Expect keys to be present for neighbours of 1+1j and 3+1j
    assert any(k[1] == 1+1j for k in m.keys())
    assert any(k[1] == 3+1j for k in m.keys())


def test_surface_basic():
    patch = {
        1+1j: "X-STAB",
        3+1j: "Z-STAB",
        2+2j: "X-STAB",
    }
    # default surface call
    m = populate_stab_to_data(patch, is_flipped=False, y_basis=False, y_switch=False, y_memory=False, distance=2, offset=0+0j)
    # expect entries for these stabs
    assert any(v.startswith('1-') or v.startswith('2-') for v in m.values())
    assert any(k[1] == 1+1j or k[1] == 2+2j for k in m.keys())


def test_lattice_surgery_basic():
    patch = {
        1+1j: "X-STAB",
        3+1j: "Z-STAB",
        0+1j: "X-STAB-BOUND-B-A",
    }
    m = populate_stab_to_data(patch, True, "AC")
    assert isinstance(m, dict)
    # expect at least one CX order string
    assert any(isinstance(v, str) and 'CX' in v for v in m.values())
