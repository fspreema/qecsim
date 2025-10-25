from qecsim.core.stabilizers import populate_stab_to_data
from qecsim.codes.surface_code_rotated import stabilizers as surf_stab
from qecsim.codes.lattice_surgery import stabilizers as ls_stab


def test_surface_equivalence_simple():
    # small patch used in existing tests
    patch = {
        1+1j: "X-STAB",
        3+1j: "Z-STAB",
        2+2j: "X-STAB",
    }
    core = populate_stab_to_data(patch, is_flipped=False, y_basis=False, y_switch=False, y_memory=False, distance=2, offset=0+0j)
    # call the per-code implementation directly (old behaviour)
    per_code = surf_stab.populate_stab_to_data(patch, is_flipped=False, y_basis=False, y_switch=False, y_memory=False, distance=2, offset=0+0j)

    assert isinstance(core, dict)
    assert isinstance(per_code, dict)
    # keys produced by per-code implementation should be present in core output
    for k in per_code.keys():
        assert k in core


def test_lattice_surgery_equivalence_simple():
    patch = {
        1+1j: "X-STAB",
        3+1j: "Z-STAB",
        0+1j: "X-STAB-BOUND-B-A",
    }
    core = populate_stab_to_data(patch, True, "AC")
    per_code = ls_stab.populate_stab_to_data(patch, True, "AC")

    assert isinstance(core, dict)
    assert isinstance(per_code, dict)
    for k in per_code.keys():
        assert k in core


def test_surface_boundary_cases():
    # exercise a few boundary labels used in surface code module
    patch = {
        0+1j: "Z-STAB-BOUND-L",
        4+1j: "Z-STAB-BOUND-R",
        2+0j: "X-STAB-BOUND-B",
        2+4j: "X-STAB-BOUND-U",
    }
    core = populate_stab_to_data(patch, is_flipped=False, y_basis=False, y_switch=False, y_memory=False, distance=3, offset=0+0j)
    per_code = surf_stab.populate_stab_to_data(patch, is_flipped=False, y_basis=False, y_switch=False, y_memory=False, distance=3, offset=0+0j)

    assert isinstance(core, dict)
    assert isinstance(per_code, dict)
    for k, v in per_code.items():
        assert k in core
        assert core[k] == v
