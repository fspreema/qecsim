from qecsim.core.stabilizers import populate_stab_to_data
from qecsim.codes import surface_code_rotated as surf_mod
from qecsim.codes import lattice_surgery as ls_mod


def _exercise_surface_cases(patch, **kwargs):
    core = populate_stab_to_data(patch, **kwargs)
    per = surf_mod.stabilizers.populate_stab_to_data(patch, **kwargs)
    assert isinstance(core, dict) or isinstance(core, tuple)
    assert isinstance(per, dict) or isinstance(per, tuple)
    # Normalize tuples (core may return (dict, dict) when y_switch is True)
    core_dict = core[0] if isinstance(core, tuple) else core
    per_dict = per[0] if isinstance(per, tuple) else per
    for k, v in per_dict.items():
        assert k in core_dict
        # only compare string values when both are strings
        if isinstance(v, str):
            assert core_dict[k] == v


def test_surface_various_modes():
    base_patch = {
        1+1j: "X-STAB",
        3+1j: "Z-STAB",
        2+2j: "X-STAB",
    }

    # default
    _exercise_surface_cases(base_patch, is_flipped=False, y_basis=False, y_switch=False, y_memory=False, distance=2, offset=0+0j)

    # flipped
    _exercise_surface_cases(base_patch, is_flipped=True, y_basis=False, y_switch=False, y_memory=False, distance=2, offset=0+0j)

    # y_basis
    _exercise_surface_cases(base_patch, is_flipped=False, y_basis=True, y_switch=False, y_memory=False, distance=2, offset=0+0j)

    # y_switch (returns tuple)
    _exercise_surface_cases(base_patch, is_flipped=False, y_basis=False, y_switch=True, y_memory=False, distance=3, offset=0+0j)


def test_surface_boundary_label_coverage():
    patch = {
        0+1j: "Z-STAB-BOUND-L",
        4+1j: "Z-STAB-BOUND-R",
        2+0j: "X-STAB-BOUND-B",
        2+4j: "X-STAB-BOUND-U",
        0+0j: "Z-STAB-BOUND-L-T",
        4+0j: "Z-STAB-BOUND-R-T",
        0+2j: "X-STAB-BOUND-A-T",
        4+2j: "X-STAB-BOUND-B-T",
    }
    _exercise_surface_cases(patch, is_flipped=False, y_basis=True, y_switch=True, y_memory=False, distance=3, offset=0+0j)


def test_lattice_surgery_many_merging_types():
    patch = {
        1+1j: "X-STAB",
        3+1j: "Z-STAB",
        0+1j: "X-STAB-BOUND-B-A",
        0+0j: "Z-STAB-BOUND-L-A",
    }
    core_ac = populate_stab_to_data(patch, True, "AC")
    per_ac = ls_mod.stabilizers.populate_stab_to_data(patch, True, "AC")
    assert isinstance(core_ac, dict) and isinstance(per_ac, dict)
    for k, v in per_ac.items():
        assert k in core_ac

    core_at = populate_stab_to_data(patch, True, "AT")
    per_at = ls_mod.stabilizers.populate_stab_to_data(patch, True, "AT")
    for k, v in per_at.items():
        assert k in core_at


def test_exhaustive_label_coverage_surface_and_ls():
    # All labels used across per-code stabilizer modules — exercise each once
    labels = [
        "X-STAB",
        "X-STAB-BOUND-A-A",
        "X-STAB-BOUND-A-C",
        "X-STAB-BOUND-A-T",
        "X-STAB-BOUND-B",
        "X-STAB-BOUND-B-A",
        "X-STAB-BOUND-B-C",
        "X-STAB-BOUND-B-T",
        "X-STAB-BOUND-R",
        "X-STAB-BOUND-R-H",
        "X-STAB-BOUND-U",
        "X-STAB-SURGERY-B",
        "X-STAB-SURGERY-M",
        "Z-STAB",
        "Z-STAB-BOUND-L",
        "Z-STAB-BOUND-L-A",
        "Z-STAB-BOUND-L-C",
        "Z-STAB-BOUND-L-T",
        "Z-STAB-BOUND-R",
        "Z-STAB-BOUND-R-A",
        "Z-STAB-BOUND-R-C",
        "Z-STAB-BOUND-R-T",
        "Z-STAB-BOUND-U",
        "Z-STAB-BOUND-U-H",
        "Z-STAB-SURGERY-L",
        "Z-STAB-SURGERY-M",
    ]

    # place each label at a distinct coordinate so their neighbor calculations don't collide
    patch = {complex(i, i): lab for i, lab in enumerate(labels, start=1)}

    # surface variations
    core = populate_stab_to_data(patch, is_flipped=False, y_basis=True, y_switch=True, y_memory=False, distance=5, offset=0+0j)
    per = surf_mod.stabilizers.populate_stab_to_data(patch, is_flipped=False, y_basis=True, y_switch=True, y_memory=False, distance=5, offset=0+0j)
    core_dict = core[0] if isinstance(core, tuple) else core
    per_dict = per[0] if isinstance(per, tuple) else per
    for k, v in per_dict.items():
        assert k in core_dict

    # lattice surgery for labels that the lattice surgery module recognizes
    ls_labels = [l for l in labels if l.startswith('X-STAB') or l.startswith('Z-STAB')]
    ls_patch = {complex(i, i): lab for i, lab in enumerate(ls_labels, start=1)}
    core_ls = populate_stab_to_data(ls_patch, True, "AC")
    per_ls = ls_mod.stabilizers.populate_stab_to_data(ls_patch, True, "AC")
    for k in per_ls.keys():
        assert k in core_ls
