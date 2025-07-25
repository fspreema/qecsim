# test_threshold_approx.py
import math
import numpy as np
import pytest

from qecsim.threshold_num.calc_threshold import threshold_approx

# helpers.py
import uuid
from collections import Counter
import sinter


import math, uuid
from collections import Counter
import numpy as np
import sinter


# ---------------------------------------------------------------------
# 1 · Build a *real* TaskStats obeying shots ≥ errors + discards
# ---------------------------------------------------------------------
def make_task_stats(*, distance, p, logical_p, shots,
                    strong_id=None, decoder="", seconds=0.0, discards=0):
    errors = int(np.random.binomial(shots, logical_p))
    return sinter.TaskStats(
        strong_id = strong_id or str(uuid.uuid4()),
        decoder   = decoder,
        json_metadata = {"distance": distance, "p": p},
        shots     = shots,
        errors    = errors,
        discards  = discards,
        seconds   = seconds,
        custom_counts = Counter(),
    )


# ---------------------------------------------------------------------
# 2 · Synthetic data that *never* makes logical_p > 1
# ---------------------------------------------------------------------
def synthetic_dataset(
        p_star,                # crossing position (linear space)
        slope1, slope2,        # log‑log slopes (usually negative)
        L_star   = 1e-3,       # logical error‑rate at the crossing
        distances = (5, 7),
        p_window  = (1e-5, 1e-1),
        n_points  = 20,
        shots     = 1_000_000,
        noise     = 0.0,       # σ of Gaussian in *log10*-space
        seed      = 0,
    ):
    """
    Returns a list[TaskStats] whose two distance‑layers
    cross exactly at (p_star, L_star).
    """
    rng = np.random.default_rng(seed)
    data = []

    logL_star = math.log10(L_star)
    logp_star = math.log10(p_star)

    slopes   = (slope1, slope2)
    ps       = np.logspace(math.log10(p_window[0]),
                           math.log10(p_window[1]),
                           n_points)

    for dist, k in zip(distances, slopes):
        # Choose intercept so the lines meet at (p_star, L_star)
        a = logL_star - k * logp_star

        for p in ps:
            logL = a + k * math.log10(p)
            if noise:
                logL += rng.normal(0.0, noise)

            L = max(min(10 ** logL, 0.999), 1e-15)   # clamp to (0,1)
            data.append(make_task_stats(distance=dist, p=p,
                                        logical_p=L, shots=shots))
    return data



def test_clean_linear():
    tgt = 4e-3
    data = synthetic_dataset(tgt, -1.0, -2.0, noise=0)
    est = threshold_approx(data, p_min=1e-3, p_max=2e-2)
    assert abs(est - tgt) / tgt < 0.02


def test_noisy_curved():
    tgt = 6e-3
    data = synthetic_dataset(tgt, -0.9, -2.1,
                             noise=0.07, shots=5000, n_points=25)
    est = threshold_approx(data, p_min=2e-3, p_max=1.5e-2)
    assert abs(est - tgt) / tgt < 0.12


def test_no_crossing():
    # parallel slopes ⇒ never meet
    data = synthetic_dataset(5e-3, -1.0, -1.0)
    with pytest.raises(RuntimeError):
        threshold_approx(data, p_min=1e-3, p_max=2e-2)
