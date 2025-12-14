import itertools
import os
import pickle

import numpy as np
import sinter

from src.codes.xzzx.builder import xzzx_code
from src.tools.thresholds.threshold_approx import threshold_approx

if __name__ == "__main__":
    # -----------------------------Calc----------------------------

    # Initlize Saving Array
    num_value: list = []
    step = 0.0125
    bias_steps = np.arange(0.0, 1 + (step / 2), step)
    all_bias_triplets = [
        [bx, by, 1 - bx - by]
        for bx, by in itertools.product(bias_steps, repeat=2)
        if 0 <= 1 - bx - by <= 1
    ]

    # Changing Bias setting
    for bx, by, bz in all_bias_triplets:
        # Creating Current Bias
        current_bias = [bx, by, bz]

        # Running simulation
        task_noisy_v = [
            sinter.Task(
                circuit=xzzx_code(
                    distance=d,
                    rounds=d,
                    state_init="Ver",
                    noise_bias=current_bias,
                    after_c_pauli_channel_prob=noise,
                ),
                json_metadata={"p": noise, "distance": d, "bias": current_bias},
            )
            for noise in [i for i in np.arange(0.005, 0.1, 0.005)]
            for d in [5, 7]
        ]

        stats_noisy_v: list[sinter.TaskStats] = sinter.collect(
            num_workers=os.cpu_count(),
            tasks=task_noisy_v,
            decoders=["pymatching"],
            max_shots=500_000,
            max_errors=10_000,
            print_progress=True,
        )

        """
        Drop in for implementation of near PM decoder
        -> Really long runtime
        """

        """
        stats_noisy_v : list[sinter.TaskStats] = sinter.collect(
            num_workers = os.cpu_count(),
            tasks=task_noisy_v,
            decoders=['bposd'],
            custom_decoders=sinter_decoders(),
            max_shots=500_000,
            max_errors=10_000,
            print_progress=True
        )"""

        # calculate threshold and save -> If Error skip and set 0
        try:
            calc_th = threshold_approx(stats_noisy_v)
            num_value.append([current_bias, calc_th])

        except Exception:
            num_value.append([current_bias, -1])

    # Saving num_values with pickle
    with open("XZZX_num_value(full_bias_low_d).pkl", "wb") as file:
        pickle.dump(num_value, file)
