import gc
import glob
import os

import joblib
import numpy as np
import sinter

from qecsim.lattice_surgery.circuit import surgery_circuit

"""
This is super INEFFICIENT!

-> For high d and really low prob for phys error, max_shots are to low to detect any error (Can be skipped)
-> Combined files are saved thorugh loop which loops until d_max = 15 (Have files we don't need in form of _12_ _14_)
-> Only save the run datasets where 2 distances are included -> Every other doesn't need to be calculated

TODO:
-> Rewrite sim run with different noise for different distances
-> Change saving loop
-> Check if distance pairs are present before running sim run
"""

if __name__ == "__main__":

    d_max = 15

    for d in [i for i in range(3,d_max,2)]:

        for r in [i for i in range(2,d,2)]:

            ###################
            # Calc d thresholds
            ###################

            task_noisy_xx = [sinter.Task(
                circuit = surgery_circuit(
                    distance = d,
                    round_num = r,
                    round_split = r,
                    round_merge = r,
                    target_state_init ="X+",
                    control_state_init ="X+",
                    flow_observable = "X -> XX",
                    noise_after_clifford_depol = noise,
                    noise_measure_flip = noise,
                    noise_after_reset = noise,
                    noise_depol_data_init = noise,
                    ),
                json_metadata={"p": noise, "distance" : d, "rounds" : r},
                )
                for noise in [i for i in np.arange(1e-5, 0.013, 5e-5)]
            ]

            stats_noisy_xx : list[sinter.TaskStats] = sinter.collect(
                num_workers = os.cpu_count(),
                tasks = task_noisy_xx,
                decoders = ["pymatching"],
                max_shots = 1_000_000,
                max_errors = 10_000,
                print_progress = True,
            )

            #Saving Stats
            fname = f"lscx_singlerounds_d{d}_r{r}.pkl"
            try:
                joblib.dump(stats_noisy_xx, fname, compress=3)
                print(f"Successfully saved {fname}")
            except Exception as e:
                print(f"Failed to save {fname}: {e}")

            # Clear RAM by deleting large variables and running garbage collection
            del stats_noisy_xx
            del task_noisy_xx
            gc.collect()

    ###########################
    # Adding all files together
    ###########################

    def combine_files(round_num):

        pattern = f"lscx_singlerounds_d*_r{round_num}.pkl"
        files = sorted(glob.glob(pattern))
        combined_stats = []
        for f in files:
            part = joblib.load(f)
            combined_stats.extend(part)
        out_name = f"lscx_singlerounds_{round_num}_combined.pkl"
        joblib.dump(combined_stats, out_name, compress=3)
        print(f"Combined and saved as {out_name}")

    for r in [i for i in range(2,d_max,2)]:
        combine_files(r)

    ######################################
    # Delete all individual distance files
    ######################################

    for fname in glob.glob("lscx_*_d*.pkl"):
        try:
            os.remove(fname)
            print(f"Deleted {fname}")
        except Exception as e:
            print(f"Failed to delete {fname}: {e}")
