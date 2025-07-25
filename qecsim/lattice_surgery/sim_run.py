import sinter
import os
import gc
import joblib
import numpy as np
import glob

from qecsim.lattice_surgery.circuit import surgery_circuit

if __name__ == "__main__":

    for d in [i for i in range(3,15,2)]:

        ###################
        # Calc d thresholds
        ###################
        """
        task_noisy_XX = [sinter.Task(
        circuit = surgery_circuit(
            distance = d,
            round_split = d,
            target_state_init ="X+", 
            control_state_init ="X+",
            flow_observable = "X -> XX",
            noise_after_clifford_depol = noise,
            noise_measure_flip = noise,
            noise_after_reset = noise,
            noise_depol_data_init = noise
            ),
        json_metadata={'p': noise, 'distance' : d, 'rounds' : d}
        ) 
        for noise in [i for i in np.arange(0.004, 0.013, 0.001)]
        ]

        stats_noisy_XX : list[sinter.TaskStats] = sinter.collect(
            num_workers = os.cpu_count(),
            tasks=task_noisy_XX,
            decoders=['pymatching'],
            max_shots=1_000_000,
            max_errors=10_000,
            print_progress=True
        )

        #Saving Stats
        fname = f"lscx_singlerounds_d{d}.pkl"
        try:
            joblib.dump(stats_noisy_XX, fname, compress=3)
            print(f"Successfully saved {fname}")
        except Exception as e:
            print(f"Failed to save {fname}: {e}")

        # Clear RAM by deleting large variables and running garbage collection
        del stats_noisy_XX
        del task_noisy_XX
        gc.collect()
        """

        #######################
        # Calc d * 2 thresholds
        #######################

        task_noisy_XX = [sinter.Task(
        circuit = surgery_circuit(
            distance = d,
            round_split = d * 2,
            target_state_init ="X+", 
            control_state_init ="X+",
            flow_observable = "X -> XX",
            noise_after_clifford_depol = noise,
            noise_measure_flip = noise,
            noise_after_reset = noise,
            noise_depol_data_init = noise
            ),
        json_metadata={'p': noise, 'distance' : d, 'rounds' : d * 2}
        ) 
        for noise in [i for i in np.arange(0.004, 0.013, 0.001)]
        ]

        stats_noisy_XX : list[sinter.TaskStats] = sinter.collect(
            num_workers = os.cpu_count(),
            tasks=task_noisy_XX,
            decoders=['pymatching'],
            max_shots=1_000_000,
            max_errors=10_000,
            print_progress=True
        )

        #Saving Stats
        fname = f"lscx_doublerounds_split_d{d}.pkl"
        try:
            joblib.dump(stats_noisy_XX, fname, compress=3)
            print(f"Successfully saved {fname}")
        except Exception as e:
            print(f"Failed to save {fname}: {e}")

        # Clear RAM by deleting large variables and running garbage collection
        del stats_noisy_XX
        del task_noisy_XX
        gc.collect()

        #########################
        # Calc d * 2/3 thresholds
        #########################

        task_noisy_XX = [sinter.Task(
        circuit = surgery_circuit(
            distance = d,
            round_split = int(d * 2/3),
            target_state_init ="X+", 
            control_state_init ="X+",
            flow_observable = "X -> XX",
            noise_after_clifford_depol = noise,
            noise_measure_flip = noise,
            noise_after_reset = noise,
            noise_depol_data_init = noise
            ),
        json_metadata={'p': noise, 'distance' : d, 'rounds' : int(d * 2/3)}
        ) 
        for noise in [i for i in np.arange(0.004, 0.013, 0.001)]
        ]

        stats_noisy_XX : list[sinter.TaskStats] = sinter.collect(
            num_workers = os.cpu_count(),
            tasks=task_noisy_XX,
            decoders=['pymatching'],
            max_shots=1_000_000,
            max_errors=10_000,
            print_progress=True
        )

        #Saving Stats
        fname = f"lscx_twothirdrounds_split_d{d}.pkl"
        try:
            joblib.dump(stats_noisy_XX, fname, compress=3)
            print(f"Successfully saved {fname}")
        except Exception as e:
            print(f"Failed to save {fname}: {e}")

        # Clear RAM by deleting large variables and running garbage collection
        del stats_noisy_XX
        del task_noisy_XX
        gc.collect()

        #--------------MERGE------------------

        #######################
        # Calc d * 2 thresholds
        #######################

        task_noisy_XX = [sinter.Task(
        circuit = surgery_circuit(
            distance = d,
            round_merge = d * 2,
            target_state_init ="X+", 
            control_state_init ="X+",
            flow_observable = "X -> XX",
            noise_after_clifford_depol = noise,
            noise_measure_flip = noise,
            noise_after_reset = noise,
            noise_depol_data_init = noise
            ),
        json_metadata={'p': noise, 'distance' : d, 'rounds' : d * 2}
        ) 
        for noise in [i for i in np.arange(0.004, 0.013, 0.001)]
        ]

        stats_noisy_XX : list[sinter.TaskStats] = sinter.collect(
            num_workers = os.cpu_count(),
            tasks=task_noisy_XX,
            decoders=['pymatching'],
            max_shots=1_000_000,
            max_errors=10_000,
            print_progress=True
        )

        #Saving Stats
        fname = f"lscx_doublerounds_merge_d{d}.pkl"
        try:
            joblib.dump(stats_noisy_XX, fname, compress=3)
            print(f"Successfully saved {fname}")
        except Exception as e:
            print(f"Failed to save {fname}: {e}")

        # Clear RAM by deleting large variables and running garbage collection
        del stats_noisy_XX
        del task_noisy_XX
        gc.collect()

        #########################
        # Calc d * 2/3 thresholds
        #########################

        task_noisy_XX = [sinter.Task(
        circuit = surgery_circuit(
            distance = d,
            round_merge = int(d * 2/3),
            target_state_init ="X+", 
            control_state_init ="X+",
            flow_observable = "X -> XX",
            noise_after_clifford_depol = noise,
            noise_measure_flip = noise,
            noise_after_reset = noise,
            noise_depol_data_init = noise
            ),
        json_metadata={'p': noise, 'distance' : d, 'rounds' : int(d * 2/3)}
        ) 
        for noise in [i for i in np.arange(0.004, 0.013, 0.001)]
        ]

        stats_noisy_XX : list[sinter.TaskStats] = sinter.collect(
            num_workers = os.cpu_count(),
            tasks=task_noisy_XX,
            decoders=['pymatching'],
            max_shots=1_000_000,
            max_errors=10_000,
            print_progress=True
        )

        #Saving Stats
        fname = f"lscx_twothirdrounds_merge_d{d}.pkl"
        try:
            joblib.dump(stats_noisy_XX, fname, compress=3)
            print(f"Successfully saved {fname}")
        except Exception as e:
            print(f"Failed to save {fname}: {e}")

        # Clear RAM by deleting large variables and running garbage collection
        del stats_noisy_XX
        del task_noisy_XX
        gc.collect()

    ###########################
    # Adding all files together
    ###########################

    def combine_files(round_type):

        pattern = f"lscx_{round_type}_d*.pkl"
        files = sorted(glob.glob(pattern))
        combined_stats = []
        for f in files:
            part = joblib.load(f)
            combined_stats.extend(part)
        out_name = f"lscx_{round_type}_combined.pkl"
        joblib.dump(combined_stats, out_name, compress=3)
        print(f"Combined and saved as {out_name}")

    combine_files("doublerounds_merge")
    combine_files("doublerounds_split")
    combine_files("twothirdrounds_merge")
    combine_files("twothirdrounds_split")

    ######################################
    # Delete all individual distance files
    ######################################

    for fname in glob.glob("lscx_*_d*.pkl"):
        try:
            os.remove(fname)
            print(f"Deleted {fname}")
        except Exception as e:
            print(f"Failed to delete {fname}: {e}")