import sinter
import os
import pickle
import itertools
import random
import numpy as np
import glob

from collections import defaultdict
from functools import lru_cache
from multiprocessing import Manager, Pool, cpu_count
from concurrent.futures import ProcessPoolExecutor

from qecsim.xzzx_code.circuit import XZZX_code
from qecsim.threshold_num.calc_threshold import threshold_approx

os.makedirs("results", exist_ok=True)

##########################
# Parallelizable Functions
##########################

def run_task_multi(single_task: sinter.Task):

    stats_noisy_v : list[sinter.TaskStats] = sinter.collect(
        num_workers = 1,
        tasks = single_task,
        decoders = ['pymatching'],
        max_shots = 1_000_000,
        max_errors = 10_000,
        print_progress = False
        )

    try: 
        # calculate threshold and save -> If Error skip and set 0 
        calc_th = threshold_approx(stats_noisy_v)
        result = [single_task.json_metadata["bias"],calc_th]

    except Exception:
        result = [single_task.json_metadata["bias"], -1]

    #Saving num_values with pickle
    filename = f"bias_{single_task.json_metadata['bias'][0]:.3f}_\
    {single_task.json_metadata['bias'][1]:.3f}_\
    {single_task.json_metadata['bias'][2]:.3f}.pkl"

    with open(os.path.join("results", filename), "wb") as file:
        pickle.dump(result, file)

#################################
# Run all jobs on different cpus:
#################################

if __name__ == "__main__":

    ###################################
    # Define list of valid bias options
    ###################################

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"

    step = 0.025
    bias_steps = np.arange(0.0, 1 + (step/2), step)

    all_bias_triplets = [(bx, by, 1 - bx - by) 
                        for bx, by in itertools.product(bias_steps, repeat = 2) 
                        if 0 <= 1 - bx - by <= 1
                        and not (bx == 0 and by == 0)]
        
    #######################
    # Precompiling Circuits
    #######################

    _manager   = Manager()          # must be created before workers spawn
    circuit_db = _manager.dict()    # proxy shared by all processes

    @lru_cache(maxsize=None)
    def compiled_xzzx(dist: int, bias_key: tuple, noise: float):
        key = (dist, bias_key, noise)
        if key not in circuit_db:   # compile once per unique key
            circuit_db[key] = XZZX_code(
                distance = dist,
                rounds = dist,
                state_init = "Ver",
                noise_bias = list(bias_key),
                after_c_pauli_channel_prob = noise,
            )
        return circuit_db[key]

    ################
    # Defining Tasks
    ################

    def make_tasks_for_bias(bias):
        return [
            sinter.Task(
                circuit = compiled_xzzx(dist = d, bias_key = bias, noise = p),
                json_metadata={'p': p, 'distance': d, 'bias': bias},
            )
            for d in (9, 11)
            for p in np.arange(0.005, 0.105, 0.005)
        ]

    ############################
    # Building tasks in parallel
    ############################

    with Pool(processes=min(64, cpu_count())) as pool:
        task_chunks = pool.map(make_tasks_for_bias, all_bias_triplets)

    tasks = [t for chunk in task_chunks for t in chunk]

    #######################################################
    # Shuffle Tasks -> low p distributed across all workers
    #######################################################

    random.shuffle(tasks)

    ##############################
    # Group tasks once per bias
    ##############################

    tasks_by_bias = defaultdict(list)
    for t in tasks:
        tasks_by_bias[tuple(t.json_metadata["bias"])].append(t)

    ###########################
    # Run groups in parallel
    ###########################

    with ProcessPoolExecutor(max_workers=cpu_count()) as ex:
        ex.map(run_task_multi, tasks_by_bias.values())
        
    #########################################
    # Merge the per‑bias pickles into one
    #########################################

    merged_results = []
    for fn in glob.glob("results/bias_*.pkl"):
        with open(fn, "rb") as f:
            merged_results.append(pickle.load(f))

    #sort so rows appear in a predictable order
    merged_results.sort(key=lambda row: row[0])

    #save as pickle
    with open("results/all_thresholds.pkl", "wb") as f:
        pickle.dump(merged_results, f)