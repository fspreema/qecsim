import sinter
import os, pickle, itertools, collections
import functools
import random
import numpy as np

from collections import namedtuple, defaultdict
from multiprocessing import Manager, Pool, cpu_count

from qecsim.xzzx_code.circuit import XZZX_code
from qecsim.threshold_num.calc_threshold import threshold_approx

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

    def compiled_xzzx(dist: int, bias_key: tuple, noise: float):
        key = (dist, bias_key, noise)
        if key not in circuit_db:   # compile once per unique key
            circuit_db[key] = XZZX_code(
                distance=dist,
                rounds=dist,
                state_init="Ver",
                noise_bias=list(bias_key),
                after_c_pauli_channel_prob=noise,
            )
        return circuit_db[key]

    ################
    # Defining Tasks
    ################

    def circuit_factory(d, b, p):
        return compiled_xzzx(d, b, p)

    def make_tasks_for_bias(bias):
        return [
            sinter.Task(
                circuit = compiled_xzzx(dist = d, bias_key = bias, noise = p),
                json_metadata={'p': p, 'distance': d, 'bias': bias},
            )
            for d in (9, 11)
            for p in np.arange(0.005, 0.1, 0.005)
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

    ################
    # Sinter Collect
    ################

    stats_all_bias = sinter.collect(
        tasks=tasks,
        decoders=['pymatching'],
        num_workers = os.cpu_count(),
        max_shots = 50_000,
        max_errors = 5_000,
        max_batch_size = 5,
        print_progress=True,
    )

    """
    As we have a full task list we have to filter out all the individual stats for the json_metadata with the correct bias
    """

    stats_by_bias = collections.defaultdict(list)
    
    for elements in stats_all_bias:
        stats_by_bias[tuple(elements.json_metadata["bias"])].append(elements)


    results : list = []

    #####################
    # Calculate threshold
    #####################

    for bias, sub_stats in stats_by_bias.items():

        try:
            calc_th = threshold_approx(sub_stats)
            results.append([bias,calc_th])

        except Exception:
            results.append([bias,-1])

    ##############################
    #Saving num_values with pickle
    ##############################

    with open("XZZX_num_value(full_bias).pkl", "wb") as file:
        pickle.dump(results, file)