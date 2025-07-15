import sinter
import os, pickle, itertools, collections
from qecsim.xzzx_code.circuit import XZZX_code
from qecsim.threshold_num.calc_threshold import threshold_approx
import numpy as np

if __name__ == "__main__":
    
    ###################################
    # Define list of valid bias options
    ###################################

    step = 0.0125
    bias_steps = np.arange(0.0, 1 + (step/2), step)

    all_bias_triplets = [[bx, by, 1 - bx - by] 
                         for bx, by in itertools.product(bias_steps, repeat = 2) 
                         if 0 <= 1 - bx - by <= 1
                         and not (bx == 0 and by == 0)]

    ##############################
    # Running full task simulation
    ##############################

    task_all_bias = [sinter.Task(
        circuit = XZZX_code(
            distance = d,
            rounds = d, 
            state_init ="Ver",
            noise_bias = current_bias,
            after_c_pauli_channel_prob = noise
            ),
        json_metadata={'p': noise, 'distance' : d, 'bias' : current_bias}
        ) 
        for noise in [i for i in np.arange(0.005, 0.1, 0.005)]
        for d in [9, 11]
        for current_bias in all_bias_triplets
        ]

    stats_all_bias : list[sinter.TaskStats] = sinter.collect(
        num_workers = os.cpu_count(),
        tasks=task_all_bias,
        decoders=['pymatching'],
        max_shots=1_000_000,
        max_errors=5_000,
        print_progress=True
    )

    #########################################################
    # calculate threshold and save -> If Error skip and set 0
    #########################################################

    """
    As we have a full task list we have to filter out all the individual stats for the json_metadata with the correct bias
    """

    stats_by_bias = collections.defaultdict(list)
    
    for elements in stats_all_bias:
        stats_by_bias[tuple(elements.json_metadata["bias"])].append(elements)


    results : list = []

    for bias, sub_stats in stats_by_bias.items():

        try:
            calc_th = threshold_approx(sub_stats)
            results.append([bias,calc_th])

        except Exception:
            results.append([bias,-1])

    #Saving num_values with pickle
    with open("XZZX_num_value(full_bias).pkl", "wb") as file:
        pickle.dump(results, file)