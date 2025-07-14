import sinter
import os
from qecsim.xzzx_code.circuit import XZZX_code
from qecsim.threshold_num.calc_threshold import threshold_approx
import pickle
import numpy as np

__all__ = ["sim_run"]

#-----------------------------Calc----------------------------

# Initlize Saving Array
num_value : list = []

# Changing Bias setting
for b_x in np.arange(0.0, 1.2, 0.2): 
    for b_y in np.arange(0.0, 1.2, 0.2): 
        for b_z in np.arange(0.0, 1.2, 0.2):

            #Creating Current Bias
            current_bias = [0,b_y,b_z]

            if current_bias == [0,0,0]:
                continue

            # Running simulation
            task_noisy_v = [sinter.Task(
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
                ]

            stats_noisy_v : list[sinter.TaskStats] = sinter.collect(
                num_workers = os.cpu_count(),
                tasks=task_noisy_v,
                decoders=['pymatching'],
                max_shots=100_000_0,
                max_errors=10_000,
                print_progress=True
            )

            # calculate threshold and save -> If Error skip and set 0 
            try:
                calc_th = threshold_approx(stats_noisy_v)
                num_value.append([current_bias,calc_th])

            except Exception:
                num_value.append([current_bias,0])

#Saving num_values with pickle
with open("XZZX_num_value(full_bias).pkl", "wb") as file:
    pickle.dump(num_value, file)