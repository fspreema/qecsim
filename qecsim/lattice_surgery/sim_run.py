import pickle
import sinter
import os, pickle

from qecsim.lattice_surgery.circuit import surgery_circuit

if __name__ == "__main__":

    task_noisy_XX = [sinter.Task(
    circuit = surgery_circuit(
        distance = d,
        round_num = d * 2,
        target_state_init ="X+", 
        control_state_init ="X+",
        flow_observable = "X -> XX",
        noise_after_clifford_depol = noise,
        noise_measure_flip = noise,
        noise_after_reset = noise,
        noise_depol_data_init = noise
        ),
    json_metadata={'p': noise, 'distance' : d, 'rounds' : 3}
    ) 
    for noise in [0.001, 0.0015, 0.002, 0.0025, 0.005, 0.010, 0.015, 0.1, 0.125, 0.2, 0.25, 0.3, 0.4, 0.5]
    for d in [3, 5, 7, 9, 11, 13, 15, 17]
    ]

    stats_noisy_XX : list[sinter.TaskStats] = sinter.collect(
        num_workers = os.cpu_count(),
        tasks=task_noisy_XX,
        decoders=['pymatching'],
        max_shots=1_000_000,
        max_errors=10_000,
        print_progress=True
    )

    #Saving Stats with pickle
    with open("lscx_doublerounds.pkl", "wb") as file:
        pickle.dump(stats_noisy_XX, file)

    task_noisy_XX = [sinter.Task(
    circuit = surgery_circuit(
        distance = d,
        round_num = d * 3,
        target_state_init ="X+", 
        control_state_init ="X+",
        flow_observable = "X -> XX",
        noise_after_clifford_depol = noise,
        noise_measure_flip = noise,
        noise_after_reset = noise,
        noise_depol_data_init = noise
        ),
    json_metadata={'p': noise, 'distance' : d, 'rounds' : 3}
    ) 
    for noise in [0.001, 0.0015, 0.002, 0.0025, 0.005, 0.010, 0.015, 0.1, 0.125, 0.2, 0.25, 0.3, 0.4, 0.5]
    for d in [3, 5, 7, 9, 11, 13, 15, 17]
    ]

    stats_noisy_XX : list[sinter.TaskStats] = sinter.collect(
        num_workers = os.cpu_count(),
        tasks=task_noisy_XX,
        decoders=['pymatching'],
        max_shots=1_000_000,
        max_errors=10_000,
        print_progress=True
    )

    #Saving Stats with pickle
    with open("lscx_triplerounds.pkl", "wb") as file:
        pickle.dump(stats_noisy_XX, file)