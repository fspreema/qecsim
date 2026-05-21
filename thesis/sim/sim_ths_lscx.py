import itertools
import math
import os
import sinter
from src.codes.lattice_surgery.builder import SurgeryBuilder
from src.core.data_models import NoiseParameters
from tqdm import tqdm

def main():
    # Get the directory where script is saved
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

    # Go up two levels to reach root
    WORKSPACE_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "../../")) 
    
    # Get thesis/data directory
    DATA_DIR = os.path.join(WORKSPACE_ROOT, "thesis", "data")
    os.makedirs(DATA_DIR, exist_ok=True)

    # Intial Setup Values
    geom = [10**(math.log10(1e-5) + i*(math.log10(1e-3) - math.log10(1e-5))/9) for i in range(10)]
    lin = [1.1e-3 + i*(0.015 - 1.1e-3)/39 for i in range(40)]
    p_values = geom + lin
    distances = [3, 5, 7, 9]
    combinations = list(itertools.product(distances, p_values))

    # Combinations to simulate
    cases = [
        ("I0", "Z0", "Z", "Z"), 
        ("Z0", "I0", "Z", "I"), 
        ("X+", "I0", "X", "X"), 
        ("I0", "X+", "I", "X")
    ]

    for curr_config in cases:
        # Build the tasks -> Usinf tqdm due to larger duration of circuit compilation
        tasks = [
            sinter.Task(
                circuit=SurgeryBuilder(
                    distance=d,
                    control_state_init=curr_config[0],
                    target_state_init=curr_config[1],
                    control_measure_basis=curr_config[2],
                    target_measure_basis=curr_config[3],
                    noise=NoiseParameters(
                        circuit_noise_prob=noise,
                    )
                ).build_circuit(),
                json_metadata={'d': d, 'p': noise, 'config': curr_config},
            )
            for d, noise in tqdm(combinations, desc=f"Building {curr_config}")
        ]

        # Set Filename
        filename = os.path.join(DATA_DIR, f"stats_surgery_{'_'.join(curr_config)}.csv")

        # Hand them off to sinter
        collected_stats = sinter.collect(
            num_workers=os.cpu_count(),
            tasks=tasks,
            decoders=['pymatching'],
            max_shots=10_000_000,
            max_errors=150_000,
            print_progress=True,
            save_resume_filepath= filename
        )
        print(f"Finished processing configuration: {curr_config}\n")

if __name__ == "__main__":
    main()