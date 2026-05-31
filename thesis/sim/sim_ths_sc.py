import itertools
import numpy as np
import os
import sinter
from src.codes.surface_code_rotated.builder import SurfaceBuilder
from src.core.data_models import NoiseParameters

def main():
    # Get the directory where script is saved
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

    # Go up two levels to reach root
    WORKSPACE_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "../../")) 
    
    # Get thesis/data directory
    DATA_DIR = os.path.join(WORKSPACE_ROOT, "thesis", "data")
    os.makedirs(DATA_DIR, exist_ok=True)

    # Initial Setup Values
    p_values = np.concatenate([
    np.geomspace(1e-5, 1e-3, 10),    # 10 points logscaling
    np.linspace(1.1e-3, 1e-2, 40),   # 40 dense linear points
    ])
    distances = [3, 5, 7, 9]
    combinations = list(itertools.product(distances, p_values))

    # Combinations to simulate
    cases = [
        ("X+", "X"), 
        ("Y+", "Y"), 
        ("Z0", "Z")
    ]

    for curr_config in cases:
        print(f"--- Starting Configuration: {curr_config} ---")
        print("Building circuit tasks...")
        
        # Build the tasks
        tasks = [
            sinter.Task(
                circuit=SurfaceBuilder(
                    distance=d,
                    state_init=curr_config[0],
                    log_obs=curr_config[1],
                    noise=NoiseParameters(
                        circuit_noise_prob=noise,
                    )
                ).build_circuit(),
                json_metadata={'d': d, 'p': noise, 'config': curr_config},
            )
            for d, noise in combinations
        ]

        # Set Filename
        filename = os.path.join(DATA_DIR, f"stats_memory_{'_'.join(curr_config)}.csv")

        # Hand them off to sinter
        sinter.collect(
            num_workers=os.cpu_count(),
            tasks=tasks,
            decoders=['pymatching'],
            max_shots=10_000_000,
            max_errors=150_000,
            print_progress=True,
            save_resume_filepath=filename
        )
        print(f"Finished processing configuration: {curr_config}\n")

if __name__ == "__main__":
    main()
