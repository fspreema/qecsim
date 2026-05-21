import functools
from concurrent.futures import ProcessPoolExecutor
from itertools import product

import numpy as np
import pandas as pd
import stim
import os
from tqdm import tqdm

from src.codes.surface_code_rotated.builder import SurfaceBuilder
from src.core.data_models import NoiseParameters, PTMCircuits
from src.tools.qem_estimator.logical_level.calc_ptm import PTMCalculator
from src.tools.qem_estimator.logical_level.logical_estimator import GeneralLogicalEstimator

__all__ = ["GetPTMThreshold"]

# Pauli alphabet for input and output states
PAULIS = ["I", "X", "Y", "Z"]

# Non Zeor and Diagonal Flows
DIAG_NON_ZERO = {(p, p) for p in ["X", "Y", "Z"]}

class GetPTMThreshold:

    def __init__(self, samples: int, only_diag: bool):
        """
        This file executes the simulation for the creation of the threshold diagram
        including the overhead factor gamma of the porbabilistic error cancellation technique

        Functionality:
        - For Any given Distance d in [3,5,...] we do the following:
            1) Construct the ideal noiseless PTM
            2) For any probabilitiy p in [0.01, 0.02,...] we do the following:
                1) Construct the noisy PTM
                2) Build the logical estimator using the noisy and the clean PTM
                3) Getting the overhead factor gamma
        """

        self.samples = samples
        self.only_diag = only_diag

    def run_simulation(self,
                       distances: list[int],
                       physical_err_probs: list[float], 
                       noise_type: str, 
                       bias: list[float],
                       save_ptm_files: bool = False,
                       output_folder: str = "ptm_matrices",) -> list[dict[str, float]]:

        # Init Result List
        results: list[dict[str, float]] = []

        if save_ptm_files and not os.path.exists(output_folder):
            os.makedirs(output_folder)

        for d in distances:
            print(f"\n--- Starting Distance d={d} ---")

            # Construct the ideal PTMS
            # -> These do only need to be calculated once per distance!
            ideal_ptm = self._build_ptm(distance=d, 
                                        physical_err_probs=0.0, 
                                        noise_type = noise_type, 
                                        bias = bias,
                                        )

            
            if save_ptm_files:
                ideal_filename = f"ptm_ideal_d{d}_{noise_type}.npy"
                np.save(os.path.join(output_folder, ideal_filename), ideal_ptm)

            # Prepare tasks for workers for all physical error probabilites at this fixed distance
            # Use partial function to freeze everything except the physical error prob.
            worker_task = functools.partial(
                self._build_ptm, 
                distance=d, 
                noise_type=noise_type, 
                bias=bias, 
            )

            # Implementing parralel execution
            with ProcessPoolExecutor() as executor:
                # We use tqdm to show the progress bar
                noisy_ptms = list(tqdm(
                    executor.map(worker_task, physical_err_probs),
                    total=len(physical_err_probs),
                    desc=f"Calculating Noisy PTMS for d={d}",
                ))

            # Calculate the overhead factor gamma
            for p, noisy_ptm in zip(physical_err_probs, noisy_ptms):
                gamma = self._get_overhead_gamma(ptm_ideal=ideal_ptm, ptm_noisy=noisy_ptm)
                results.append({"distance": d, "physical_error_probability": p, "gamma": gamma})

                if save_ptm_files:
                    noisy_filename = f"ptm_noisy_d{d}_p{p:.2e}_{noise_type}.npy"
                    np.save(os.path.join(output_folder, noisy_filename), noisy_ptm)

        return results

    def _build_ptm(self, 
                   physical_err_probs: float,
                   distance: int,
                   noise_type: str, 
                   bias: list[float],
                   ) -> np.ndarray:

        # Input checks
        assert noise_type in ["CircuitNoise", "BiasNoise"], "Invalid noise type. Supported types are 'CircuitNoise' and 'BiasNoise'."

        #########################
        # Construct Noise Class #
        #########################

        if noise_type == "CircuitNoise":

            noise_class = NoiseParameters(circuit_noise_prob= physical_err_probs)

        else:
            raise ValueError(f"Invalid noise type: {noise_type}. Supported types are 'CircuitNoise' and 'BiasNoise'.")

        ############
        # Circuits #
        ############

        # Create Mapping for input states
        input_to_init_state = {
            "X": ["X+", "X-"],
            "Y": ["Y+", "Y-"],
            "Z": ["Z0", "Z1"],
            "I": ["I0", "I1"],
        }

        circuits_memory: dict[str, tuple[dict[str, stim.Circuit], list[int]]] = {}

        if self.only_diag:
            total_combinations = DIAG_NON_ZERO
        else:
            # We combine the products so the progress bar tracks all 16 combinations
            total_combinations = list(product(PAULIS, PAULIS))

        for p_out, p_in in total_combinations:

            # Creating Label
            label_basis = f"{p_in}->{p_out}"

            # As all ...->I flows have no observable, sampling is not possible
            if p_out == "I":
                continue

            # Creating inner dict
            inner_dict: dict[str, stim.Circuit] = {}

            # Inner loop for the 4 initial states
            for init_state in input_to_init_state[p_in]:

                builder = SurfaceBuilder(distance=distance, 
                                         state_init= init_state, 
                                         log_obs= p_out, 
                                         noise=noise_class)

                circuit = builder.build_circuit()

                # It doesnt matter which records we get as the logical observable stays the same
                # in this loop
                curr_meas_rec = builder.get_logical_meas_rec(observable_index=0)

                # Adding the circuit to the inner dict
                inner_dict[init_state] = circuit

            # Adding to your existing dict
            circuits_memory[label_basis] = (inner_dict, curr_meas_rec)

        ############################
        # Calculate the PTM-Matrix #
        ############################

        ptm_calculator = PTMCalculator(PTMCircuits(circuits=circuits_memory), 
                                       samples= self.samples, 
                                       pauli_channel_2_used= (noise_type == "BiasNoise"))
        # Doesnt matter as reuqirements for bot are the same
        ptm_mtx = ptm_calculator.calc_ptm(sparse_ideal_ptm= False, sparse_noisy_ptm= True)

        return ptm_mtx

    @staticmethod
    def _get_overhead_gamma(ptm_ideal: np.ndarray, ptm_noisy: np.ndarray):

        # Set II entry on both to one
        ptm_ideal[0, 0] = 1
        ptm_noisy[0, 0] = 1

        # Building Class
        cls = GeneralLogicalEstimator(ptm_ideal=ptm_ideal, ptm_noisy=ptm_noisy)
        cls.setup_estimator()

        return cls.gamma

# Run Simulation
if __name__ == "__main__":
    # Settings
    ds = [3, 5, 7, 9]
    ps = np.concatenate([
    np.geomspace(1e-5, 1e-3, 10),    # 10 points logscaling
    np.linspace(1.1e-3, 1e-2, 40),   # 40 dense linear points
    ])
    
    # Executing Simulation
    sim = GetPTMThreshold(samples=10_000, only_diag=True)
    

    # --- Experiment 1: Standard Depolarizing ---
    results_std = sim.run_simulation(
        distances=ds,
        physical_err_probs=ps,
        noise_type="CircuitNoise",
        save_ptm_files=True,
        output_folder="thesis/data/ptm_matrices_memory",
        bias=[0, 0, 0]
    )
    pd.DataFrame(results_std).to_csv("thesis/data/gamma_memory.csv")

    """
        # --- Experiment 2: Z-Biased Noise ---
        results_bias = sim.run_simulation(
            distances=ds,
            physical_err_probs=ps,
            noise_type="BiasNoise",
            bias=[0.01, 0.01, 0.98],
            samples=1000
        )
        pd.DataFrame(results_bias).to_csv("gamma_biased.csv")
    """
