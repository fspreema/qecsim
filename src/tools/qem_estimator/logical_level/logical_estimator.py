from itertools import product

import cvxpy as cp
import numpy as np
import stim

# Define the 1-qubit Pauli basis
I = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

PAULIS = [I, X, Y, Z]
PAULI_LABELS = ["I", "X", "Y", "Z"]

class GeneralLogicalEstimator:
    def __init__(self, ptm_noisy: np.ndarray, ptm_ideal: np.ndarray):
        
        # Initialize Class Variables
        self.ptm_noisy = ptm_noisy
        self.ptm_ideal = ptm_ideal

        # Get Mtx Dimensions
        self.n_rows, _ = self.ptm_noisy.shape
        self.num_qubits = int(np.log2(self.n_rows) / 2)

    def setup_estimator(self):
        """
        This method sets up the estimator by performing the following steps:
        """

        # 1) Get Noise Vector
        self.noise_vec = self._stack_vectors_of_mtrx(self._get_inv_noise_mtx())

        # 2) Get n-qubit Pauli PTMs and labels for each entry
        ptm_computation_basis_mtx, self.pauli_labels = self._get_ptm_computation_basis()

        # 3) Transfrom the basis PTMs to vectors -> Add these as columns to a matrix A 
        #   for the optimization problem (A * x = b)
        self.vec_pauli_mtx = np.zeros((self.noise_vec.shape[0], len(ptm_computation_basis_mtx)))
        for curr_index in range(len(ptm_computation_basis_mtx)):
            self.vec_pauli_mtx[:, curr_index] = self._stack_vectors_of_mtrx(ptm_computation_basis_mtx[curr_index])

        # 4) Solve Optimization Problem to get Quasi-Probabilities
        self.quasi_probabilities = self.solve_optimization()

        # 5) Compute Sign Array
        self.sgn_array = np.sign(self.quasi_probabilities)

        # 6) Get Sampling Probabilities and Gamma
        self.sampling_prob, self.gamma = self.get_sampling_probabilities_and_gamma()

    def sample_circuit(self,
                       circuit: stim.Circuit,
                       logical_obs_rec_pos: list[int],
                       meas_basis: str,
                       shots: int) -> float:
        """
        Perform Monte Carlo Sampling to get a mitigated estimate for the logical observable.
        Steps for each shot:

            - Sample Noisy Circuit and get the raw logical observable estimate
            - Sample from the Quasi-Probability distribution to get a correction term
            - Check if Pauli commutes or anticommutes with
              current measurement basis (bit_flip_factor: -1 or 1)
            - mitgiated result is calculated by:
                sgn_mtx(sampled index) * bitflip_fac * total sampling overhead *
                real_meassured_circuit_result
            - This is done n shots and the average is taken
        """

        summed_mitigated_result = 0
        number_of_operations = len(self.sampling_prob)

        # Build Sampler and sample the circuit to get the real measurement results for each shot
        sampler = circuit.compile_sampler()
        results_samples = sampler.sample(shots=shots)

        # Sample from the probability dist shots times
        sampled_index = np.random.choice(number_of_operations, size=shots, p=self.sampling_prob)


        for shot_idx in range(shots):

            # Get real measurement result for curr shot
            meas_result = 0
            for curr_obs_pos in logical_obs_rec_pos:
                meas_result ^= results_samples[shot_idx, curr_obs_pos].astype(np.int8)

            # Convert 0/1 measurement result to +1/-1
            conv_meas_res = 1 - 2 * meas_result

            # Convert Sampled Index to Pauli String
            pauli_string = self.pauli_labels[sampled_index[shot_idx]]

            # Check if Pauli commutes with measurement basis
            if self._does_commute(meas_basis, pauli_string):
                bitflip_factor = 1
            else:
                bitflip_factor = -1

            # Get current sgn of sampled index
            sgn_idx = self.sgn_array[sampled_index[shot_idx]]

            # Calculate mitigated result
            mitigated_result = bitflip_factor * self.gamma * conv_meas_res * sgn_idx

            # Add to summed result
            summed_mitigated_result += mitigated_result

        return summed_mitigated_result / shots

    def _unitary_to_ptm(self, unitary: np.ndarray) -> np.ndarray:
        """
        This helper method takes a 1 qubit unitary matrix and converts it to its PTM representation.

        Given a unitary U, the PTM representation is given by:
        R_{ij} = (1/d) * Tr(P_i * U * P_j * U^\dagger)
        where P_i and P_j are the Pauli operators and d is the dimension of the Hilbert space (d = 2^n for n qubits)
        """

        # Initializing the PTM matrix
        ptm = np.zeros((4, 4), dtype=float)

        # Compute the PTM elements using the formula above
        for i, p_i in enumerate(PAULIS):
            for j, p_j in enumerate(PAULIS):
                ptm[i, j] = np.real(0.5 * np.trace(p_i @ unitary @ p_j @ unitary.conj().T))

        return ptm

    def _get_ptm_computation_basis(self) -> tuple[list[np.ndarray], list[str]]:
        """
        This helper method returns the basis gates for the given measurement basis.

        -> For 1 qubit: I, X, Y, Z & Clifford gates
        -> For 2 qubit: tensor products of the above gates
        -> For n qubits: tensor products of the above gates
        """

        single_qubit_ptm = [self._unitary_to_ptm(gate) for gate in PAULIS]

        if self.num_qubits == 1:
            # Return single qubit gates as the basis for the computation
            # and the corresponding label ordering
            return single_qubit_ptm, PAULI_LABELS

        elif self.num_qubits > 1:
            multi_qubit_ptm = []
            labels = []
            
            # Get all combinations with itertools.product
            for ptm_combination in product(single_qubit_ptm, repeat=int(self.num_qubits)):
                # For each combination, compute the tensor product and append to the multi_qubit_ptm list
                combined_ptm = None
                for curr_ptm in ptm_combination:
                    if combined_ptm is None:
                        combined_ptm = curr_ptm
                    else:
                        combined_ptm = np.kron(combined_ptm, curr_ptm)
                multi_qubit_ptm.append(combined_ptm)
            
            # Generate corresponding labels
            for label_combination in product(PAULI_LABELS, repeat=int(self.num_qubits)):
                labels.append(''.join(label_combination))

            return multi_qubit_ptm, labels

        else:
            raise ValueError("Number of qubits must be greater than 0")

    
    def _get_inv_noise_mtx(self) -> np.ndarray:
        """
        This helper method returns the inverse of the noisy matrix as any operation
        can be seen as: R_noisy = R_noise * R_ideal.

        -> This reduced the above expression to:
            R_noise = R_noisy * R_ideal^-1
        """

        # Building the Transpose of the noiseless PTM
        ptm_ideal_inv = np.linalg.inv(self.ptm_ideal)
        ptm_noise = np.matmul(self.ptm_noisy, ptm_ideal_inv)

        # We force the noise mtx to be diagonal as non-pauli noise is not supported
        ptm_noise = np.diag(np.diag(ptm_noise))

        # Invert noise mtx to get pauli fidelities
        ptm_noise_inv = np.linalg.inv(ptm_noise)

        return ptm_noise_inv

    def solve_optimization(self) -> np.ndarray:
        """
        Solves the optimization problem to find the quasi-probabilities (find minimal L1 norm solution)

        A * x = b

        where A is the matrix of Pauli fidelities, 
        x is the vector of quasi-probabilities, and 
        b is the vector of the noise ptm
        """

        # Define the optimization variable
        quasi_prob = cp.Variable(self.vec_pauli_mtx.shape[1])

        # Define the objective function (L1 norm)
        objective = cp.Minimize(cp.norm(quasi_prob, 1))

        # Define the constraints (A * x = b)
        constraints = [self.vec_pauli_mtx @ quasi_prob == self.noise_vec]

        # Formulate and solve the optimization problem
        prob = cp.Problem(objective, constraints)
        prob.solve()

        # Check if a solution was found
        if prob.status not in ["optimal", "optimal_inaccurate"]:
            raise RuntimeError(f"Optimization failed to find an optimal solution. Status: {prob.status}")

        return quasi_prob.value
    
    def get_sampling_probabilities_and_gamma(self) -> tuple[np.ndarray, float]:
        """
        This method computes the sampling probabilities and gamma from the quasi-probabilities.
        """

        # Compute gamma as the sum of the absolute values of the quasi-probabilities
        gamma = np.sum(np.abs(self.quasi_probabilities))

        # Compute sampling probabilities by normalizing the absolute values of the quasi-probabilities
        sampling_probabilities = np.abs(self.quasi_probabilities) / gamma

        return sampling_probabilities, gamma

    def _stack_vectors_of_mtrx(self, mtx: np.ndarray) -> np.ndarray:
        """
        Function which takes a given mtx and stacks its columns on 
        top of each other to create a vector
        """

        # What shape is needed?
        shape = self.n_rows ** 2

        reshaped_mtx = mtx.reshape(shape, order="F")

        return reshaped_mtx
    
    @staticmethod
    def _does_commute(meas_basis:str, pauli_correction: str) -> bool:
        """
        Helper which checks if the measurement basis f.ex IX commutes
        with the current corrective pauli correction f.ex. ZX
        """

        # Single-qubit anti-commutation pairs
        anticomm_pairs = {
            ("X", "Y"), ("Y", "X"),
            ("X", "Z"), ("Z", "X"),
            ("Y", "Z"), ("Z", "Y"),
        }

        anticomm_parity = 0

        # Deconstruct the measurement basis and pauli correction into their single-qubit components
        for meas_op, pauli_op in zip(meas_basis, pauli_correction, strict=True):
            if (meas_op, pauli_op) in anticomm_pairs:
                anticomm_parity += 1

        # Do Commute if equal number of anti-commutations (including 0)
        return anticomm_parity % 2== 0
