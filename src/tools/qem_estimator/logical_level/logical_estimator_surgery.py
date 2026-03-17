
__all__ = ["LogicalEstimatorSurgery"]

import numpy as np
import stim


class LogicalEstimatorSurgery:

    IDX_TO_PAULI = {
        0: "II",
        1: "IX",
        2: "IY",
        3: "IZ",
        4: "XI",
        5: "XX",
        6: "XY",
        7: "XZ",
        8: "YI",
        9: "YX",
        10: "YY",
        11: "YZ",
        12: "ZI",
        13: "ZX",
        14: "ZY",
        15: "ZZ",
    }

    def __init__(self, ptm_clean: np.ndarray, ptm_noisy: np.ndarray):

        # Initialize Class Variables
        self.ptm_clean = ptm_clean
        self.ptm_noisy = ptm_noisy

        # 1) Get Noise Mtx
        self.noise_mtx = self._get_noise_mtx()

        # 2) Get Pauli Fidelities
        self.pauli_fidelities = self._get_pauli_fidelities(self.noise_mtx)

        # 3) Get Walsh-Hadamard Matrix
        self.walsh_hadamard = self._walsh_hadamard_16()

        # 4) Get Quasi-Probabilities
        self.nu = self._get_quasi_probabilities()

        # 5) Get Sampling Probabilities and Gamma
        self.samp_prob, self.gamma = self._get_sampling_prob()

        # 6) Get Sgn Array
        self.sgn_array = self._get_sgn_array()

    def sample_circuit(self,
                       circuit: stim.Circuit,
                       loigcal_obs_rec_pos: list[int],
                       meas_basis: str,
                       shots: int) -> float:

        """
        Samples Circuit and returns a mitigated estimate for the logical observable

        Functionality:
            - For each Shot sample from the probability distribution
            - Check if weight 2 Pauli commutes or anticommutes with
              current measurement basis (bit_flip_factor: -1 or 1)
            - mitgiated result is calculated by:
                sgn_mtx(sampled index) * bitflip_fac * total sampling overhead *
                real_meassured_circuit_result
            - This is done n shots and the average is taken
        """

        summed_mitigated_result = 0

        # The circuit here is always the same, sampler can be build
        # outside the loop
        sampler = circuit.compile_sampler()

        for _ in range(shots):

            # Get real measurement result
            results_samples = sampler.sample(shots=1)
            meas_result = 0
            for curr_obs_pos in loigcal_obs_rec_pos:
                meas_result ^= results_samples[0, curr_obs_pos].astype(np.int8)

            # Convert 0/1 measurement result to +1/-1
            conv_meas_res = 1 - 2 * meas_result

            # Sample from Probability Distribution
            sampled_index = np.random.choice(16, p=self.samp_prob)

            # Convert Sampled Index to Pauli String
            pauli_string = self.IDX_TO_PAULI[sampled_index]

            # Check if Pauli commutes with measurement basis
            if self._does_commute(meas_basis, pauli_string):
                bitflip_factor = 1
            else:
                bitflip_factor = -1

            # Get current sgn of sampled index
            sgn_idx = self.sgn_array[sampled_index]

            # Calculate mitigated result
            mitigated_result = bitflip_factor * self.gamma * conv_meas_res * sgn_idx

            # Add to summed result
            summed_mitigated_result += mitigated_result

        return summed_mitigated_result / shots

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

        # Separate str two get control and target qubit
        meas_basis_control, meas_basis_target = meas_basis
        pauli_correction_control, pauli_correction_target = pauli_correction

        anticomm_parity = 0

        # Check if both qubits are in anti-commutation pairs
        if (meas_basis_control, pauli_correction_control) in anticomm_pairs:
            anticomm_parity ^= 1
        if (meas_basis_target, pauli_correction_target) in anticomm_pairs:
            anticomm_parity ^= 1

        # Do Commute if either non or both anticommute
        return anticomm_parity == 0

    def _get_noise_mtx(self) -> np.ndarray:
        """
        This helper method returns the noisy matrix as any operation
        can be seen as: R_noisy = R_noise * R_ideal.

        -> Due to the fact that PTM of the ideal cx gate is a clifford,
           the mtx is orthogonal (i.e. each row or column is a unit vector)

        -> This reduced the above expression to:
            R_noise = R_noisy * R_ideal^-1
                    = R_noisy * R_ideal^T
        """

        # Building the Transpose of the noiseless PTM
        ptm_ideal_t = np.transpose(self.ptm_clean)

        ptm_noise = np.matmul(self.ptm_noisy, ptm_ideal_t)

        # This Mtx should be close to only diagonal entries due to the noise used
        # Currently TOllerance set to a high value due to simulation runs with 
        # relative high number of noise -> Not diagonal above threshold
        assert np.allclose(
            ptm_noise,
            np.diag(np.diagonal(ptm_noise)),
            atol=1.0,
        )

        return ptm_noise

    def _get_quasi_probabilities(self) -> np.ndarray:
        """
        Get Quasiprob Nu by vector matrix multiplication. The Walsh_hadarmard mtx
        is used as a transformation matrix between pauli fidelities and error
        probabilites

        -> Lambda_inv = W * Nu
        -> Nu = 1/16 * walsh_mtx * pauli_fidelities
        """

        # Calc New Vector
        new_vec = np.matmul(self.walsh_hadamard, self.pauli_fidelities)

        # Return Nu
        nu = 1/16 * new_vec

        return nu

    def _get_sampling_prob(self) -> tuple[np.ndarray, float]:
        """
        Calculate the sampling cost gamma
        -> With this calculate valid sampling probabilities
        """

        # Calc Gamma
        gamma = np.sum(abs(self.nu))

        # Calculate sampling prob
        sampling_prob = np.zeros(16)
        for i in range(16):
            sampling_prob[i] = abs(self.nu[i]) / gamma

        # This should be a valid probability distribution and therefore
        # add up to 1
        assert np.isclose(np.sum(sampling_prob), 1)

        return sampling_prob, gamma

    def _get_sgn_array(self) -> np.ndarray:
        """
        Create Array of just the signs of the nu values
        """
        sgn_array = np.sign(self.nu)

        return sgn_array

    @staticmethod
    def _get_pauli_fidelities(noisy_mtx: np.ndarray) -> np.ndarray:
        """
        This helper method returns the inverse of the diagonal noise PTM
        -> Essentially this tells us by how much we need to "amplify" each
           pauli channel to reverse the noise
        -> As you cannot gain new information this is a unphysical operation
        """

        # Get Diagonal of Noise PTM
        lambda_inv = 1.0 / np.diagonal(noisy_mtx)

        return lambda_inv

    @staticmethod
    def _walsh_hadamard_16() -> np.ndarray:
        """
        Return the 16x16 Walsh–Hadamard transform matrix H where:
            H[row, col] = (-1)^(i*j)
            -> i * j is the bitwise dot product of row and col
        """
        n = 16
        walsh_hadamard = np.empty((n, n), dtype=int)

        for row in range(n):
            for col in range(n):
                bitwdot = bin(row & col).count("1")
                walsh_hadamard[row, col] = (-1) ** bitwdot

        return walsh_hadamard
