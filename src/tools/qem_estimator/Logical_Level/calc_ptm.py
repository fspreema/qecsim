import numpy as np
import pymatching
import stim

from src.core.data_models import PTMCircuits

__all__ = ["PTMCalculator"]


class PTMCalculator:
    def __init__(
        self,
        circuits: PTMCircuits,
        samples: int = 1_000,
    ):
        """
        Building PTM by building each possible memory circuit and measurement combination

        What we do:
            * Determine the real noisy measurement value of the observable
            * Construct the DEM and Matcher for the current circuit
            * Run QEC and check whether the noisy observable should be flipped or not
            * Flip the observable i.e. the noisy measurement value
            * Build up PTM out of the expectation values

        Arguments:
            *circuits: Data model containing all circuits required for PTM calculation
            *samples: Number of samples to take for each circuit
            *only_diagonal: Only the diagonal Values are calculated

        """

        # Set perliminary attributes
        self.samples = samples
        self.circuits = circuits.circuits

    def calc_ptm(self):
        # Calculate the expectation values for each basis combination
        exp_vals_per_basis_comb = self._calc_entries()

        # Build the PTM Matrix out of the expectation values
        ptm_matrix = self._build_mtx_for_ptm(exp_vals_per_basis_comb)

        return ptm_matrix

    def _build_mtx_for_ptm(self, exp_vals_per_basis_comb: dict[str, float]) -> np.ndarray:
        """
        Converts the dictionary of expentation values into the PTM Matrix

        Functionality:
            1) Determine size of Matrix by number of entires
            2) Map Pauli Indices into Numbers
            3) Map Pauli strings from dict to matrix entires
            4) Populate Matrix by [input_index, output_index] = exp_val
        """

        # Determine the size of the PTM Matrix
        first_key = next(iter(exp_vals_per_basis_comb))
        input_pauli_str = first_key.split("->")[0]
        n_qubits = len(input_pauli_str)

        if n_qubits == 1:
            len_mtx = 3
        elif n_qubits == 2:
            len_mtx = 15

        # Initialize empty PTM Matrix
        ptm_matrix = np.zeros((len_mtx, len_mtx))

        # Define mapping from Pauli String to Matrix Indices
        for pauli_comb, exp_val in exp_vals_per_basis_comb.items():
            input_pauli, output_pauli = pauli_comb.split("->")
            input_index = self._map_pauli_string_to_indices(input_pauli)
            output_index = self._map_pauli_string_to_indices(output_pauli)

            # Populate the PTM Matrix
            ptm_matrix[input_index, output_index] = exp_val

        return ptm_matrix

    @staticmethod
    def _map_pauli_string_to_indices(pauli_str: str) -> int:
        """
        Given a already converted pauli number (e.g. 0 for II, 1 for IX, etc.) this function returns
        the corresponding input and output pauli indices
        """

        # Single Pauli Case
        if len(pauli_str) == 1:
            pauli_to_num = {"X": 0, "Y": 1, "Z": 2}
            return pauli_to_num[pauli_str]

        elif len(pauli_str) == 2:
            pauli_to_num = {"I": 0, "X": 1, "Y": 2, "Z": 3}
            pauli_1, pauli_2 = pauli_str

            return_index = pauli_to_num[pauli_1] * 4 + pauli_to_num[pauli_2]
            return return_index

        else:
            raise ValueError(f"Invalid Pauli String: {pauli_str}")

    @staticmethod
    def _get_logical_meas_from_samples(
        samples_from_sampler: np.ndarray,
        rec_pos: list[int],
    ) -> np.ndarray:
        """
        Returns full logical measurement by xoring the needed measurements
        -> This is Vectorized to speed up the process as much as possible
        """

        final_meas: np.ndarray = np.empty(shape=(samples_from_sampler.shape[0],), dtype=np.bool_)

        for curr_rec in rec_pos:
            # XOR current measurement with the final measurement (Vectorized)
            final_meas ^= samples_from_sampler[:, curr_rec]

        return final_meas

    @staticmethod
    def _xor_meas_and_decoder(
        decoder_prediction: np.ndarray,
        logical_state_meas: np.ndarray,
    ) -> np.ndarray:
        """
        Returns the final measurement for all shots (vectorized).
        -> XOR logical measurement with the decoder prediction
        """

        # Takin all decoder Predictions from observable 0
        decoder_prediction_col = decoder_prediction[:, 0]

        # XOR Measurement with Decoder Prediction (Vectorized)
        # -> Converts the boolean values into +1/-1 values for the final logical measurement
        xor_result = decoder_prediction_col.astype(np.int8) ^ logical_state_meas.astype(np.int8)
        final_log_meas = 1 - 2 * xor_result

        return final_log_meas

    @staticmethod
    def _get_dem_and_matcher(
        circuit: stim.Circuit,
    ) -> tuple[pymatching.Matching, stim.DetectorErrorModel]:
        """
        Creates for a given stim Circuit the detector error model as well as the pymachting matcher
        """

        dem = circuit.detector_error_model(decompose_errors=True)
        matcher = pymatching.Matching.from_detector_error_model(dem)

        return matcher, dem

    def _calc_entries(self):
        """
        We can't directly use the DEM as we need the
        raw measurements to infer what logical state we have

        -> We use compile sampler to infer the logical state
        -> Convert into a dem to run the matching
        """

        # Init Save Dict for the expectation values of each basis combination
        exp_vals_per_basis_comb: dict[str, float] = {}

        # Run through diagonal circuits
        for basis_combination, circuit_and_meas_rec in self.circuits.items():
            circuit, meas_rec = circuit_and_meas_rec

            # Buid the normal measurement smaples and sample n shots
            sampler = circuit.compile_sampler()
            results_samples = sampler.sample(shots=self.samples)

            # Getting the individual logical states after measurement (Vectorized)
            logical_states_noisy = self._get_logical_meas_from_samples(
                results_samples,
                meas_rec,
            )

            # Build DEM and Matcher
            matcher, dem = self._get_dem_and_matcher(circuit)

            # Converting the measurement sample into DEM sample and continue as usual with decoding
            detectors_observable_states = circuit.compile_m2d_converter().convert(
                measurements=results_samples,
                append_observables=True,
            )
            detectors = detectors_observable_states[:, : dem.num_detectors]

            # Decode all sample round in one Batch decode
            is_flipped_prediction = matcher.decode_batch(detectors)

            # Calculate final logical states (Vectorized)
            # This returns an array of +1/-1 values
            final_logical_states = self._xor_meas_and_decoder(
                is_flipped_prediction,
                logical_states_noisy,
            )

            # Take the Average of the logical state over all samples
            average_logical_state = np.average(final_logical_states)

            # Save the average logical state for the current basis combination to build the PTM
            exp_vals_per_basis_comb[basis_combination] = average_logical_state

        return exp_vals_per_basis_comb
