import numpy as np
import pymatching
import stim

from src.core.data_models import PTMCircuits

__all__ = ["PTMCalculator"]


class PTMCalculator:
    BASIS_TO_PAULI = {
        "X": ["X+", "X-"],
        "Y": ["Y+", "Y-"],
        "Z": ["Z0", "Z1"],
    }

    NON_ZERO_FLOWS = [
        "II->II",
        "XI->XX",
        "IX->IX",
        "XX->XI",
        "ZI->ZI",
        "IZ->ZZ",
        "ZZ->IZ",
        "ZX->ZX",
        "YZ->XY",
        "YX->YI",
        "YI->YX",
        "YY->XZ",
        "IY->ZY",
        "XY->YZ",
        "XZ->YY",
        "ZY->IY",
    ]

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
        """

        # Set perliminary attributes
        self.samples = samples
        self.circuits = circuits.circuits

    def calc_ptm(self, only_non_zero: bool = False) -> np.ndarray:
        # Calculate the expectation values for each basis combination
        exp_vals_per_basis_comb = self._calc_entries_for_surgery(only_non_zero=only_non_zero)

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
        len_mtx = 0

        if n_qubits == 1:
            len_mtx = 3
        elif n_qubits == 2:
            len_mtx = 16

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

        final_meas: np.ndarray = np.zeros(shape=(samples_from_sampler.shape[0],), dtype=np.int8)

        for curr_rec in rec_pos:
            # XOR current measurement with the final measurement (Vectorized)
            final_meas ^= samples_from_sampler[:, curr_rec].astype(np.int8)

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

    @staticmethod
    def _get_sgn(curr_basis_combination: str) -> dict[str, int]:
        """
        For a given Basis determine if the current combination is an identity combination
        or not i.e. if a sign flip is needed for the expectation value calculation or not
        """

        # Split sign
        input_pauli = curr_basis_combination.split("->")[0]
        control_basis, target_basis = input_pauli

        # Convention used later: we want a SUM for identity and a DIFFERENCE for non-identity.
        # term = a + sgn*b  ->  sgn = +1 for sum, sgn = -1 for difference.
        return {
            "control": 1 if control_basis == "I" else -1,
            "target": 1 if target_basis == "I" else -1,
        }

    @staticmethod
    def _split_label(init_state_label: str) -> tuple[str, str]:
        """
        Splits a label such as "X+,Z0" into its control and target state components
        """

        if "," in init_state_label:
            parts = [p for p in init_state_label.split(",")]

        else:
            raise ValueError(
                f"Invalid init state label: {init_state_label!r}. "
                f"Expected format 'ControlState,TargetState'.",
            )

        return parts[0], parts[1]

    @staticmethod
    def _get_init_pairing_value(
        pair: tuple[int, int],
        average_logical_state: dict[str, float],
    ) -> float:
        """
        Given a pair of two states (e.g. (1, 0))) this function returns the corresponding
        measurement result
        """

        # Init parameters
        state_order = {"X+": 0, "X-": 1, "Y+": 0, "Y-": 1, "Z0": 0, "Z1": 1}
        control_state, target_state = pair

        for curr_init_state_label, exp_val in average_logical_state.items():
            curr_control_state, curr_target_state = PTMCalculator._split_label(
                curr_init_state_label,
            )

            if (state_order[curr_control_state] == control_state) and (
                state_order[curr_target_state] == target_state
            ):
                return exp_val

        raise ValueError(
            f"Could not find matching init state for pair: {pair}.",
        )

    def _calc_entries_for_surface_patch(self):
        pass

    def _calc_entries_for_surgery(self, only_non_zero: bool = False) -> dict[str, float]:
        """
        We can't directly use the DEM as we need the
        raw measurements to infer what logical state we have

        -> We use compile sampler to infer the logical state
        -> Convert into a dem to run the matching
        """

        # Init Save Dict for the expectation values of each basis combination
        exp_vals_per_basis_comb: dict[str, float] = {}

        # Run through diagonal circuits
        for basis_combination, circuit_dict_and_meas_recs in self.circuits.items():
            if only_non_zero and basis_combination not in self.NON_ZERO_FLOWS:
                continue

            circuit_dict, meas_rec = circuit_dict_and_meas_recs

            # For each circuit we need to calculate the expectation value
            # For Identity we take the sum of the states!
            # For X Y Z we calculate measurement(+eigenstate)
            # - measurement(-eigenstate) / total number of samples

            ####################### EXAMPLE 1 ################################
            #       IX                                                       #
            #       IX+ -> Z0 X+ (0+)                                        #
            #           -> Z1 X+ (1+)                                        #
            #                                                                #
            #       IX- -> Z0 X- (0-)                                        #
            #           -> Z1 X- (1-)                                        #
            #                                                                #
            #   We calc at the end (0+ + 1+) - (0- + 1-) / total samples     #
            ##################################################################

            ####################### EXAMPLE 2 ################################
            #       ZX                                                       #
            #       ZX+ -> Z0 X+ (0+)                                        #
            #           -> Z1 X+ (1+)                                        #
            #                                                                #
            #       ZX- -> Z0 X- (0-)                                        #
            #           -> Z1 X- (1-)                                        #
            #                                                                #
            #   We calc at the end (0+ - 1+) - (0- - 1-) / total samples     #
            ##################################################################

            average_logical_state: dict[str, int] = {}

            for curr_init_state_label, curr_circuit in circuit_dict.items():
                # Buid the normal measurement smaples and sample n shots
                sampler = curr_circuit.compile_sampler()
                results_samples = sampler.sample(shots=self.samples)

                # Getting the individual logical states after measurement (Vectorized)
                logical_states_noisy = self._get_logical_meas_from_samples(
                    results_samples,
                    meas_rec,
                )

                # Build DEM and Matcher
                matcher, dem = self._get_dem_and_matcher(circuit=curr_circuit)

                # Converting the measurement sample into DEM sample and continue as usual
                # with decoding
                detectors_observable_states = curr_circuit.compile_m2d_converter().convert(
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
                average_logical_state[curr_init_state_label] = np.average(final_logical_states)

            # Calculate Sign for the current basis combination
            sgn_dict = self._get_sgn(basis_combination)
            target_sgn = sgn_dict["target"]
            control_sgn = sgn_dict["control"]

            # Get the individual Terms for construction
            term1 = self._get_init_pairing_value(
                pair=(0, 0),
                average_logical_state=average_logical_state,
            ) + target_sgn * self._get_init_pairing_value(
                pair=(1, 0),
                average_logical_state=average_logical_state,
            )
            term2 = self._get_init_pairing_value(
                pair=(0, 1),
                average_logical_state=average_logical_state,
            ) + target_sgn * self._get_init_pairing_value(
                pair=(1, 1),
                average_logical_state=average_logical_state,
            )

            # Get Exp. Value and add to final Value Dict
            exp_val = term1 + control_sgn * term2
            exp_vals_per_basis_comb[basis_combination] = exp_val

        return exp_vals_per_basis_comb
