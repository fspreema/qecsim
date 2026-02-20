import stim

from src.core.base_geometry import BaseGeometry


class CircuitChunk:
    def __init__(self, chunk_circuit: stim.Circuit, geometry: BaseGeometry):
        """
        Circuit Chunk with automated flow extraction including measurements for detectors.
        """
        self.circuit = chunk_circuit
        self.num_measurements = chunk_circuit.num_measurements
        self.coords = geometry.coords
        self.q2i = geometry.q2i
        self.i2q = geometry.i2q
        self.flows = list(self._extract_flows())

    def get_annihilation_flows(self):
        """
        Get all annihilation flows from the circuit chunk.

        Returns:
            A list of flow dictionaries corresponding to annihilation flows.
        """

        annihilation_flows = [flow for flow in self.flows if flow["type"] == "annihilation"]

        return annihilation_flows

    def get_creation_flows(self):
        """
        Get all creation flows from the circuit chunk.

        Returns:
            A list of flow dictionaries corresponding to creation flows.
        """

        creation_flows = [flow for flow in self.flows if flow["type"] == "creation"]

        return creation_flows

    @staticmethod
    def get_pauli_indices(current_flow: stim.Flow):
        """
        Get the pauli string indices from the flow.

        Returns:
            A List of indexes comprising the pauli string (Type not included).
        """

        string_repr = str(current_flow)
        left, right = string_repr.split("->")

        # Extract pauli string indices from the left side of the flow
        pauli_indices = []

        for element in left.split():
            if any(pauli in element for pauli in "XYZ"):
                # Find position of the Pauli operator
                for char_pos, char in enumerate(element):
                    if char in "XYZ":
                        pauli_indices.append(char_pos)

        for element in right.split():
            if any(pauli in element for pauli in "XYZ"):
                # Find position of the Pauli operator
                for char_pos, char in enumerate(element):
                    if char in "XYZ":
                        pauli_indices.append(char_pos)

        return pauli_indices

    @staticmethod
    def _count_flow_weight(current_flow: stim.Flow):
        """
        Count the number of non Identity Pauli strings in flow evolution.

        Returns:
            A tuple (num_non_identity_paulis_before, num_non_identity_paulis_after)
        """

        num_non_identity_paulis_before = 0
        num_non_identity_paulis_after = 0

        string_repr = str(current_flow)
        left, right = string_repr.split("->")

        num_non_identity_paulis_before += sum(left.count(c) for c in "XYZ")
        num_non_identity_paulis_after += sum(right.count(c) for c in "XYZ")

        return num_non_identity_paulis_before, num_non_identity_paulis_after

    @staticmethod
    def _classify_flow_type(current_flow: stim.Flow):
        """
        Classify the flow type based on its pauli string before and after.

        Returns:
            A string indicating the flow type: "creation", "annihilation", "passthrough", or "other"
        """

        string_repr = str(current_flow)
        left, right = string_repr.split("->")

        has_pauli_before = any(c in left for c in "XYZ")
        has_pauli_after = any(c in right for c in "XYZ")

        if not has_pauli_before and has_pauli_after:
            return "creation"
        elif has_pauli_before and not has_pauli_after:
            return "annihilation"
        elif has_pauli_before and has_pauli_after:
            return "passthrough"
        else:
            return "other"

    def _get_measurements_for_flow(self, current_flow: stim.Flow):
        """
        Get all measurements that are part of the flow.

        Returns:
            A list of measurement indices.
        """

        try:
            (included_measurements,) = self.circuit.solve_flow_measurements([current_flow])

            converted_measurements = []

            # Converting into detector compatible measurements recs
            for curr_measurement in included_measurements:
                converted_measurements.append(-self.num_measurements + curr_measurement)

            return converted_measurements

        except (RuntimeError, ValueError) as e:
            print(f"Could not get measurements for flow {current_flow}: {e}")
            return []

    def _check_locality(self, current_flow: stim.Flow, threshold: complex = 2 + 2j):
        """
        Check if all coordinates in the flow are within a certain threshold.

        Args:
            threshold (complex): The maximum allowed distance for both the real and imaginary parts
            separately.
            The threshold is applied independently to the real and imaginary components of the
            coordinates.
        """

        # Get the index of all corresponding pauli strings
        current_indices = self.get_pauli_indices(current_flow)

        currently_checking_coords: complex | None = None

        # Convert into Coordinates
        for index in current_indices:
            coord = self.i2q[index]

            if currently_checking_coords is None:
                currently_checking_coords = coord

            elif coord != currently_checking_coords:
                if abs(coord.real - currently_checking_coords.real) > threshold.real or (
                    abs(coord.imag - currently_checking_coords.imag) > threshold.imag
                ):
                    return False

        return True

    def _extract_flows(self, allowed_weights: list[int] = None):
        """
        Uses the flow generator function in order to get all flows up to a certain weight.
        """

        # Initializing allowed weights if None
        if allowed_weights is None:
            allowed_weights = [2, 4]

        try:
            # Retrieving possible flows from the circuit
            possible_flows = self.circuit.flow_generators()

            for current_flow in possible_flows:
                # Filter for flows with correct weights and creation/annihilation properties
                weight_incoming = self._count_flow_weight(current_flow)[0]
                weight_outgoing = self._count_flow_weight(current_flow)[1]
                flow_type = self._classify_flow_type(current_flow)

                # Filter based on flow type and allowed weights:
                # - Creation: 1 -> X..X
                # - Annihilation: X..X -> 1
                if flow_type == "creation" and weight_incoming == 0:
                    # Check if outgoing weight is in allowed range
                    if weight_outgoing in allowed_weights:
                        # Check Locality
                        if not self._check_locality(current_flow):
                            continue

                        # Get measurement indices for the flow
                        measurements_in_flow = self._get_measurements_for_flow(current_flow)

                        # Yield the flow dictionary
                        yield {
                            "flow": current_flow,
                            "measurements": measurements_in_flow,
                            "type": flow_type,
                            "weight": weight_outgoing,
                        }

                elif flow_type == "annihilation" and weight_outgoing == 0:
                    # Check if incoming weight is in allowed range
                    if weight_incoming in allowed_weights:
                        # Check Locality
                        if not self._check_locality(current_flow):
                            continue

                        # Get measurement indices for the flow
                        measurements_in_flow = self._get_measurements_for_flow(current_flow)

                        # Yield the flow dictionary
                        yield {
                            "flow": current_flow,
                            "measurements": measurements_in_flow,
                            "type": flow_type,
                            "weight": weight_incoming,
                        }

        except (RuntimeError, ValueError) as e:
            if "current_flow" in locals():
                print(f"Could not extract flow {current_flow}: {e}")
            else:
                print(f"Could not extract flows from circuit: {e}")


class CompileChunk:
    def __init__(self, chunk_circuits: list[CircuitChunk]):
        """
        Compilation of circuit by adding necessary Detector instructions
        """
        self.chunks = chunk_circuits
        self.compiled_circuit = None

    def add_chunk(self, chunk: CircuitChunk):
        """
        Add a CircuitChunk to the compilation process.

        Args:
            chunk (CircuitChunk): The CircuitChunk to be added.
        """
        self.chunks.append(chunk)

    def compile(self) -> stim.Circuit:
        """
        Compile the circuit by adding Detector instructions based on measurement outcomes.

        Logic:
        - Each chunk should have creation and annihilation flows
        - each of these types are used and the measurements are extracted
        - The detector compares measurements from creation with those from annihilation
        """

        compiled_circ = stim.Circuit()
        previous_chunk = None
        measurement_offset = 0

        for i, curr_chunk in enumerate(self.chunks):
            # Adding Circuit to compiled circuit
            compiled_circ += curr_chunk.circuit
            measurement_offset = curr_chunk.num_measurements

            if i == 0:
                previous_chunk = curr_chunk

            # Second chunk only has creation flows (No matcher needed!)
            # -> Need to be compatible with the initlized state and therefore reset
            #    needs to be included
            if i == 1:
                creation_flows = [flow for flow in curr_chunk.flows if flow["type"] == "creation"]

                for flow in creation_flows:
                    measurements = flow["measurements"]

                    # Check if measurements exist
                    if not measurements:
                        continue

                    # Add Detector instruction for creation flows
                    compiled_circ.append_operation(
                        "DETECTOR",
                        [stim.target_rec(m) for m in measurements],
                    )

                # Adding previous Chunk
                previous_chunk = curr_chunk
                measurement_offset = curr_chunk.num_measurements

            # Subsequent chunks have both creation and annhilation flows
            elif i < len(self.chunks) - 1:
                # Cross reference flows to find matches
                curr_matches = self._flow_matcher(previous_chunk, curr_chunk)

                for a_flow, c_flow in curr_matches:
                    measurement_prev = [
                        meas_idx - measurement_offset for meas_idx in a_flow["measurements"]
                    ]
                    measurement_curr = c_flow["measurements"]

                    # Check if measurements exist
                    if not measurement_prev or not measurement_curr:
                        continue

                    # Add Detector instruction comparing creation and annihilation flows
                    compiled_circ.append_operation(
                        "DETECTOR",
                        [stim.target_rec(m) for m in measurement_curr]
                        + [stim.target_rec(m) for m in measurement_prev],
                    )

                # Shift coords after Chunk is finished
                compiled_circ.append_operation("SHIFT_COORDS", arg=(0, 0, 1))

                # Update previous chunk and measurement offset
                previous_chunk = curr_chunk
                measurement_offset += curr_chunk.num_measurements

            # Last chunk only has annihilation flows
            else:
                """
                Detectors may are included in the last round but this heavily depends on the 
                measurements of the data qubits.
                """

                annihilation_flows = [
                    flow
                    for flow in (previous_chunk.flows + curr_chunk.flows)
                    if flow["type"] == "annihilation"
                ]

                for flow in annihilation_flows:
                    measurements = flow["measurements"]

                    # Check if measurements exist
                    if not measurements:
                        continue

                    # Add Detector instruction for annihilation flows
                    compiled_circ.append_operation(
                        "DETECTOR",
                        [stim.target_rec(m) for m in measurements],
                    )

        return compiled_circ

    @staticmethod
    def _flow_matcher(prev_chunk: CircuitChunk, curr_chunk: CircuitChunk):
        """
        Match creation and annihilation flows based on pauli string indices.

        - X.... -> 1 matches 1 -> X.... etc.

        Args:
            prev_chunk (CircuitChunk): The previous CircuitChunk containing annihilation flows.
            curr_chunk (CircuitChunk): The current CircuitChunk containing creation flows.

        Returns:
            list[tuple[dict, dict]]: A list of tuples, each containing a pair of matching
            annihilation and creation flow dictionaries.
        """

        # Flow pairings
        flow_pairings = []

        # Get annihilation flow indices from previous chunk
        flow_annihilation_indices = prev_chunk.get_annihilation_flows()
        flow_creation_indices = curr_chunk.get_creation_flows()

        # Match creation and annihilation flows based on their indices
        for a_flow in flow_annihilation_indices:
            for c_flow in flow_creation_indices:
                if prev_chunk.get_pauli_indices(a_flow["flow"]) == curr_chunk.get_pauli_indices(
                    c_flow["flow"],
                ):
                    # Yield the flow dictionary
                    flow_pairings.append((a_flow, c_flow))

        return flow_pairings

    def _extract_coords_from_pauli(self, pauli_string: str):
        """
        Extract coordinates from a given pauli string.

        Args:
            pauli_string (str): The pauli string to extract coordinates from.
        """
