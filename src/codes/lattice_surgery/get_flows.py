import stim

from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry

__all__ = ["SurgeryFlowObservables"]


class SurgeryFlowObservables:
    def __init__(
        self,
        geometry: SurgeryGeometry,
    ):
        # Initialize Geometry
        self.geometry = geometry

    def get_observable_from_flow(
        self,
        flow_circuit: stim.Circuit,
        flow_type: str,
    ) -> stim.Circuit:
        """
        For a given circuit without final measurement and reset, the needed Measurements can be
        determined by the flow type.

        -> These are given as a list of qubit indices to be measured.
        -> These are then applied to the circuit as additional measurement
        records to the same observable

        Returns:
            stim.Circuit: Corrected Circuit with the correct logical observable included
        """

        # Init Return Circuit
        return_circuit = stim.Circuit()

        # Get Pauli Strings for Flow
        pauli_start, pauli_end = self._get_pauli_strings_for_flows(flow_type=flow_type)

        # Construct Full Pauli Strings
        start_string = self._construct_pauli_string(logical_operator_strings=pauli_start)
        end_string = self._construct_pauli_string(logical_operator_strings=pauli_end)

        # Determine Flow Measurements
        (included_measurements,) = flow_circuit.solve_flow_measurements(
            [
                stim.Flow(f"{start_string} -> {end_string}"),
            ],
        )

        # Print Flow Info
        print("Logical Flow from Circuit:")
        print(start_string, "->", end_string)

        # Calculating target rec pos
        rec_pos = []

        # Adding Measurements to the Observable if solution exists
        try:
            for index in included_measurements:
                current_rec_tar = flow_circuit.num_measurements - index
                rec_pos.append(-current_rec_tar)

        except TypeError:
            print(
                "No Logical Observable for given Flow Type, "
                "circuit returned without correct observable.",
            )
            print(
                "These are the available flows for the current circuit:",
            )
            self._debug_print_available_flows(flow_circuit=flow_circuit)

        # Adding measurements to the logical observable
        return_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(k) for k in rec_pos], 0)

        return return_circuit

    def _get_pauli_strings_for_flows(
        self,
        flow_type: str,
    ) -> list[list[str]]:
        """
        Returns the Pauli string from the creation and the end of the flow type
        -> i.e. start -> end of flow
        -> This describes how the pauli string should begin and propagate through the circuit
        """

        # Define Dictionary for Flow Types
        flow_dict = {
            # Non Mixed Logical Strings
            "XI -> XX": [["c_x"], ["c_x", "t_x"]],
            "XX -> XI": [["c_x", "t_x"], ["c_x"]],
            "IX -> IX": [["t_x"], ["t_x"]],
            "IZ -> ZZ": [["t_z"], ["c_z", "t_z"]],
            "ZZ -> IZ": [["c_z", "t_z"], ["t_z"]],
            "ZI -> ZI": [["c_z"], ["c_z"]],
            # Mixed Logical Strings
            "ZX -> ZX": [["c_z", "t_x"], ["c_z", "t_x"]],
            # Y Logical Strings
            "YZ -> XY": [["t_z_shifted"], ["c_x_shifted", "t_y"]],
            "YX -> YI": [["t_x_shifted"], ["c_y"]],
            "YI -> YX": [[], ["t_x_shifted", "c_y"]],
            "YY -> XZ": [[], ["c_x_shifted", "t_z_shifted"]],
            "IY -> ZY": [[], ["c_z_shifted", "t_y"]],
            "XY -> YZ": [["c_x_shifted"], ["c_y", "t_z_shifted"]],
        }

        return flow_dict.get(flow_type, [[], []])

    def _construct_pauli_string(self, logical_operator_strings: list[str]) -> str:
        """
        Returns the full pauli string from the logical operator strings
        Args:
            logical_operator_strings (list[str]): List of logical operator strings
                e.g. ["c_x", "t_z", "t_x_shifted"]
        """

        # If no logical operator strings are given, return identity
        if logical_operator_strings == []:
            return "1"

        # Get Logical Strings
        log_strings = self.geometry.get_logical_strings()
        log_strings_shifted = self.geometry.get_logical_strings(
            shift_cx_for_y=True,
            shift_tz_for_y=True,
            shift_tx_for_y=True,
            shift_cz_for_y=True,
        )

        terms = []

        # Identify basis based on key endings
        for key in logical_operator_strings:
            if key.endswith("_shifted"):
                # Determine basis
                basis = "X" if key.endswith("_x_shifted") else "Z"

                # remove _shifted for lookup
                key_clean = key[:-8]

                # Get indices from shifted logical strings
                indices = log_strings_shifted.get(key_clean, [])

                terms.extend(f"{basis}{i}" for i in indices)

            elif key.endswith("_y"):
                # Get indices
                indices_x = log_strings[key]["x_string"]
                index_y = log_strings[key]["y_corner"]
                indices_z = log_strings[key]["z_string"]

                # Format terms
                terms.extend(f"Z{i}" for i in indices_z)
                terms.extend(f"Y{i}" for i in index_y)
                terms.extend(f"X{i}" for i in indices_x)

            else:
                basis = "X" if key.endswith("_x") else "Z"
                # Get indices
                indices = log_strings.get(key, [])
                # Format terms
                terms.extend(f"{basis}{i}" for i in indices)

        return "*".join(terms)

    def _debug_print_available_flows(self, flow_circuit: stim.Circuit):
        """
        print all available flows for the current circuit provided
        """

        for flows in flow_circuit.flow_generators():
            print(flows)
