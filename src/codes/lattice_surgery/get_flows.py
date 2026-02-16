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
        curr_flow: str,
    ) -> stim.Circuit:
        """
        For a given circuit without final measurement and reset, the needed Measurements can be
        determined by the flow type.

        -> These are given as a list of qubit indices to be measured.
        -> These are then applied to the circuit as additional measurement
        records to the same observable

        If We are in the Y Basis the Flows are calculated by the X and Z flows combined.
        -> So X and Z flows are calculated seperatly and then combined into one observable

        Returns:
            stim.Circuit: Corrected Circuit with the correct logical observable included
        """

        # Init Return Circuit
        return_circuit = stim.Circuit()

        # Get Pauli Strings for Flow
        if "Y" in curr_flow:
            # Get Info out of Flow Dictionary
            flow_pauli_strings = self._get_pauli_strings_from_flows(flow_type=curr_flow)
            pauli_start_x, pauli_end_x = flow_pauli_strings[0]
            pauli_start_z, pauli_end_z = flow_pauli_strings[1]

            # Construct Full Pauli Strings
            start_string_x = self._construct_pauli_string(logical_operator_strings=pauli_start_x)
            end_string_x = self._construct_pauli_string(logical_operator_strings=pauli_end_x)
            start_string_z = self._construct_pauli_string(logical_operator_strings=pauli_start_z)
            end_string_z = self._construct_pauli_string(logical_operator_strings=pauli_end_z)

            # Determine Flow sign -> YY -> -XZ
            flow_sign = ""
            if curr_flow == "YY -> XZ":
                flow_sign = "-"

            # Determine Flow Measurements
            (included_measurements_x,) = flow_circuit.solve_flow_measurements(
                [
                    stim.Flow(f"{start_string_x} -> {flow_sign}{end_string_x}"),
                ],
            )

            (included_measurements_z,) = flow_circuit.solve_flow_measurements(
                [
                    stim.Flow(f"{start_string_z} -> {end_string_z}"),
                ],
            )

            # Full list of measurements
            if included_measurements_x is None or included_measurements_z is None:
                full_measurements = None
            else:
                full_measurements = included_measurements_x + included_measurements_z

        else:
            pauli_start, pauli_end = self._get_pauli_strings_from_flows(flow_type=curr_flow)

            # Construct Full Pauli Strings
            start_string = self._construct_pauli_string(logical_operator_strings=pauli_start)
            end_string = self._construct_pauli_string(logical_operator_strings=pauli_end)

            # Determine Flow Measurements
            (full_measurements,) = flow_circuit.solve_flow_measurements(
                [
                    stim.Flow(f"{start_string} -> {end_string}"),
                ],
            )

        # Calculating target rec pos
        rec_pos = []

        # Adding Measurements to the Observable if solution exists
        try:
            for index in full_measurements:
                current_rec_tar = flow_circuit.num_measurements - index
                rec_pos.append(-current_rec_tar)

        except TypeError:
            # No Flow Found although Flows calculated here are valid!
            print(
                "No Flow Found for "
                + curr_flow
                + " although this flow is valid! Check if flow_circuit is correct.",
            )

        # Adding measurements to the logical observable
        return_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(k) for k in rec_pos], 0)

        return return_circuit

    @staticmethod
    def _get_pauli_strings_from_flows(
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
            "XI -> XX": [["c_x_shifted"], ["c_x_shifted", "t_x_shifted"]],
            "XX -> XI": [["c_x_shifted", "t_x_shifted"], ["c_x_shifted"]],
            "IX -> IX": [["t_x_shifted"], ["t_x_shifted"]],
            "IZ -> ZZ": [["t_z_shifted"], ["c_z_shifted", "t_z_shifted"]],
            "ZZ -> IZ": [["c_z_shifted", "t_z_shifted"], ["t_z_shifted"]],
            "ZI -> ZI": [["c_z_shifted"], ["c_z_shifted"]],
            # Mixed Logical Strings
            "ZX -> ZX": [["c_z_shifted", "t_x_shifted"], ["c_z_shifted", "t_x_shifted"]],
        }

        # Define Y flow dict by just adding XZ flows
        flow_dict_y = {
            "YZ -> XY": [
                [["c_x_shifted"], ["c_x_shifted", "t_x_shifted"]],
                [["c_z_shifted", "t_z_shifted"], ["t_z_shifted"]],
            ],
            "YX -> YI": [
                [["c_x_shifted", "t_x_shifted"], ["c_x_shifted"]],
                [["c_z_shifted"], ["c_z_shifted"]],
            ],
            "YI -> YX": [
                [["c_x_shifted"], ["c_x_shifted", "t_x_shifted"]],
                [["c_z_shifted"], ["c_z_shifted"]],
            ],
            "YY -> XZ": [
                [["c_x_shifted", "t_x_shifted"], ["c_x_shifted"]],
                [["c_z_shifted", "t_z_shifted"], ["t_z_shifted"]],
            ],
            "IY -> ZY": [
                [["t_x_shifted"], ["t_x_shifted"]],
                [["t_z_shifted"], ["c_z_shifted", "t_z_shifted"]],
            ],
            "XY -> YZ": [
                [["c_x_shifted", "t_x_shifted"], ["c_x_shifted"]],
                [["t_z_shifted"], ["c_z_shifted", "t_z_shifted"]],
            ],
        }

        # Determine which flow dict to use
        if flow_type in flow_dict_y:
            return flow_dict_y[flow_type]
        elif flow_type in flow_dict:
            return flow_dict[flow_type]
        else:
            return [[], []]

    def _construct_pauli_string(self, logical_operator_strings: list[str]) -> str:
        """
        Returns the full pauli string from the logical operator strings
        Args:
            logical_operator_strings (list[str]): List of logical operator strings
                e.g. ["c_x", "t_z", "t_x_shifted"]
        """

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

    @staticmethod
    def _get_flow_generators(
        flow_circuit: stim.Circuit,
        must_have: list[str] | None = None,
    ) -> None:
        """
        Returns all available flows or if must_have is specified only those containing all the
        strings listed
        """

        available_flows = flow_circuit.flow_generators()
        for flow in available_flows:
            if must_have is not None:
                if all(item in str(flow) for item in must_have):
                    print(flow)
            else:
                print(flow)
