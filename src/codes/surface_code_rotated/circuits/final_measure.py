import stim

from src.codes.surface_code_rotated.data_geometry import MasterGeometry, MasterPairings

__all__ = ["FinalMeasureCircuit"]


class FinalMeasureCircuit:
    def __init__(
        self,
        master_geometry: MasterGeometry,
        master_pairings: MasterPairings,
        type: str,
    ):
        """
        Initialize the Final Measurement Circuit

        Parameters:
            geometry : SurfaceGeometry
            pairings : SurfacePairings
            noise : NoiseParameters
                -> Depending on the type the correct pairings need to be loaded in (Surface/Memory)!
            type : str, default "standard"
                -> What type of final measurement, i.e. standard (x,z basis) or y-basis measurement
        """

        if type not in {"standard", "log_h"}:
            raise ValueError(
                f"Unknown final measurement circuit type: {type}. "
                f"Type must be one of 'standard', 'log_h'.",
            )

        self.type = type

        # Initialize Geometry and Pairings depending on the type
        if self.type == "standard":
            # Get Geometry and Pairings for standard final measurement
            self.geometry = master_geometry.geometry_std
            self.pairings = master_pairings.pairings_std

        elif self.type == "log_h":
            # Get Geometry and Pairings for logical H final measurement
            self.geometry = master_geometry.geometry_std
            self.pairings = master_pairings.pairings_log_h

        # Defining Non deterministic pairing flag
        self.non_deterministic_pairing = (
            False
            if (
                (self.geometry.state_init in {"0", "1"} and self.geometry.obs == "Z")
                or (self.geometry.state_init in {"+", "-"} and self.geometry.obs == "X")
                or (self.geometry.state_init in {"+i", "-i"} and self.geometry.obs == "Y")
            )
            else True
        )

    def build_final_measurement_circuit(self) -> stim.Circuit:
        final_circuit = stim.Circuit()

        # Building Circuit
        final_circuit += self._apply_final_measurement()
        final_circuit += self._apply_observables()

        return final_circuit

    def build_observable_meas_rec(self) -> list[int]:
        # Defining Record List
        rec_list: list[int] = []

        # Flatten Observale indice (Needed for Y Obs)
        # [[...], [...], [.]] -> [.....]
        if self.geometry.obs == "Y":
            qubit_idx = [
                qubit_idx
                for qubit_type_list in self.geometry.get_logical_observables(self.geometry.obs)
                for qubit_idx in qubit_type_list
            ]
        else:
            qubit_idx = self.geometry.get_logical_observables(self.geometry.obs)

        # Return Measurement Records of the logical operator for later decoding
        rec_list = [
            -len(self._what_qubits_measured()) + k
            for k in self._get_rec_targets(
                indicies=qubit_idx,
                measured_qubits=self._what_qubits_measured(),
            )
        ]

        return rec_list

    def _apply_final_measurement(self) -> stim.Circuit:
        ##########################
        # Adding final measurement
        ##########################

        final_circuit = stim.Circuit()

        if self.geometry.obs == "Z":
            final_circuit.append("MZ", self.geometry.data_idx)

        elif self.geometry.obs == "X":
            final_circuit.append("MX", self.geometry.data_idx)

        elif self.geometry.obs == "Y" and self.non_deterministic_pairing:
            idx_mx, idx_my, idx_mz = self.geometry.get_logical_observables(
                "Y",
                fixed_coord=(self.geometry.distance * 2 - 1),
            )

            final_circuit.append("MX", idx_mx)
            final_circuit.append("MY", idx_my)
            final_circuit.append("MZ", idx_mz)

        return final_circuit

    @staticmethod
    def _get_rec_targets(indicies: list[int], measured_qubits: list[int]) -> list[int]:
        """
        Finds the record targets for given qubit indicies

        Parameters:
            indicies : list[int]
                -> List of qubit indicies to find rec targets for
            measured_qubits : list[int]
                -> List of all measured qubits to find correct record targets from
        """

        tar_rec = []

        for rec_pos, index in enumerate(measured_qubits):
            if index in indicies:
                tar_rec.append(rec_pos)

        return tar_rec

    def _what_qubits_measured(self) -> list[int]:
        """
        Finds which qubits are measured depending on the observable type

        Returns:
            list[int] : List of measured qubit indicies
        """

        if self.geometry.obs in {"X", "Z"}:
            measured_qubits = self.geometry.data_idx
        else:
            measured_qubits = (
                self.geometry.get_logical_observables("Y")[0]
                + self.geometry.get_logical_observables("Y")[1]
                + self.geometry.get_logical_observables("Y")[2]
            )

        return measured_qubits

    def _apply_observables(self) -> stim.Circuit:
        # Defining Observable Circuit
        observable_circuit = stim.Circuit()

        # Defining Logical Observables
        if self.non_deterministic_pairing:
            # Non-deterministic observables are added by Pauli Indexes
            if self.geometry.obs in {"X", "Z"}:
                log_indices = self.geometry.get_logical_observables(self.geometry.obs)

                observable_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [f"{self.geometry.obs}{index}" for index in log_indices],
                    0,
                )

            elif self.geometry.obs in {"Y"}:
                log_x, log_y, log_z = self.geometry.get_logical_observables(
                    "Y",
                    fixed_coord=(self.geometry.distance * 2 - 1),
                )

                observable_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [f"X{index}" for index in log_x]
                    + [f"Y{index}" for index in log_y]
                    + [f"Z{index}" for index in log_z],
                    0,
                )

        # As the logical Y obsesrvable inside the y basis is measured by stabilizers in the reverse
        # switch, we do not need to add any observable here
        elif not self.non_deterministic_pairing and self.geometry.obs in {"X", "Z"}:
            # Determinstic observables implemented by real rec targets
            observable_circuit.append(
                "OBSERVABLE_INCLUDE",
                [
                    stim.target_rec(-len(self._what_qubits_measured()) + i)
                    for i in self._get_rec_targets(
                        indicies=self.geometry.get_logical_observables(self.geometry.obs),
                        measured_qubits=self._what_qubits_measured(),
                    )
                ],
                0,
            )

        return observable_circuit
