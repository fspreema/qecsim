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
            # -> This includes non-FT emasurement of Y OBS
            self.geometry = master_geometry.geometry_std
            self.pairings = master_pairings.pairings_std

        elif self.type == "log_h":
            # Get Geometry and Pairings for logical H final measurement
            self.geometry = master_geometry.geometry_std
            self.pairings = master_pairings.pairings_log_h

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
            final_circuit.append("MZ", self.geometry.get_logical_observables("Z"))

        elif self.geometry.obs == "X":
            final_circuit.append("MX", self.geometry.get_logical_observables("X"))

        elif self.geometry.obs == "Y":
            idx_mx, idx_my, idx_mz = self.geometry.get_logical_observables(
                "Y",
                fixed_coord=(self.geometry.distance * 2 - 1),
            )

            final_circuit.append("MX", idx_mx)
            final_circuit.append("MY", idx_my)
            final_circuit.append("MZ", idx_mz)

            print(idx_mx, idx_my, idx_mz)

        return final_circuit

    def _apply_observables(self) -> stim.Circuit:
        # Defining Observable Circuit
        observable_circuit = stim.Circuit()

        if self.geometry.obs == "Y":
            flat_obs_idx = [
                qubit_idx
                for qubit_type_list in self.geometry.get_logical_observables(self.geometry.obs)
                for qubit_idx in qubit_type_list
            ]
        else:
            flat_obs_idx = self.geometry.get_logical_observables(self.geometry.obs)

        observable_circuit.append(
            "OBSERVABLE_INCLUDE",
            [
                stim.target_rec(-i) for i in range(1, len(flat_obs_idx) + 1)
            ],
            0,
        )

        return observable_circuit
