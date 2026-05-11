import stim

from src.codes.surface_code_rotated.data_geometry import MasterGeometry, MasterPairings
from src.core.cx_builder import cx_builder
from src.core.measurement_tracker import MeasurementTracker

Coord = complex

__all__ = ["SurfaceRepetitionCircuit"]


class SurfaceRepetitionCircuit:
    def __init__(
        self,
        master_geometry: MasterGeometry,
        master_pairings: MasterPairings,
        type: str,
        tracker: MeasurementTracker,
    ):
        """
        Initialize the Surface Repetition Circuit

        Parameters:
            geometry : SurfaceGeometry
            pairings : SurfacePairings
            noise : NoiseParameters
                -> Depending on the type the correct pairings need to be loaded in (Surface/Memory)!
            type : str, default "standard"
                -> What type of repetition, i.e. standard (x,z basis) or y-basis repetition/ memory

        """

        if type not in {"standard", "y_memory", "y_basis", "log_h", "non_ft_init"}:
            raise ValueError(
                f"Unknown repetition circuit type: {type}. "
                f"Type must be one of 'standard', 'y_memory', 'y_basis', 'log_h', 'non_ft_init'.",
            )

        self.type = type
        self.tracker = tracker

        # Initialize Geometry and Pairings depending on the type
        if self.type == "y_memory":
            # Get Geometry and Pairings for y-basis memory round
            self.geometry = master_geometry.geometry_ybasis
            self.pairings = master_pairings.pairings_ymemory

            # Define Indexes for memory round
            self.stab_idx = self.geometry.stab_x_memory + self.geometry.stab_z_memory
            self.stab_x_idx = self.geometry.stab_x_memory
            self.stab_z_idx = self.geometry.stab_z_memory

            # Set number of rounds
            self.rounds = self.geometry.distance - 1

        elif self.type == "y_basis":
            # Get Geometry and Pairings for standard repetition or y-basis repetition
            self.geometry = master_geometry.geometry_ybasis
            self.pairings = master_pairings.pairings_ybasis

            # Define Indexes for standard repetition
            self.stab_idx = self.geometry.stab_x_idx + self.geometry.stab_z_idx
            self.stab_x_idx = self.geometry.stab_x_idx
            self.stab_z_idx = self.geometry.stab_z_idx

            # Set number of rounds
            self.rounds = int((self.geometry.distance - 1) / 2)

        elif self.type == "log_h":
            # Get Geometry and Pairings for logical H repetition round
            self.geometry = master_geometry.geometry_std
            self.pairings = master_pairings.pairings_log_h

            # Define Indexes for memory round
            self.stab_idx = self.geometry.stab_x_idx + self.geometry.stab_z_idx
            self.stab_x_idx = self.geometry.stab_x_idx
            self.stab_z_idx = self.geometry.stab_z_idx

            # Set number of rounds
            self.rounds = self.geometry.distance - 1

        elif self.type in {"standard", "non_ft_init"}:
            # Get Geometry and Pairings for standard repetition or y-basis repetition
            self.geometry = master_geometry.geometry_std
            self.pairings = master_pairings.pairings_std

            # Define Indexes for standard repetition
            self.stab_idx = self.geometry.stab_x_idx + self.geometry.stab_z_idx
            self.stab_x_idx = self.geometry.stab_x_idx
            self.stab_z_idx = self.geometry.stab_z_idx

            # Set number of rounds
            self.rounds = self.geometry.distance - 1

    def build_circuit(self) -> stim.Circuit:
        circuit = stim.Circuit()

        # If normal repetition (Non y-basis)
        if self.type in {"standard", "log_h", "y_basis", "non_ft_init"}:
            circuit += self._adding_repetition_rounds()

        elif self.type == "y_memory":
            circuit += self._y_basis_memory_prep_circuit()
            circuit += self._adding_repetition_rounds()
            circuit += self._y_basis_add_non_det_obs()[0]

        return circuit

    def rec_list(self):
        if self.type == "y_memory":
            return self._y_basis_add_non_det_obs()[1]
        return []
    
    def _get_detectors(
        self,
        measured_qubits: list[int],
        patch_type: str,
        qubits_for_detectors: list[int] | None = None,
    ) -> stim.Circuit:
        # Initialization Circuit
        detector_circuit = stim.Circuit()

        self.tracker.add_measurements_to_tracker(
            measured_qubits=measured_qubits,
            qubits_for_detectors=qubits_for_detectors,
            patch_type=patch_type,
        )

        det_record_pairings = self.tracker.get_records_for_detectors(
            patch_type=patch_type,
        )
        for curr_pairing in det_record_pairings:
            detector_circuit.append("DETECTOR", curr_pairing)

        # Shifting Coords
        detector_circuit.append("SHIFT_COORDS")

        return detector_circuit

    def _y_basis_memory_prep_circuit(self) -> stim.Circuit:
        # Init reset Circuit
        y_memory_prep_circ = stim.Circuit()

        y_memory_prep_circ.append("R", self.stab_idx)
        y_memory_prep_circ.append("TICK")

        if self.geometry.state_init == "-i":
            # Get logical y string
            x_idx, y_idx, z_idx = self.geometry.get_logical_observables("Y")

            # Append logical flip of Y Observable
            y_memory_prep_circ.append("X", x_idx)
            y_memory_prep_circ.append("Y", y_idx)
            y_memory_prep_circ.append("Z", z_idx)
            y_memory_prep_circ.append("TICK")

        y_memory_prep_circ.append("H", self.stab_x_idx)
        y_memory_prep_circ.append("TICK")

        # 2) CX Operations
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.pairings.stab_to_data,
            circuit=y_memory_prep_circ,
        )

        # 3) Basis/ Measurement
        y_memory_prep_circ.append("H", self.stab_x_idx)
        y_memory_prep_circ.append("TICK")
        y_memory_prep_circ.append("M", self.stab_x_idx + self.stab_z_idx)
        y_memory_prep_circ.append("SHIFT_COORDS", arg=(0, 0, 1))
        y_memory_prep_circ.append("TICK")

        return y_memory_prep_circ

    def _y_basis_add_non_det_obs(self) -> tuple[stim.Circuit, list[int]]:
        observable_circ = stim.Circuit()

        if self.geometry.obs == "X":
            # Getting corresponding logical string and rec
            log_x = self.geometry.get_logical_observables(
                "X",
                fixed_coord=(self.geometry.distance * 2 - 1),
            )

            observable_circ.append("MX", log_x)
            observable_circ.append("OBSERVABLE_INCLUDE", [f"X{index}" for index in log_x], 0)

            # For later decoding we need the measurement record postiions of the logical operator
            rec_list = [-i - 1 for i in range(len(log_x))]

            return observable_circ, rec_list

        elif self.geometry.obs == "Z":
            # Getting corresponding logical string and rec
            log_z = self.geometry.get_logical_observables(
                "Z",
                fixed_coord=(self.geometry.distance * 2 - 1),
            )

            observable_circ.append("MZ", log_z)
            observable_circ.append("OBSERVABLE_INCLUDE", [f"Z{index}" for index in log_z], 0)

            # For later decoding we need the measurement record postiions of the logical operator
            rec_list = [-i - 1 for i in range(len(log_z))]

            return observable_circ, rec_list

        return observable_circ, []

    def _adding_repetition_rounds(self) -> stim.Circuit:
        repetition_circ = stim.Circuit()

        # We need to do *3 rounds as we currently only do non-ft and therefore need
        # d rounds prioir and d rounds after which are nosieless. To ensure that for every
        # d we still have d noisy rounds we do 3*d rounds in total...
        for curr_round in range((self.rounds * 3) - 1):

            # -----BUILDING-REPETITION-CIRC------
            repetition_circ.append("R", self.stab_idx)
            repetition_circ.append("TICK")

            # 1) Reset/ Basis
            repetition_circ.append("H", self.stab_x_idx)
            repetition_circ.append("TICK")

            # 2) CX Operations
            cx_builder(
                q2i=self.geometry.q2i,
                stab_to_data=self.pairings.stab_to_data,
                circuit=repetition_circ,
                excluded_index=self.geometry.y_index if self.type == "y_basis" else None,
            )

            # 3) Basis/ Measurement
            repetition_circ.append("H", self.stab_x_idx)
            repetition_circ.append("TICK")
            repetition_circ.append("M", self.stab_idx)

            # 4) Adding Detectors
            repetition_circ += self._get_detectors(
                measured_qubits=self.stab_idx,
                patch_type="STD_PATCH",
                qubits_for_detectors= []
                if curr_round == 0
                else self.stab_idx,
            )

            repetition_circ.append("TICK")

        return repetition_circ
