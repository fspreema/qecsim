import stim

from src.codes.lattice_surgery.data_geometry import MasterPairings
from src.codes.lattice_surgery.measurement_tracker import MeasurementTracker
from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.core.cx_builder import cx_builder

Coord = complex

__all__ = ["SurgerySplit"]


class SurgerySplit:
    def __init__(
        self,
        geometry: SurgeryGeometry,
        master_pairings: MasterPairings,
        split_type: str,
        tracker: MeasurementTracker,
    ):
        # Preliminary Setup
        self.geometry = geometry
        self.master_pairings = master_pairings
        self.tracker = tracker
        self.log_strings = self.geometry.get_logical_strings()

        # Getting Specific Splitting Type Info
        self.split_type = split_type
        if split_type == "AT":
            self.x_stab_index_untouched_circ = self.geometry.control_x_stb_idx
            self.z_stab_index_untouched_circ = self.geometry.control_z_stb_idx
            self.combined_x_stab_merging_lattices, self.combined_z_stab_merging_lattices = (
                self.geometry.get_combined_xz_stabs_merging_lattice(split_type="AT")
            )
        elif split_type == "AC":
            self.x_stab_index_untouched_circ = self.geometry.target_x_stb_idx
            self.z_stab_index_untouched_circ = self.geometry.target_z_stb_idx
            self.combined_x_stab_merging_lattices, self.combined_z_stab_merging_lattices = (
                self.geometry.get_combined_xz_stabs_merging_lattice(split_type="AC")
            )
        else:
            raise ValueError("No valid splitting Type in Function selected!")

    def build_circuit(self) -> stim.Circuit:
        # Init return Circuit
        return_circuit = stim.Circuit()

        # Add Inital Split Round
        return_circuit += self._init_split_circuit()
        return_circuit += self._init_detectors()

        # Add Repetition Rounds
        return_circuit += self._repeat_split_circuit()
        return_circuit += self._repeat_detectors()

        # Add Final Split Round
        return_circuit += self._final_split_circuit()
        return_circuit += self._final_detectors()

        return return_circuit

    def _init_detectors(self) -> stim.Circuit:
        # Define Detctor Circuit
        detector_circuit = stim.Circuit()

        # Adding Detectors for Ancilla Measurements
        # and Tracking Measurements
        self.tracker.add_measurements_to_tracker(
            measured_qubits=self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx,
            patch_type=f"Ancilla_Split_{self.split_type}",
        )
        detector_pairings = self.tracker.get_records_for_detectors(
            patch_type=f"Ancilla_Split_{self.split_type}",
        )
        for curr_pairing in detector_pairings:
            detector_circuit.append("DETECTOR", curr_pairing)

        # Adding Shift Coords
        detector_circuit.append("SHIFT_COORDS")

        # Adding Detectors for Control and Target Stabilizer Measurements
        self.tracker.add_measurements_to_tracker(
            measured_qubits=self.geometry.control_target_all_stab_idx,
            patch_type=f"Control_&_Target_Split_{self.split_type}",
        )
        detector_pairings = self.tracker.get_records_for_detectors(
            patch_type=f"Control_&_Target_Split_{self.split_type}",
        )
        for curr_pairing in detector_pairings:
            detector_circuit.append("DETECTOR", curr_pairing)

        # Adding Shift Coords
        detector_circuit.append("SHIFT_COORDS")

        return detector_circuit

    def _repeat_detectors(self) -> stim.Circuit:
        # Define Detctor Circuit
        detector_circuit = stim.Circuit()

        for _ in range(self.geometry.distance - 2):
            # Adding Detectors for Ancilla Measurements
            self.tracker.add_measurements_to_tracker(
                measured_qubits=self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx,
                patch_type=f"Ancilla_Split_{self.split_type}",
            )

            detector_pairings = self.tracker.get_records_for_detectors(
                patch_type=f"Ancilla_Split_{self.split_type}",
            )
            for curr_pairing in detector_pairings:
                detector_circuit.append("DETECTOR", curr_pairing)

            # Adding Shift Coords
            detector_circuit.append("SHIFT_COORDS")

            # Adding Detectors for Control and Target Stabilizer Measurements
            self.tracker.add_measurements_to_tracker(
                measured_qubits=self.geometry.control_target_all_stab_idx,
                patch_type=f"Control_&_Target_Split_{self.split_type}",
            )

            detector_pairings = self.tracker.get_records_for_detectors(
                patch_type=f"Control_&_Target_Split_{self.split_type}",
            )
            for curr_pairing in detector_pairings:
                detector_circuit.append("DETECTOR", curr_pairing)

            # Adding Shift Coords
            detector_circuit.append("SHIFT_COORDS")

        return detector_circuit

    def _final_detectors(self) -> stim.Circuit:
        # Define Detctor Circuit
        detector_circuit = stim.Circuit()

        # Tracking Measurements and
        # Adding Detectors for Ancilla Measurements
        self.tracker.add_measurements_to_tracker(
            measured_qubits=self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx,
            patch_type=f"Ancilla_Split_{self.split_type}",
        )
        detector_pairings = self.tracker.get_records_for_detectors(
            patch_type=f"Ancilla_Split_{self.split_type}",
        )
        for curr_pairing in detector_pairings:
            detector_circuit.append("DETECTOR", curr_pairing)

        # Adding Shift Coords
        detector_circuit.append("SHIFT_COORDS")

        # Tracking Measurements and
        # Adding Detectors for Control and Target Stabilizer Measurements
        self.tracker.add_measurements_to_tracker(
            measured_qubits=self.geometry.control_target_all_stab_idx,
            patch_type=f"Control_&_Target_Split_{self.split_type}",
        )
        detector_pairings = self.tracker.get_records_for_detectors(
            patch_type=f"Control_&_Target_Split_{self.split_type}",
        )
        for curr_pairing in detector_pairings:
            detector_circuit.append("DETECTOR", curr_pairing)

        # Adding Shift Coords
        detector_circuit.append("SHIFT_COORDS")

        return detector_circuit

    def _init_split_circuit(self) -> stim.Circuit:
        # Define initial split Circuit
        split_init_circuit = stim.Circuit()

        # Adding resets from merge ac & at
        split_init_circuit.append("TICK")
        split_init_circuit.append(
            "R",
            self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
        )

        split_init_circuit.append("TICK")
        split_init_circuit.append("H", self.geometry.combined_x_stab_idx_filtered)

        # CX Operations for Ancilla
        split_init_circuit.append("TICK")
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_init_circuit,
        )

        # Retreive Boundary + Normal Stabilizers Ancilla
        # I.e. Basis switch and measurement of ancillas
        split_init_circuit.append("H", self.geometry.anc_x_stb_idx)
        split_init_circuit.append("TICK")

        split_init_circuit.append("M", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_init_circuit.append("TICK")

        # Adding Reset and basis preparation for Ancilla for next round
        split_init_circuit.append("R", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_init_circuit.append("TICK")

        split_init_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
        split_init_circuit.append("TICK")

        # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_init_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Retreive Boundary + Normal Stabilizers from Target and Control
        # (Basis Change + Measurement):
        split_init_circuit.append(
            "H",
            self.geometry.control_x_stb_idx + self.geometry.target_x_stb_idx,
        )
        split_init_circuit.append("TICK")

        split_init_circuit.append("M", self.geometry.control_target_all_stab_idx)

        return split_init_circuit

    def _repeat_split_circuit(self) -> stim.Circuit:
        # Implementing Repeat Block
        split_repeat_circuit = stim.Circuit()

        # Reset Ancilla and prepare measurement basis
        split_repeat_circuit.append("TICK")
        split_repeat_circuit.append("R", self.geometry.control_target_all_stab_idx)

        split_repeat_circuit.append("TICK")
        split_repeat_circuit.append("H", self.geometry.combined_x_stab_idx_filtered)

        split_repeat_circuit.append("TICK")

        # CX Operations for Ancilla
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_repeat_circuit,
        )

        # Retreive Boundary + Normal Stabilizers Ancilla:
        # I.e. Basis switch and measurement of ancillas
        split_repeat_circuit.append("H", self.geometry.anc_x_stb_idx)
        split_repeat_circuit.append("TICK")

        split_repeat_circuit.append("M", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_repeat_circuit.append("TICK")

        # Adding Reset and basis preparation for Ancilla for next round
        split_repeat_circuit.append("R", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_repeat_circuit.append("TICK")

        split_repeat_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
        split_repeat_circuit.append("TICK")

        # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_repeat_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Retreive Boundary + Normal Stabilizers from Target and Control
        # (Basis Change + Measurement):
        split_repeat_circuit.append(
            "H",
            self.geometry.control_x_stb_idx + self.geometry.target_x_stb_idx,
        )
        split_repeat_circuit.append("TICK")
        split_repeat_circuit.append("M", self.geometry.control_target_all_stab_idx)

        return split_repeat_circuit * (self.geometry.distance - 2)

    def _final_split_circuit(self) -> stim.Circuit:
        # Adding Final Circ
        """
        In this Section we add the Conditional X_L and Z_L depending on the XX and ZZ Measurements
        """

        split_final_circuit = stim.Circuit()

        # Adding Reset and basis preparation
        split_final_circuit.append("TICK")
        split_final_circuit.append("R", self.geometry.control_target_all_stab_idx)

        split_final_circuit.append("TICK")
        split_final_circuit.append("H", self.geometry.combined_x_stab_idx_filtered)

        split_final_circuit.append("TICK")

        # CX Operations for Ancilla
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_final_circuit,
        )

        # Retreive Boundary + Normal Stabilizers Ancilla:
        # I.e. Basis switch and measurement of ancillas
        split_final_circuit.append("H", self.geometry.anc_x_stb_idx)
        split_final_circuit.append("TICK")

        split_final_circuit.append("M", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_final_circuit.append("TICK")

        split_final_circuit.append("R", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_final_circuit.append("TICK")
        split_final_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
        split_final_circuit.append("TICK")

        # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_final_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Retreive Boundary + Normal Stabilizers from Target and Control
        # (Basis Change + Measurement):
        split_final_circuit.append(
            "H",
            self.geometry.control_x_stb_idx + self.geometry.target_x_stb_idx,
        )
        split_final_circuit.append("TICK")
        split_final_circuit.append("M", self.geometry.control_target_all_stab_idx)

        # Adding Conditional Operations depending on non-deterministic measurements
        # Corrections appliead after all splits/ merges i.e. in AT split
        if self.split_type in {"AT"}:
            split_final_circuit.append("TICK")
            split_final_circuit += self._get_conditional_operations()

        return split_final_circuit

    def _get_conditional_operations(self):
        """
        We need logical X or Z corrections depending on the non determinstic measurement outcome of
        the newly introduzed stabilizers on the merge (i.e. the old boundary stabilizers)

        -> Use measurement tracker to find the correct record targets which where tagged
        """

        # init return Circuit
        conditional_operations_circuit = stim.Circuit()

        ####################################
        # Adding Conditional CZ/CX-Operators
        ####################################

        # Getting Measurement recs from tracker
        logical_obs_rec_tar_ac = self.tracker.get_tagged_measurements(
            tag="AC_non_deterministic_measurements",
        )

        logical_obs_rec_tar_at = self.tracker.get_tagged_measurements(
            tag="AT_non_deterministic_measurements",
        )

        for records in logical_obs_rec_tar_at:
            for data in self.log_strings["t_z"]:
                conditional_operations_circuit.append("CZ", [stim.target_rec(records + 8), data])

        ###################################
        # Measuring Ancilla in the Z Basis
        ###################################

        # Measuring Data
        conditional_operations_circuit.append("TICK")
        conditional_operations_circuit.append("MZ", self.geometry.anc_data_idx)
        conditional_operations_circuit.append("TICK")

        # Adding the conditional Gate on Control (XORing two measurements)
        # 1) Z measurements on data Ancilla
        for rec_tar, index in enumerate(self.geometry.anc_data_idx):
            if index in self.log_strings["a_z"]:
                for data in self.log_strings["c_x"]:
                    conditional_operations_circuit.append(
                        "CX",
                        [stim.target_rec(-len(self.geometry.anc_data_idx) + rec_tar), data],
                    )

        # 2) XOR from AC non deterministic ZZ measurement
        for records in logical_obs_rec_tar_ac:
            for data in self.log_strings["c_x"]:
                conditional_operations_circuit.append(
                    "CX",
                    [stim.target_rec(records + 8), data],
                )

        return conditional_operations_circuit
