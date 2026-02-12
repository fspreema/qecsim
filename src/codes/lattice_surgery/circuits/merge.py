import stim

from src.codes.lattice_surgery.data_geometry import MasterPairings
from src.codes.lattice_surgery.measurement_tracker import MeasurementTracker
from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.core.cx_builder import cx_builder

Coord = complex

__all__ = ["SurgeryMerge"]


class SurgeryMerge:
    def __init__(
        self,
        geometry: SurgeryGeometry,
        master_pairings: MasterPairings,
        merging_type: str,
        tracker: MeasurementTracker,
    ):
        # Preliminary Setup
        self.geometry = geometry
        self.master_pairings = master_pairings
        self.tracker = tracker

        # Getting Specific Merging Type Info
        self.merging_type = merging_type

        if merging_type == "AC":
            self.stab_to_data_curr_merg = self.master_pairings.ac_merge_pairings.get_schedule()
            self.stab_to_data_untouched_circ = (
                self.master_pairings.std_pairings.get_specific_region(
                    patch_coords=self.geometry.coords_target,
                )
            )
            self.x_stab_index_untouched_circ = self.geometry.target_x_stb_idx
            self.z_stab_index_untouched_circ = self.geometry.target_z_stb_idx
            self.combined_x_stab = self.geometry.combined_x_stab_idx_filtered
            self.combined_x_stab_merging_lattices, self.combined_z_stab_merging_lattices = (
                self.geometry.get_combined_xz_stabs_merging_lattice(merging_type="AC")
            )
            self.non_det_stab_indices = (
                self.geometry.surgery_z_l_stb_idx + self.geometry.surgery_z_m_stb_idx
            )

        elif merging_type == "AT":
            self.stab_to_data_curr_merg = self.master_pairings.at_merge_pairings.get_schedule()
            self.stab_to_data_untouched_circ = (
                self.master_pairings.std_pairings.get_specific_region(
                    patch_coords=self.geometry.coords_control,
                )
            )
            self.x_stab_index_untouched_circ = self.geometry.control_x_stb_idx
            self.z_stab_index_untouched_circ = self.geometry.control_z_stb_idx
            self.combined_x_stab = self.geometry.combined_x_stab_idx_filtered_AT
            self.combined_x_stab_merging_lattices, self.combined_z_stab_merging_lattices = (
                self.geometry.get_combined_xz_stabs_merging_lattice(merging_type="AT")
            )
            self.non_det_stab_indices = (
                self.geometry.surgery_x_b_stb_idx + self.geometry.surgery_x_m_stb_idx
            )

        else:
            raise ValueError("No valid merging Type in Function selected!")

    def build_circuit(self):
        # Init return circuit
        return_circuit = stim.Circuit()

        # Initial Merge Circuit
        return_circuit += self._initial_merge_circuit()
        # Repeat Merge Circuit
        return_circuit += self._repeat_merge_circuit()

        return return_circuit

    def _initial_merge_circuit(self):
        # Init Circuit
        merge_init_circuit = stim.Circuit()

        # Adding reset from initial round and from AC split round
        merge_init_circuit.append("TICK")
        merge_init_circuit.append("R", self.geometry.control_target_all_stab_idx)

        merge_init_circuit.append("TICK")
        merge_init_circuit.append("H", self.combined_x_stab)

        merge_init_circuit.append("TICK")

        # CX Operations
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_to_data_curr_merg | self.stab_to_data_untouched_circ,
            circuit=merge_init_circuit,
        )

        # Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement):
        # I.e. return to Z-Basis where needed and Measure
        merge_init_circuit.append("H", self.combined_x_stab_merging_lattices)
        merge_init_circuit.append("TICK")

        merge_init_circuit.append(
            "M",
            self.combined_z_stab_merging_lattices + self.combined_x_stab_merging_lattices,
        )

        # Updating Measurement Tracker
        self.tracker.add_measurements(
            measured_qubits=self.combined_z_stab_merging_lattices
            + self.combined_x_stab_merging_lattices,
        )

        merge_init_circuit.append("TICK")

        merge_init_circuit.append(
            "R",
            self.combined_z_stab_merging_lattices + self.combined_x_stab_merging_lattices,
        )
        merge_init_circuit.append("TICK")

        # Continue CX-Implementation for untouched lattice (As AC/AT-Lattice already has a full run)
        """
        Adding needed H-Gates for X-Stabs which are shared between
        merged lattice and untouched lattice
        """

        if self.merging_type == "AT":
            merge_init_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
            merge_init_circuit.append("TICK")

        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_to_data_untouched_circ,
            circuit=merge_init_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Retreive Boundary + Normal Stabilizers from Target and Control
        # (Basis Change + Measurement):
        merge_init_circuit.append("H", self.x_stab_index_untouched_circ)
        merge_init_circuit.append("TICK")

        merge_init_circuit.append(
            "M",
            self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
        )

        # Updating Measurement Tracker
        self.tracker.add_measurements(
            measured_qubits=self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
        )

        return merge_init_circuit

    def _repeat_merge_circuit(self):
        # Defining Repeat Circuit
        merge_round_circuit = stim.Circuit()

        # Reinitializing Stabilizers and add basis change where needed
        merge_round_circuit.append("TICK")
        merge_round_circuit.append(
            "R",
            self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
        )

        merge_round_circuit.append("TICK")
        merge_round_circuit.append("H", self.combined_x_stab)

        merge_round_circuit.append("TICK")

        # CX Operations for Ancilla qubits
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_to_data_curr_merg | self.stab_to_data_untouched_circ,
            circuit=merge_round_circuit,
        )

        # Retreive Boundary + Normal Stabilizers Ancilla (Basis change and Measurement):
        # I.e. return to Z-Basis where needed and Measure
        merge_round_circuit.append("H", self.combined_x_stab_merging_lattices)
        merge_round_circuit.append("TICK")

        merge_round_circuit.append(
            "M",
            self.combined_z_stab_merging_lattices + self.combined_x_stab_merging_lattices,
        )
        merge_round_circuit.append("TICK")

        # Adding Detectors for Combined Patches
        self.tracker.add_measurements_to_detector_dict(
            measured_qubits=self.combined_z_stab_merging_lattices
            + self.combined_x_stab_merging_lattices,
            patch_type=f"Merge_{self.merging_type}",
        )

        det_record_pairings = self.tracker.get_records_for_detectors(
            patch_type=f"Merge_{self.merging_type}",
        )
        for curr_pairing in det_record_pairings:
            merge_round_circuit.append("DETECTOR", curr_pairing)

        # Shifting Coords
        merge_round_circuit.append("SHIFT_COORDS")

        # Adding Reset for Ancilla Qubits of Stabilizers
        merge_round_circuit.append(
            "R",
            self.combined_z_stab_merging_lattices + self.combined_x_stab_merging_lattices,
        )
        merge_round_circuit.append("TICK")

        # Continue CX-Implementation for untouched lattice (As AC/AT-Lattice already has a full run)

        """
        Adding needed H-Gates for X-Stabs which are shared between 
        merged lattice and untouched lattice
        """

        if self.merging_type == "AT":
            merge_round_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
            merge_round_circuit.append("TICK")

        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_to_data_untouched_circ,
            circuit=merge_round_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Retreive Boundary + Normal Stabilizers from Target and Control
        # (Basis Change + Measurement):
        merge_round_circuit.append("H", self.x_stab_index_untouched_circ)
        merge_round_circuit.append("TICK")

        merge_round_circuit.append(
            "M",
            self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
        )

        # Adding Detectors for Combined Patches
        self.tracker.add_measurements_to_detector_dict(
            measured_qubits=self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
            patch_type=f"Merge_{self.merging_type}_untouched",
        )

        det_record_pairings = self.tracker.get_records_for_detectors(
            patch_type=f"Merge_{self.merging_type}_untouched",
        )
        for curr_pairing in det_record_pairings:
            merge_round_circuit.append("DETECTOR", curr_pairing)

        # Shifting Coords
        merge_round_circuit.append("SHIFT_COORDS")

        # Updating Measurement Trackers
        for curr_round in range(self.geometry.distance - 1):
            # Updating Non Deterministic Measurement Tracker
            self.tracker.add_measurements(
                measured_qubits=self.combined_z_stab_merging_lattices
                + self.combined_x_stab_merging_lattices,
                specific_qubits=self.non_det_stab_indices,
                tag=f"{self.merging_type}_non_deterministic_measurements"
                if curr_round == self.geometry.distance - 2
                else None,
            )

            # Updating Determisntic Measurement Tracker
            self.tracker.add_measurements(
                measured_qubits=self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
            )

        return merge_round_circuit * (self.geometry.distance - 1)
