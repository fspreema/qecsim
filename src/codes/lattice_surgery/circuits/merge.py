import stim

from src.codes.lattice_surgery.data_geometry import MasterPairings
from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.core.cx_builder import cx_builder
from src.core.measurement_tracker import MeasurementTracker

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
            self.non_det_stab_indices = self.geometry.non_det_stab_indices_ac

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
            self.non_det_stab_indices = self.geometry.non_det_stab_indices_at

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

    def _get_detectors(
        self,
        measured_qubits: list[int],
        patch_type: str,
        qubits_for_detectors: list[int] | None = None,
        tag: str | None = None,
        tagged_qubits: list[int] | None = None,
    ) -> stim.Circuit:
        # Initialize Circuit
        detector_circuit = stim.Circuit()

        self.tracker.add_measurements_to_tracker(
            measured_qubits=measured_qubits,
            qubits_for_detectors=qubits_for_detectors,
            patch_type=patch_type,
            tag=tag,
            tagged_qubits=tagged_qubits,
        )

        det_record_pairings = self.tracker.get_records_for_detectors(
            patch_type=patch_type,
        )
        for curr_pairing in det_record_pairings:
            detector_circuit.append("DETECTOR", curr_pairing)

        # Shifting Coords
        detector_circuit.append("SHIFT_COORDS")

        return detector_circuit

    def _initial_merge_circuit(self) -> stim.Circuit:
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
        merge_init_circuit.append("TICK")

        # Adding Detectors for Combined Patches
        merge_init_circuit += self._get_detectors(
            measured_qubits=self.combined_z_stab_merging_lattices
            + self.combined_x_stab_merging_lattices,
            patch_type=f"Merge_{self.merging_type}",
        )

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

        # Adding Detectors for untouched Patches
        merge_init_circuit += self._get_detectors(
            measured_qubits=self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
            patch_type=f"Merge_{self.merging_type}_untouched",
        )

        return merge_init_circuit

    def _repeat_merge_circuit(self) -> stim.Circuit:
        # Defining Repeat Circuit
        merge_round_circuit = stim.Circuit()


        ### CHECK IF RANGE NOT OFF BY ONE AKA DISTANCE - 1
        for curr_round in range(self.geometry.distance):
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

            # Adding Detectors for Combined Patches & Labeling Non Deterministic Measurements
            merge_round_circuit += self._get_detectors(
                measured_qubits=self.combined_z_stab_merging_lattices
                + self.combined_x_stab_merging_lattices,
                patch_type=f"Merge_{self.merging_type}",
                tagged_qubits=self.non_det_stab_indices,
                tag=f"{self.merging_type}_non_deterministic_measurements"
                if curr_round == self.geometry.distance - 1
                else None,
            )

            # Adding Reset for Ancilla Qubits of Stabilizers
            merge_round_circuit.append(
                "R",
                self.combined_z_stab_merging_lattices + self.combined_x_stab_merging_lattices,
            )
            merge_round_circuit.append("TICK")

            # Continue CX-Implementation for untouched lattice
            # (As AC/AT-Lattice already has a full run)

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

            # Adding Detectors for untouched Patches
            merge_round_circuit += self._get_detectors(
                measured_qubits=self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
                patch_type=f"Merge_{self.merging_type}_untouched",
            )

        return merge_round_circuit
