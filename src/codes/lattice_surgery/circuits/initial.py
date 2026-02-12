import stim

from src.codes.lattice_surgery.data_geometry import MasterPairings
from src.codes.lattice_surgery.measurement_tracker import MeasurementTracker
from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.core.cx_builder import cx_builder

Coord = complex

__all__ = ["SurgeryInitialization"]


class SurgeryInitialization:
    def __init__(
        self,
        geometry: SurgeryGeometry,
        master_pairings: MasterPairings,
        measurement_tracker: MeasurementTracker,
    ):
        # Preliminary Setup
        self.geometry = geometry
        self.stab_to_data = master_pairings.std_pairings.get_schedule()
        self.measurement_tracker = measurement_tracker

    def build_circuit(self) -> stim.Circuit:
        circuit = stim.Circuit()

        # Adding Initializations
        circuit += self._adding_ancilla_initializations()
        circuit += self._adding_control_target_initializations()

        # Adding Init Repeat Block
        circuit += self._adding_repeat_block()

        return circuit

    def _adding_ancilla_initializations(self):
        # Initialization Circuit
        anc_init_circuit = stim.Circuit()

        # Resetting Ancilla Data Qubits in X Basis and Stabilizers in Z Basis
        # -> This needs to be in the inital circuit as this is the reset
        #    that needs to be inside the flow_circuit! (At least the Data reset)
        anc_init_circuit.append("TICK")
        anc_init_circuit.append("RZ", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        anc_init_circuit.append("RX", self.geometry.anc_data_idx)
        anc_init_circuit.append("TICK")

        # Adding h gate for X stabilizers -> Filtering out double coords
        anc_init_circuit.append("H", self.geometry.combined_x_stab_idx_filtered)
        anc_init_circuit.append("TICK")

        # CX Operations for Ancilla qubits
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_to_data,
            circuit=anc_init_circuit,
        )

        # Basis Change and Measurement of Ancilla Stabilizers
        anc_init_circuit.append("TICK")
        anc_init_circuit.append("H", self.geometry.anc_x_stb_idx)
        anc_init_circuit.append("TICK")
        anc_init_circuit.append("M", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        anc_init_circuit.append("TICK")

        # Adding Detectors -> Only X Type Stabilizers for Ancilla Init
        # -> Init. in X Basis
        self.measurement_tracker.add_measurements_to_detector_dict(
            measured_qubits=self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx,
            selected_qubits=self.geometry.anc_x_stb_idx,
            patch_type="Ancilla",
        )

        det_record_pairings = self.measurement_tracker.get_records_for_detectors(
            patch_type="Ancilla",
        )
        for curr_pairing in det_record_pairings:
            anc_init_circuit.append("DETECTOR", curr_pairing)

        # Shifting Coords
        anc_init_circuit.append("SHIFT_COORDS")

        return anc_init_circuit

    def _adding_control_target_initializations(self):
        # Initialization Circuit
        ct_init_circuit = stim.Circuit()

        # Resetting Boundary Stabilizers which are shared with Target/ Control
        ct_init_circuit.append("R", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        ct_init_circuit.append("TICK")
        ct_init_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
        ct_init_circuit.append("TICK")

        # CX Operations for Control and Target qubits
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_to_data,
            circuit=ct_init_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Basis Change and Measurement of Control and Target Stabilizers
        ct_init_circuit.append("TICK")
        ct_init_circuit.append(
            "H",
            self.geometry.control_x_stb_idx + self.geometry.target_x_stb_idx,
        )
        ct_init_circuit.append("TICK")
        ct_init_circuit.append("M", self.geometry.control_target_all_stab_idx)

        # Adding Detectors -> Selected qubits depend on Basis
        control_stabs = (
            self.geometry.control_x_stb_idx
            if self.geometry.control_state_init in {"X+", "X-"}
            else self.geometry.control_z_stb_idx
        )
        target_stabs = (
            self.geometry.target_x_stb_idx
            if self.geometry.target_state_init in {"X+", "X-"}
            else self.geometry.target_z_stb_idx
        )

        self.measurement_tracker.add_measurements_to_detector_dict(
            measured_qubits=self.geometry.control_target_all_stab_idx,
            selected_qubits=control_stabs + target_stabs,
            patch_type="Control_&_Target",
        )

        det_record_pairings = self.measurement_tracker.get_records_for_detectors(
            patch_type="Control_&_Target",
        )
        for curr_pairing in det_record_pairings:
            ct_init_circuit.append("DETECTOR", curr_pairing)

        # Shifting Coords
        ct_init_circuit.append("SHIFT_COORDS")

        return ct_init_circuit

    def _adding_repeat_block(self):
        rep_init_circuit = stim.Circuit()

        # Adding reset from initial round
        rep_init_circuit.append("TICK")
        rep_init_circuit.append("R", self.geometry.control_target_all_stab_idx)

        rep_init_circuit.append("TICK")
        rep_init_circuit.append("H", self.geometry.combined_x_stab_idx_filtered)
        rep_init_circuit.append("TICK")

        # CX Operations for Ancilla qubits
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_to_data,
            circuit=rep_init_circuit,
        )

        # Retreive Boundary + Normal Stabilizers Ancilla
        # I.e. Basis change and measurement of ancilla
        rep_init_circuit.append("H", self.geometry.anc_x_stb_idx)
        rep_init_circuit.append("TICK")

        rep_init_circuit.append(
            "M",
            self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx,
        )
        rep_init_circuit.append("TICK")

        # Adding Full Ancilla Detectors
        self.measurement_tracker.add_measurements_to_detector_dict(
            measured_qubits=self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx,
            patch_type="Ancilla",
        )

        det_record_pairings = self.measurement_tracker.get_records_for_detectors(
            patch_type="Ancilla",
        )

        for curr_pairing in det_record_pairings:
            rep_init_circuit.append("DETECTOR", curr_pairing)

        # Shifting Coords
        rep_init_circuit.append("SHIFT_COORDS")

        rep_init_circuit.append(
            "R",
            self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx,
        )
        rep_init_circuit.append("TICK")
        rep_init_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
        rep_init_circuit.append("TICK")

        # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.stab_to_data,
            circuit=rep_init_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Retreive Boundary + Normal Stabilizers from Target and Control
        # (Basis Change + Measurement):
        rep_init_circuit.append(
            "H",
            self.geometry.control_x_stb_idx + self.geometry.target_x_stb_idx,
        )
        rep_init_circuit.append("TICK")
        rep_init_circuit.append("M", self.geometry.control_target_all_stab_idx)

        # Adding Full Control and Target Detectors
        self.measurement_tracker.add_measurements_to_detector_dict(
            measured_qubits=self.geometry.control_target_all_stab_idx,
            patch_type="Control_&_Target",
        )

        det_record_pairings = self.measurement_tracker.get_records_for_detectors(
            patch_type="Control_&_Target",
        )

        for curr_pairing in det_record_pairings:
            rep_init_circuit.append("DETECTOR", curr_pairing)

        # Shifting Coords
        rep_init_circuit.append("SHIFT_COORDS")

        return rep_init_circuit * (self.geometry.distance - 1)
