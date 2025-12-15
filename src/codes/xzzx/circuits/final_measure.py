import stim

from src.codes.xzzx.xzzx_geom import XZZXGeometry

__all__ = ["FinalMeasureCircuit"]


class FinalMeasureCircuit:
    def __init__(self, geometry: XZZXGeometry):
        self.geometry = geometry

    def build_final_measurement_circuit(self) -> stim.Circuit:
        self.final_circuit = stim.Circuit()

        # Building Circuit
        self.final_circuit += self._apply_final_measurement()
        self.final_circuit += self._apply_detectors()
        self.final_circuit += self._apply_observables()

        return self.final_circuit

    def _apply_final_measurement(self) -> stim.Circuit:
        ##########################
        # Adding final measurement
        ##########################

        final_circuit = stim.Circuit()

        final_circuit.append("SHIFT_COORDS", arg=(0, 0, 1))

        ###########################################
        # 6) Implementing Final Measurement-Round
        # -> Detectors compromised by last round stab. measurements
        #    & Data parity checks from final measurement
        ###########################################

        final_circuit.append("MZ", self.geometry.data_z_idx)
        final_circuit.append("MX", self.geometry.data_x_idx)

        return final_circuit

    def _apply_detectors(self) -> stim.Circuit:
        # Defining Det Circuit
        det_circuit = stim.Circuit()

        # -> Defining Data to measurement indexing
        index_to_rec_data: dict[int, int] = {
            q: i
            for i, q in enumerate(reversed(self.geometry.data_z_idx + self.geometry.data_x_idx))
        }
        index_to_rec_ancilla: dict[int, int] = {
            q: i for i, q in enumerate(reversed(self.geometry.stab_idx))
        }

        """
        Why do we only check the Z stabilizers in the last measurement round and 
        not also the x stabilizers as we did before
        -> We need to meassure the X stabilizers in the X basis but we already 
        meassured in the z basis (XZ do not commute)
        """

        for q, qtype in self.geometry.qubit_coords.items():
            if (
                self.geometry.state_init == "XZZX-VER"
                and qtype in {"STAB-VER", "STAB-BOUND-A-Ver", "STAB-BOUND-B-Ver"}
                or (
                    self.geometry.state_init == "XZZX-HOR"
                    and qtype in {"STAB-HOR", "STAB-BOUND-L-Hor", "STAB-BOUND-R-Hor"}
                )
            ):
                # Getting neighbor qubits of the ancilla
                neighbors = self.geometry.get_neighbors(q, qtype)

                # Defining the current record targets
                current_record = [-index_to_rec_data[neighbor] - 1 for neighbor in neighbors]

                # Defining the last record targets (Normal Detectors from last round)
                ancilla_index = self.geometry.q2i[q]
                last_record = [
                    -index_to_rec_ancilla[ancilla_index]
                    - 1
                    - len(self.geometry.data_z_idx + self.geometry.data_x_idx),
                ]

                # Combining the record targets
                final_record = current_record + last_record

                # Appending Detector
                det_circuit.append(
                    "DETECTOR",
                    [stim.target_rec(i) for i in final_record],
                    arg=(q.real, q.imag, 1),
                )

        return det_circuit

    def _apply_observables(self) -> stim.Circuit:
        observable_circuit = stim.Circuit()

        ########################################################
        # 7) Defining logical Operators (Final Measurement-Round)
        ########################################################

        if self.geometry.state_init in {"XZZX-VER"}:
            # Control stabilized by x logical
            tar_rec = []

            for rec_pos, index in enumerate(self.geometry.data_z_idx + self.geometry.data_x_idx):
                if index in self.geometry.get_logical_indexes("VER"):
                    tar_rec.append(rec_pos)

            observable_circuit.append(
                "OBSERVABLE_INCLUDE",
                [
                    stim.target_rec(-len(self.geometry.data_z_idx + self.geometry.data_x_idx) + k)
                    for k in tar_rec
                ],
                0,
            )

        elif self.geometry.state_init in {"XZZX-HOR"}:
            # Control stabilized by x logical
            tar_rec = []

            for rec_pos, index in enumerate(self.geometry.data_z_idx + self.geometry.data_x_idx):
                if index in self.geometry.get_logical_indexes("HOR"):
                    tar_rec.append(rec_pos)

            observable_circuit.append(
                "OBSERVABLE_INCLUDE",
                [
                    stim.target_rec(-len(self.geometry.data_z_idx + self.geometry.data_x_idx) + k)
                    for k in tar_rec
                ],
                0,
            )

        return observable_circuit
