import stim

from src.codes.surface_code_rotated.data_geometry import MasterGeometry, MasterPairings
from src.core.cx_builder import cx_builder

__all__ = ["SurfaceInitialization"]


class SurfaceInitialization:
    def __init__(self, master_geometry: MasterGeometry, master_pairings: MasterPairings, type: str):
        if type not in {"standard", "y_basis", "log_h"}:
            raise ValueError(
                f"Invalid type '{type}' for SurfaceInitialization. Must be 'standard'"
                f", 'y_basis' or 'log_h'.",
            )

        self.type = type

        # Initialize Geometry and Pairings depending on the type
        if type == "standard":
            self.geometry = master_geometry.geometry_std
            self.pairings = master_pairings.pairings_std
        elif type == "y_basis":
            self.geometry = master_geometry.geometry_ybasis
            self.pairings = master_pairings.pairings_ybasis
        elif type == "log_h":
            self.geometry = master_geometry.geometry_std
            self.pairings = master_pairings.pairings_log_h

    def build_circuit(self) -> stim.Circuit:
        circuit = stim.Circuit()

        # Adding Resets
        circuit += self._adding_resets()

        # Adding Transversal H if Logical H init
        if self.type == "log_h":
            circuit += self._adding_transversal_h()

        # Adding Initializations
        circuit += self._adding_initializations()

        return circuit

    def _adding_resets(self):
        # Reset Circuit
        reset_circuit = stim.Circuit()
        reset_circuit.append("TICK")
        reset_circuit.append("R", self.geometry.stab_x_idx + self.geometry.stab_z_idx)
        reset_circuit.append("TICK")

        return reset_circuit

    def _adding_transversal_h(self):
        # Transversal H Circuit
        h_circuit = stim.Circuit()

        h_circuit.append("H", self.geometry.data_idx)
        h_circuit.append("TICK")

        return h_circuit

    def _adding_initializations(self):
        # Init Circuit
        init_circuit = stim.Circuit()

        # 1) Reset/ Basis
        init_circuit.append("H", self.geometry.stab_x_idx)
        init_circuit.append("TICK")

        # 2) CX Operations
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.pairings.stab_to_data,
            circuit=init_circuit,
            excluded_index=self.geometry.y_index if self.type == "y_basis" else None,
        )

        # -------Continue-Circuit------------

        # 3) Basis/ Measurement
        init_circuit.append("H", self.geometry.stab_x_idx)
        init_circuit.append("TICK")
        init_circuit.append("M", self.geometry.stab_x_idx + self.geometry.stab_z_idx)
        init_circuit.append("SHIFT_COORDS", arg=(0, 0, 1))
        init_circuit.append("TICK")

        return init_circuit
