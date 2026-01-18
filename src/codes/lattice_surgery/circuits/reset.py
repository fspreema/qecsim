import stim

from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.codes.surface_code_rotated.circuits.initial import SurfaceInitialization
from src.codes.surface_code_rotated.circuits.repetition import SurfaceRepetitionCircuit
from src.codes.surface_code_rotated.circuits.reset import SurfaceReset
from src.codes.surface_code_rotated.circuits.y_switch import YSwitchCircuit
from src.codes.surface_code_rotated.data_geometry import (
    MasterGeometry,
    MasterPairings,
)
from src.codes.surface_code_rotated.get_stab_pairings import SurfacePairings
from src.codes.surface_code_rotated.surface_geom import SurfaceGeometry

__all__ = ["SurgeryReset"]


class SurgeryReset:
    def __init__(self, geometry: SurgeryGeometry):
        # Perliminary Setup
        self.geometry = geometry
        self.control_state_init = geometry.control_state_init
        self.target_state_init = geometry.target_state_init

        # Getting Y-Basis Stab to Data Information

    def build_circuit(self) -> stim.Circuit:
        circuit = stim.Circuit()

        # Adding Coords
        circuit += self._adding_coords()

        # Adding Resets
        reset_circuit, flow_circuit = self._adding_resets()
        circuit += reset_circuit

        # Adding Y-Creation Flows if relevant
        if self.control_state_init in {"+i", "-i"} or self.target_state_init in {"+i", "-i"}:
            circuit += self._getting_y_observable(flow_circuit)

        return circuit

    def _adding_coords(self):
        coord_circuit = stim.Circuit()

        # Appending Coords
        for q, i in self.geometry.q2i.items():
            coord_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

        return coord_circuit

    def _adding_resets(self):
        # Init Reset Circuit
        reset_circuit = stim.Circuit()
        flow_circuit = stim.Circuit()

        # 1. Reset Ancilla Data Qubits into X Basis and Stabs in Z Basis of Ancilla Patch
        reset_circuit.append("RX", self.geometry.anc_data_idx)

        # 2. Control Patch
        if self.control_state_init in {"+i", "-i"}:
            return_circuit, flow_creation_circuit = self._build_y_patch(
                patch_type="control",
                state_init=self.control_state_init,
            )
            reset_circuit += return_circuit
            flow_circuit += flow_creation_circuit
        else:
            reset_circuit += self._build_std_patch(
                patch_type="control",
                state_init=self.control_state_init,
            )

        # 3. Target Patch
        if self.target_state_init in {"+i", "-i"}:
            return_circuit, flow_creation_circuit = self._build_y_patch(
                patch_type="target",
                state_init=self.target_state_init,
            )
            reset_circuit += return_circuit
            flow_circuit += flow_creation_circuit
        else:
            reset_circuit += self._build_std_patch(
                patch_type="target",
                state_init=self.target_state_init,
            )

        reset_circuit.append("TICK")

        return reset_circuit, flow_circuit

    def _build_std_patch(self, patch_type: str, state_init: str) -> stim.Circuit:
        circ = stim.Circuit()
        log_strings = self.geometry.get_logical_strings()

        if patch_type == "control":
            data_idx = self.geometry.control_data_idx
            log_z = log_strings["c_z"]
            log_x = log_strings["c_x"]
        else:  # target
            data_idx = self.geometry.target_data_idx
            log_z = log_strings["t_z"]
            log_x = log_strings["t_x"]

        # 1. Physical Reset
        if state_init in {"+", "-"}:
            circ.append("RX", data_idx)
        else:
            circ.append("R", data_idx)

        # 2. Logical Operators for state prep
        if state_init == "-":
            circ.append("Z", log_z)
        elif state_init == "1":
            circ.append("X", log_x)

        return circ

    def _build_y_patch(self, patch_type: str, state_init: str) -> stim.Circuit:
        if patch_type == "control":
            offset = 0 + (self.geometry.distance * 2) * 1j
        else:  # target
            offset = (self.geometry.distance * 2) + 0j

        # Create Y-basis Geometry for the patch
        patch_geom_y = SurfaceGeometry(
            distance=self.geometry.distance,
            offset=offset,
            starting_stabilizer_x=False,
            y_basis=True,
        )
        # OVERRIDE q2i with the global surgery q2i
        patch_geom_y.q2i = self.geometry.q2i
        patch_geom_y.state_init = state_init

        # Create Y-basis Pairings
        patch_pairings_y = SurfacePairings(
            patch=patch_geom_y.coords,
            distance=self.geometry.distance,
            y_basis=True,
            offset=offset,
        )

        # Create Y-Switch Pairings
        patch_pairings_y_switch = SurfacePairings(
            patch=patch_geom_y.coords,
            distance=self.geometry.distance,
            y_basis=True,
            y_switch=True,
            offset=offset,
        )

        # Wrap in Master Objects
        master_geom = MasterGeometry(
            geometry_std=None,  # Not used for y_basis type
            geometry_ybasis=patch_geom_y,
        )
        master_pairings = MasterPairings(
            pairings_std=None,
            pairings_ybasis=patch_pairings_y,
            pairings_ymemory=None,
            pairings_yswitch=patch_pairings_y_switch,
            pairings_log_h=None,
        )

        patch_circuit = stim.Circuit()

        # 1. Data Qubits Reset
        resetter = SurfaceReset(master_geom, type="y_basis")
        patch_circuit += resetter._adding_resets()

        # 2. Stabilizer Resets + Initialization Circuit
        initer = SurfaceInitialization(master_geom, master_pairings, type="y_basis")
        patch_circuit += initer.build_circuit()

        # 3. Y-Memory Repetition
        rep = SurfaceRepetitionCircuit(master_geom, master_pairings, type="y_basis")
        repetition_circuit = rep.build_circuit()
        patch_circuit += repetition_circuit

        # 4. Y-Switch
        switch = YSwitchCircuit(master_geom, master_pairings)
        switch_circuit = switch.build_circuit()
        patch_circuit += switch_circuit

        # Defining Flow Circuit for Logical Y Creation
        flow_circuit = stim.Circuit()
        flow_circuit += repetition_circuit
        flow_circuit += switch_circuit

        return patch_circuit, flow_circuit

    def _getting_y_observable(self, circuit: stim.Circuit):
        # Init return circuit
        observable_circuit = stim.Circuit()

        # Getting Logical Strings
        log_strings = self.geometry.get_logical_strings()

        # 1) Control Flow:
        if self.control_state_init in {"+i", "-i"}:
            # Getting Logical Y Pauli Strings for Control
            y_corner = [log_strings["c_y"]["y_corner"]][0]
            logical_xyz_string = "*".join(
                [f"Z{idz}" for j, idz in enumerate(log_strings["c_y"]["z_string"])]
                + [f"Y{y_corner}"]
                + [f"X{idx}" for j, idx in enumerate(log_strings["c_y"]["x_string"])],
            )

            logical_creation = f"{1} -> {logical_xyz_string}"

            (logical_creation_rec,) = circuit.solve_flow_measurements(
                [stim.Flow(logical_creation)],
            )

            # Adding the Observable
            rec_pos = []

            for index_creation in logical_creation_rec:
                current_rec_crea = circuit.num_measurements - index_creation
                rec_pos.append(-current_rec_crea)

            observable_circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_rec(k) for k in rec_pos],
                0,
            )

        # 2) Target Flow:
        if self.target_state_init in {"+i", "-i"}:
            # Getting Logical Y Pauli Strings for Control
            y_corner = [log_strings["t_y"]["y_corner"]][0]
            logical_xyz_string = "*".join(
                [f"Z{idz}" for j, idz in enumerate(log_strings["t_y"]["z_string"])]
                + [f"Y{y_corner}"]
                + [f"X{idx}" for j, idx in enumerate(log_strings["t_y"]["x_string"])],
            )

            logical_creation = f"{1} -> {logical_xyz_string}"

            (logical_creation_rec,) = circuit.solve_flow_measurements(
                [stim.Flow(logical_creation)],
            )

            # Adding the Observable
            rec_pos = []

            for index_creation in logical_creation_rec:
                current_rec_crea = circuit.num_measurements - index_creation
                rec_pos.append(-current_rec_crea)

            observable_circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_rec(k) for k in rec_pos],
                0,
            )

        return observable_circuit
