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

# Translation Layer needed between Surgery and Surface Classes
STATE_TRANSLATOR = {"Y+": "+i", "Y-": "-i"}

class SurgeryReset:
    def __init__(self, geometry: SurgeryGeometry):
        # Perliminary Setup
        self.geometry = geometry
        self.control_state_init = geometry.control_state_init
        self.target_state_init = geometry.target_state_init

        # Y Basis Fault Tolerant Setting
        # Currently only non-fault tolerant Y Basis fully functional
        # Fault tolerant every flow valid except YY -> XZ!
        self.fault_tolerant_y = False

    def build_circuit(self) -> stim.Circuit:
        return_circuit = stim.Circuit()

        # Adding Coords
        return_circuit += self._adding_coords()

        # Adding Resets
        return_circuit += self._adding_resets(patch_type="control")
        return_circuit += self._adding_resets(patch_type="target")

        return return_circuit

    def _adding_coords(self) -> stim.Circuit:
        coord_circuit = stim.Circuit()

        # Appending Coords
        for q, i in self.geometry.q2i.items():
            coord_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

        return coord_circuit

    def _adding_resets(self, patch_type: str) -> stim.Circuit:
        # Init Reset Circuit
        reset_circuit = stim.Circuit()
        reset_circuit.append("TICK")

        if patch_type == "control":
            initial_state = self.control_state_init
        else:
            initial_state = self.target_state_init

        # 2. Control Patch
        if initial_state in {"Y+", "Y-"}:
            return_circuit = self._y_patch_builder(
                patch_type=patch_type,
                state_init=initial_state,
                fault_tolerant=self.fault_tolerant_y,
            )
            reset_circuit += return_circuit
            reset_circuit += self._add_y_state_flip(patch_type)
        else:
            reset_circuit += self._build_std_patch(
                patch_type=patch_type,
                state_init=initial_state,
            )

        return reset_circuit

    def _build_std_patch(self, patch_type: str, state_init: str) -> stim.Circuit:
        circ = stim.Circuit()
        log_strings = self.geometry.get_logical_strings()

        if patch_type == "control":
            # Control Patch
            data_idx = self.geometry.control_data_idx
            full_stabs = self.geometry.control_x_stb_idx + self.geometry.control_z_stb_idx
            log_z = log_strings["c_z"]
            log_x = log_strings["c_x"]
        else:
            # Target Patch
            data_idx = self.geometry.target_data_idx
            full_stabs = self.geometry.target_x_stb_idx + self.geometry.target_z_stb_idx
            log_z = log_strings["t_z"]
            log_x = log_strings["t_x"]

        # 1. Physical Qubit Reset
        if state_init in {"X+", "X-"}:
            circ.append("RX", data_idx)
        # Reset in Z Basis for I and Z init State
        else:
            circ.append("R", data_idx)

        # 2. Ancilla Qubit Reset -> Always Z Basis
        circ.append("RZ", full_stabs)

        # 2. Logical Operators for state prep
        if state_init == "X-":
            circ.append("Z", log_z)
        elif state_init in {"Z1", "I1"}:
           circ.append("X", log_x)

        return circ

    def _y_patch_builder(
        self,
        patch_type: str,
        state_init: str,
        fault_tolerant: bool,
    ) -> stim.Circuit:
        """
        Builds a Y-Basis Patch Circuit

        - Fault Tolerant True:
            Full Inplace Access of Y Basis Patch proposed by Gidney (i.e. Y-Switch etc.)
        - Fault Tolerant False:
            Simple Reset in the exact Basis where the paulis of the logical operator lie.
        """

        if fault_tolerant:
            return self._ft_y_patch(patch_type, STATE_TRANSLATOR[state_init])
        else:
            return self._non_ft_y_patch(patch_type)

    def _non_ft_y_patch(
        self,
        patch_type: str,
    ) -> stim.Circuit:

        # Get logical string
        log_strings = self.geometry.get_logical_strings()

        # Initialize Measurement Circuit
        reset_circuit = stim.Circuit()

        # Get Patch specific Info
        if patch_type == "control":
            anc_qubits = self.geometry.control_x_stb_idx + self.geometry.control_z_stb_idx
            key_prefix = "c"
        elif patch_type == "target":
            anc_qubits = self.geometry.target_x_stb_idx + self.geometry.target_z_stb_idx
            key_prefix = "t"
        else:
            raise ValueError("patch_type must be one of: 'control', 'target'")

        y_info = log_strings[f"{key_prefix}_y"]
        y_corner = y_info["y_corner"]
        x_string = y_info["x_string"]
        z_string = y_info["z_string"]

        # Add Measurements
        reset_circuit.append("RZ", anc_qubits)
        reset_circuit.append("RY", y_corner)
        reset_circuit.append("RX", x_string)
        reset_circuit.append("RZ", z_string)

        return reset_circuit

    def _add_y_state_flip(self, patch_type:str):

        # Get logical string
        log_strings = self.geometry.get_logical_strings()

        # Init Circuit
        y_flip_circuit = stim.Circuit()

        if patch_type == "control":
            flipped_x = log_strings["c_x"]
            state_init = self.geometry.control_state_init
        else:
            flipped_x = log_strings["t_x"]
            state_init = self.geometry.target_state_init

        # Add Logical Flip if needed
        if state_init == "Y-":
            y_flip_circuit.append("X", flipped_x)

        return y_flip_circuit

    def _ft_y_patch(self, patch_type: str, state_init: str) -> stim.Circuit:
        """
        Implements the Fault Tolerant Y Basis Patch Reset + Initialization + Y-Memory + Y-Switch
        as proposed by Gidney. (arXiv:2302.07395)
        """

        if patch_type == "control":
            # Control Patch Offset
            offset = 0 + (self.geometry.distance * 2) * 1j
        else:
            # Target patch Offset
            offset = (self.geometry.distance * 2) + 0j

        # Create Y-basis Geometry for the patch
        # Logical Observable can be arbitrary (I think)
        patch_geom_y = SurfaceGeometry(
            distance=self.geometry.distance,
            state_init=state_init,
            logical_observable="Y",
            y_basis=True,
            offset=offset,
            q2i=self.geometry.q2i,
        )

        # OVERRIDE q2i with the global surgery q2i
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

        # Initalizing Circuits
        patch_circuit = stim.Circuit()
        flow_circuit = stim.Circuit()
        return_circuit = stim.Circuit()

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
        flow_circuit += repetition_circuit
        flow_circuit += switch_circuit

        # Adding Observables
        observable_circuit = self._getting_y_observable(flow_circuit)

        # Building return Circuit
        return_circuit += patch_circuit + observable_circuit

        return return_circuit

    def _getting_y_observable(self, circuit: stim.Circuit) -> stim.Circuit:
        # Init return circuit
        observable_circuit = stim.Circuit()

        # Getting Logical Strings
        log_strings = self.geometry.get_logical_strings()

        ######################## TO-DO #############################
        # -> Check documentation for better handling of these flows
        # -> This is a pure eyesore
        ############################################################

        # 1) Control Flow:
        if self.control_state_init in {"Y+", "Y-"}:
            # Getting Logical Y Pauli Strings for Control
            y_corner = log_strings["c_y"]["y_corner"][0]
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

            # Try to map to negative recs if not, unmodified circuit gets returned
            try:
                for index_creation in logical_creation_rec:
                    current_rec_crea = circuit.num_measurements - index_creation
                    rec_pos.append(-current_rec_crea)
            except Exception:
                print("No valid flow found control, skipping observable inclusion.")

            # Adding Measurement Rec Correction given by Flow
            observable_circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_rec(k) for k in rec_pos],
                0,
            )

        # 2) Target Flow:
        if self.target_state_init in {"Y+", "Y-"}:
            # Getting Logical Y Pauli Strings for Target
            y_corner = log_strings["t_y"]["y_corner"][0]
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

            # Try to map to negative recs if not, unmodified circuit gets returned
            try:
                for index_creation in logical_creation_rec:
                    current_rec_crea = circuit.num_measurements - index_creation
                    rec_pos.append(-current_rec_crea)
            except Exception:
                print("No valid flow found for target, skipping observable inclusion.")

            # Adding Measurement Rec Correction given by Flow
            observable_circuit.append(
                "OBSERVABLE_INCLUDE",
                [stim.target_rec(k) for k in rec_pos],
                0,
            )

        return observable_circuit

    def _debug_print_available_flows(self, flow_circuit: stim.Circuit, must_have: list[str] = None):
        """
        Returns all available flows or if must_have is specified only those containing all the
        strings listed
        """

        available_flows = flow_circuit.flow_generators()
        for flow in available_flows:
            if must_have is not None:
                if all(item in str(flow) for item in must_have):
                    print(flow)
            else:
                print(flow)
