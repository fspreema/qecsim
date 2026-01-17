import stim

from src.codes.lattice_surgery.circuits.final_measure import SurgeryFinalMeasure
from src.codes.lattice_surgery.circuits.initial import SurgeryInitialization
from src.codes.lattice_surgery.circuits.merge import SurgeryMerge
from src.codes.lattice_surgery.circuits.reset import SurgeryReset
from src.codes.lattice_surgery.circuits.split import SurgerySplit
from src.codes.lattice_surgery.data_geometry import MasterPairings
from src.codes.lattice_surgery.get_stab_pairings import LatticeSurgeryPairings
from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.core.base_class_builder import BaseClassBuilder
from src.core.data_models import NoiseParameters

__all__ = ["SurgeryBuilder"]


class SurgeryBuilder(BaseClassBuilder):
    def __init__(
        self,
        distance: int,
        control_state_init: str,
        target_state_init: str,
        flow_observable: str,
        noise: NoiseParameters = None,
    ):
        """
        Returns the full lattice surgery circuit

        Arguments:
            -> target_state_init:
            In which basis should the target lattice be initlized?
            -> control_state_init:
            In which basis should the control lattice be initilized?
            -> flow_observable:
            What Observable should stim track?
            -> distance/rounds:
            Distance and rounds of the circuit (If not given any value distance = rounds)
            -> noise_depol_data_init:
            Depolaization noise after initlization of the data qubits with given probability
            -> noise_measure_flip:
            X_Error before measurement to simulate faulty measurement with given probability
            -> noise_after_reset:
            X_Error after reset of ancilla in order so simulate faulty reset with given probability
            -> noise_after_clifford_depol:
            Depolarization gate after each Clifford (CX and H) with given probability

        Returns:
            -> Fully implemented CX-Gate in stim.Circuit format
        """

        # Init Parameters
        self.distance = distance
        self.control_state_init = control_state_init
        self.target_state_init = target_state_init
        self.flow_observable = flow_observable
        self.noise = noise

        # Initilize Geometry
        self.geometry = SurgeryGeometry(distance=distance)

        # Initialize Pairings of all types (Standard, AC merging, AT merging, Surgery)
        self.master_pairings = MasterPairings(
            std_pairings=LatticeSurgeryPairings(
                qubit_coords=self.geometry.coords_ancilla
                | self.geometry.coords_control
                | self.geometry.coords_target,
                merging=False,
            ),
            ac_merge_pairings=LatticeSurgeryPairings(
                qubit_coords=self.geometry.coords_surgery,
                merging=True,
                merging_type="AC",
            ),
            at_merge_pairings=LatticeSurgeryPairings(
                qubit_coords=self.geometry.coords_surgery,
                merging=True,
                merging_type="AT",
            ),
        )

        super().__init__(noise=noise)

    def build_circuit(self) -> stim.Circuit:
        # Defining flow_circuit -> For Solving Flow Observables at the end
        flow_circuit = stim.Circuit()

        # Defining return circuit
        return_circuit = stim.Circuit()

        # Adding Reset Circuit
        if self.control_state_init not in {"Y+", "Y-"} and self.target_state_init not in {
            "Y+",
            "Y-",
        }:
            reset_circuit = self._adding_reset_circuit()
        else:
            reset_circuit, y_setup_circuit = self._adding_reset_circuit()
            flow_circuit += y_setup_circuit

        # Adding Initialization Circuit
        flow_circuit += self._adding_initialization_circuit()

        # Adding Merging Ancilla Control Circuit
        flow_circuit += self._adding_merge(type="AC")

        # Adding Splitting Ancilla Control Circuit
        flow_circuit += self._adding_split(type="AC")

        # Adding Merging Ancilla Target Circuit
        flow_circuit += self._adding_merge(type="AT")

        # Adding Splitting Ancilla Target Circuit
        flow_circuit += self._adding_split(type="AT")

        # Adding Reset to Return Circuit
        return_circuit += reset_circuit

        # Adding Flow Circuit to Return Circuit
        return_circuit += flow_circuit

        # Adding Solve Flow Observables
        return_circuit += self._adding_solve_flow_observables()

        # Adding Final Measurement Circuit
        return_circuit += self._adding_final_measurement_circuit()

        # Adding Noise Model if applicable
        return_circuit = self._apply_noise(circuit=return_circuit)

        return return_circuit

    def _adding_reset_circuit(self) -> stim.Circuit:
        # Adding Reset Circuit (If Y initilization is selected, this includes the
        # Y creation flow circuit)
        reset_circuit_builder = SurgeryReset(
            geometry=self.geometry,
        )

        return reset_circuit_builder.build_circuit()

    def _adding_initialization_circuit(self) -> stim.Circuit:
        # Adding Initialization Circuit
        init_circuit_builder = SurgeryInitialization(
            geometry=self.geometry,
            master_pairings=self.master_pairings,
        )
        return init_circuit_builder.build_circuit()

    def _adding_merge(self, type: str) -> stim.Circuit:
        # Adding Merge Circuit
        merge_circuit_builder = SurgeryMerge(
            geometry=self.geometry,
            master_pairings=self.master_pairings,
            merging_type=type,
        )
        return merge_circuit_builder.build_circuit()

    def _adding_split(self, type: str) -> stim.Circuit:
        # Adding Split Circuit
        split_circuit_builder = SurgerySplit(
            geometry=self.geometry,
            master_pairings=self.master_pairings,
            splitting_type=type,
        )
        return split_circuit_builder.build_circuit()

    def _adding_final_measurement_circuit(self) -> stim.Circuit:
        # Adding Final Measurement Circuit
        final_measure_circuit_builder = SurgeryFinalMeasure(
            geometry=self.geometry,
        )
        return final_measure_circuit_builder.build_circuit()

    def _adding_solve_flow_observables(self) -> stim.Circuit:
        pass
