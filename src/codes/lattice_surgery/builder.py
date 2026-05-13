import stim

from src.codes.lattice_surgery.circuits.final_measure import SurgeryFinalMeasure
from src.codes.lattice_surgery.circuits.initial import SurgeryInitialization
from src.codes.lattice_surgery.circuits.merge import SurgeryMerge
from src.codes.lattice_surgery.circuits.pauli_observables import SurgeryPauliObservables
from src.codes.lattice_surgery.circuits.reset import SurgeryReset
from src.codes.lattice_surgery.circuits.split import SurgerySplit
from src.codes.lattice_surgery.data_geometry import MasterPairings
from src.codes.lattice_surgery.get_flows import SurgeryFlowObservables
from src.codes.lattice_surgery.get_stab_pairings import LatticeSurgeryPairings
from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.core.base_class_builder import BaseClassBuilder
from src.core.data_models import NoiseParameters
from src.core.measurement_tracker import MeasurementTracker

__all__ = ["SurgeryBuilder"]

# Currently fixed as FT also fixed to False
FT_INIT = False
FT_MEASUREMENT = False

class SurgeryBuilder(BaseClassBuilder):

    def __init__(
        self,
        distance: int,
        control_state_init: str,
        target_state_init: str,
        control_measure_basis: str,
        target_measure_basis: str,
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
        self.control_measure_basis = control_measure_basis
        self.target_measure_basis = target_measure_basis
        self.noise = noise
        self.return_circuit = stim.Circuit()

        # Setting Up tick information for noise model construction
        # Used for construction of the noise model
        # -> Num Ticks Until d rounds have passed
        # -> Num Tikcs Until only d rounds are left (i.e. d noisy rounds have passed)
        self.tick_dict : dict[int, str] = {}

        # Initilize Geometry
        self.geometry = SurgeryGeometry(
            distance=distance,
            control_state_init=control_state_init,
            target_state_init=target_state_init,
        )

        # Initialize Pairings of all types (Standard, AC merging, AT merging, Surgery)
        self.master_pairings = MasterPairings(
            std_pairings=LatticeSurgeryPairings(
                qubit_coords=self.geometry.get_coords(specific_coord=None),
                merging=False,
            ),
            ac_merge_pairings=LatticeSurgeryPairings(
                qubit_coords=self.geometry.get_coords(specific_coord=None),
                merging=True,
                merging_type="AC",
            ),
            at_merge_pairings=LatticeSurgeryPairings(
                qubit_coords=self.geometry.get_coords(specific_coord=None),
                merging=True,
                merging_type="AT",
            ),
        )

        # Initiate Measurement Tracker
        self.tracker = MeasurementTracker(geometry=self.geometry)

        # Initlize Current Flow
        self.curr_flow = (
            f"{self.geometry.control_state_init[0]}{self.geometry.target_state_init[0]}"
            " -> " + f"{self.control_measure_basis[0]}{self.target_measure_basis[0]}"
        )

    def build_circuit(self) -> stim.Circuit:
        # Defining flow_circuit -> For Solving Flow Observables at the end
        flow_circuit = stim.Circuit()

        # Adding Reset Circuit
        reset_circuit = self._adding_reset_circuit()

        # Adding Initialization Circuit
        flow_circuit += self._adding_initialization_circuit()

        # If not valid flow, we need to XOR rec measurements with pauli observables
        # to get deterministic outcomes
        if not self._valid_flow():
            flow_circuit += self._add_pauli_observables(flow_type="incoming_flow")

        # Adding Merging Ancilla Control Circuit
        flow_circuit += self._adding_merge(merge_type="AC")

        # Adding Splitting Ancilla Control Circuit
        flow_circuit += self._adding_split(split_type="AC")

        # Adding Merging Ancilla Target Circuit
        flow_circuit += self._adding_merge(merge_type="AT")

        # Adding Splitting Ancilla Target Circuit
        flow_circuit += self._adding_split(split_type="AT")

        # If not valid flow, we need to XOR rec measurements with pauli observables
        # to get deterministic outcomes
        if not self._valid_flow():
            flow_circuit += self._add_pauli_observables(flow_type="outgoing_flow")

        # Adding Reset to Return Circuit
        self.return_circuit += reset_circuit

        # Adding Flow Circuit to Return Circuit
        self.return_circuit += flow_circuit

        # Adding Solve Flow Observables if valid cx flow is selected
        if self._valid_flow():
            self.return_circuit += self._adding_solve_flow_observables(flow_circuit=flow_circuit)

        # Adding Final Measurement Circuit
        self.return_circuit += self._adding_final_measurement_circuit()

        # Adding Noise Model if applicable
        self.return_circuit = self.apply_noise(
            input_circuit=self.return_circuit,
            distance=self.distance,
            geometry=self.geometry,
            noise=self.noise,
            ft_init=FT_INIT,
            ft_meas=FT_MEASUREMENT,
            num_tick_first_noise=self._get_noisy_tick(beginning=True),
            num_tick_last_noise=self._get_noisy_tick(beginning=False),
        )

        return self.return_circuit

    def _adding_reset_circuit(self) -> stim.Circuit:
        # Building Reset Circuit
        reset_circuit_builder = SurgeryReset(
            geometry=self.geometry,
        )
        self.reset_circuit = reset_circuit_builder.build_circuit()

        # Adding TICK info
        self.tick_dict["reset"] = self.reset_circuit.num_ticks

        return self.reset_circuit

    def _adding_initialization_circuit(self) -> stim.Circuit:
        # Count all previous measurements form init and reset circuits
        if isinstance(self.reset_circuit, tuple):
            num_measurements_reset = (
                self.reset_circuit[0].num_measurements + self.reset_circuit[1].num_measurements
            )
        else:
            num_measurements_reset = self.reset_circuit.num_measurements

        # Updating Measurement Tracker
        self.tracker.add_previous_measurements(
            count=num_measurements_reset,
        )

        # Adding Initialization Circuit
        init_circuit_builder = SurgeryInitialization(
            geometry=self.geometry,
            master_pairings=self.master_pairings,
            measurement_tracker=self.tracker,
        )
        self.init_circuit = init_circuit_builder.build_circuit()

        # Adding TICK info
        self.tick_dict["init_single"] = 14
        self.tick_dict["init_repeat_per_round"] = (
            (self.init_circuit.num_ticks - self.tick_dict["init_single"]) //
            ((2 * self.distance) - 1)
        )

        return self.init_circuit

    def _adding_merge(self, merge_type: str) -> stim.Circuit:
        # Adding Merge Circuit
        merge_circuit_builder = SurgeryMerge(
            geometry=self.geometry,
            master_pairings=self.master_pairings,
            merging_type=merge_type,
            tracker=self.tracker,
        )
        merge_circuit = merge_circuit_builder.build_circuit()

        # Adding TICK info
        if merge_type == "AT":
            self.tick_dict["merge_full_AT"] = merge_circuit.num_ticks

        else:
            self.tick_dict["merge_full_AC"] = merge_circuit.num_ticks
        

        return merge_circuit

    def _adding_split(self, split_type: str) -> stim.Circuit:
        # Adding Split Circuit
        split_circuit_builder = SurgerySplit(
            geometry=self.geometry,
            master_pairings=self.master_pairings,
            split_type=split_type,
            tracker=self.tracker,
        )
        split_circuit = split_circuit_builder.build_circuit()

        # Adding TICK info
        self.tick_dict["split_init"] = 14
        if split_type == "AC":
            self.tick_dict["split_repeat_AC"] = (
                (split_circuit.num_ticks - self.tick_dict["split_init"])
            )
        else:
            self.tick_dict["split_repeat_per_round_AT"] = (
                (split_circuit.num_ticks - self.tick_dict["split_init"]) // ((self.distance * 2) - 1)
            )

        return split_circuit

    def _adding_final_measurement_circuit(self) -> stim.Circuit:
        # Adding Final Measurement Circuit
        final_measure_circuit_builder = SurgeryFinalMeasure(
            geometry=self.geometry,
            curr_flow=self.curr_flow,
            control_measure_basis=self.control_measure_basis,
            target_measure_basis=self.target_measure_basis,
        )
        final_measure_circuit = final_measure_circuit_builder.build_circuit()

        return final_measure_circuit

    def _adding_solve_flow_observables(self, flow_circuit: stim.Circuit) -> stim.Circuit:
        # Adding Solve Flow Observables Circuit
        flow_getter = SurgeryFlowObservables(
            geometry=self.geometry,
        )
        return flow_getter.get_observable_from_flow(
            flow_circuit=flow_circuit,
            curr_flow=self.curr_flow,
        )

    def _add_pauli_observables(self, flow_type: str) -> stim.Circuit:
        # Adding Pauli Observables for both control and target patches
        observable_getter = SurgeryPauliObservables(
            geometry=self.geometry,
            control_measure_basis=self.control_measure_basis,
            target_measure_basis=self.target_measure_basis,
            control_state_init=self.geometry.control_state_init,
            target_state_init=self.geometry.target_state_init,
            flow_type=flow_type,
        )

        return observable_getter.build_circuit()

    def _valid_flow(self) -> bool:
        # Define valid flows
        valid_flows = [
            "XI -> XX",
            "XX -> XI",
            "IX -> IX",
            "IZ -> ZZ",
            "ZZ -> IZ",
            "ZI -> ZI",
            "ZX -> ZX",
            "YI -> YX",
            "YY -> XZ",
            "IY -> ZY",
            "XY -> YZ",
            "XZ -> YY",
            "YZ -> XY",
            "YX -> YI",
            "ZY -> IY",
        ]

        if self.curr_flow not in valid_flows:
            return False
        else:
            return True
        
    def _get_noisy_tick(self, beginning: bool) -> int:
    
        """
        Get the number of ticks until noise should be applied for the first time 
        (if beginning = True) or when it should be applied
        for the last time (if beginning = False)
        """

        if beginning:
            num_tick = self.tick_dict.get("reset") \
                    + self.tick_dict.get("init_single") \
                    + self.tick_dict.get("init_repeat_per_round") * 1\
                    + 1 \
            
        else:
            num_tick = self.tick_dict.get("reset") \
                    + self.tick_dict.get("init_single") \
                    + (self.tick_dict.get("init_repeat_per_round") * ((2 * self.distance) - 1))\
                    + self.tick_dict.get("merge_full_AT")\
                    + self.tick_dict.get("merge_full_AC")\
                    + self.tick_dict.get("split_init") * 2\
                    + self.tick_dict.get("split_repeat_AC")\
                    + (self.tick_dict.get("split_repeat_per_round_AT") * ((self.distance * 1) - 1))\
                    - 150 \

        return num_tick
                    
