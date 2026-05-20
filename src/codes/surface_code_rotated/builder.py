import stim

from src.codes.surface_code_rotated.circuits.final_measure import FinalMeasureCircuit
from src.codes.surface_code_rotated.circuits.initial import SurfaceInitialization
from src.codes.surface_code_rotated.circuits.repetition import SurfaceRepetitionCircuit
from src.codes.surface_code_rotated.circuits.reset import SurfaceReset
from src.codes.surface_code_rotated.circuits.y_rev_switch import YRevSwitchCircuit
from src.codes.surface_code_rotated.circuits.y_switch import YSwitchCircuit
from src.codes.surface_code_rotated.data_geometry import MasterGeometry, MasterPairings
from src.codes.surface_code_rotated.get_flows import YBasisGetCircuitFlows
from src.codes.surface_code_rotated.get_stab_pairings import SurfacePairings
from src.codes.surface_code_rotated.surface_geom import SurfaceGeometry
from src.core.base_class_builder import BaseClassBuilder
from src.core.data_models import NoiseParameters
from src.core.measurement_tracker import MeasurementTracker

Coord = complex
Label = str
Index = int
Pair = tuple[Coord, Coord]

__all__ = ["SurfaceBuilder"]


# Both currently set to false as detector construciton of ft y-basis is not
# correct and therefore currently not working!
FT_INIT = False
FT_MEAS = False

# Currently also not tested properly!
LOGICAL_H = False


class SurfaceBuilder(BaseClassBuilder):
    def __init__(
        self,
        distance: int,
        state_init: str,
        log_obs: str,
        logical_h: bool = False,
        noise: NoiseParameters = None,
        ft_init: bool = False,
        ft_measurements: bool = False,
    ):
        """
        Generates Rotated-Surface-Code

        Args:
            Distance (int): Distance of the surface code i.e. lattice size
            Rounds (int): Number of syndrome measurement rounds per Shot
            noise (float): Probability for x y and z error in Pauli-Channel
            state_init (str): Initial State of the logical qubit: "0", "1", "+", "-", "+i", "-i"
            log_obs (str): Logical Observable to be measured: "X", "Y", "Z"
            logical_h (bool): Whether to add transversal H section for logical H state preparation
            -> THIS ALSO DOES NOT REALLY WORK AS THE QUBITS NEED TO BE ROTATED AFTER THE 
                TRANSVERSAL H (Currently just flipping X and Z type)
            ft_init (bool): If True -> Perform inplace Y-access (Gidney) DOES NOT WORK CURRENTLY!!

        Information:
            Logical Operator is Z Operator and pre-Defined!
            Logical State: 0 -> Only z Stabilizer detectors in the first round as x detectors
                                are non deterministc for the first round
                                (STILL: COMPLETE MEASUREMENT)
                            -> In theory we don not even need to meassure the X stabilizers
                                at all because we do not have phase errors
                                (Would result in global phases which can be ignored)

            Noise-Model: Analog to Stims Circuit i.e. Full Noise Model implemented
                        -> Before Round Depolarization Data
                        -> Before Measurement Flip Probability
                        -> After Clifford Depolarization
                        -> After Reset Flip Propability

        Returns:
            stim.Circuit: Comiled Circuit in Stim format
        """
        # Init Paramters
        self.distance = distance
        self.state_init = state_init
        self.log_obs = log_obs
        self.logical_h = LOGICAL_H
        self.noise = noise
        self.ft_init = FT_INIT
        self.ft_measurements = FT_MEAS

        # Setting Up tick information for noise model construction
        # Used for construction of the noise model
        # -> Num Ticks Until d rounds have passed
        # -> Num Tikcs Until only d rounds are left (i.e. d noisy rounds have passed)
        self.tick_dict : dict[int, str] = {}

        # Setting Up Builder Parameters
        self.y_sections_required = state_init in {"Y+", "Y-"} and self.ft_init

        # Initialize Geometries for standard and y-basis
        self.master_geometry = MasterGeometry(
            geometry_std=SurfaceGeometry(
                distance=distance,
                state_init=state_init,
                logical_observable=log_obs,
            ),
            geometry_ybasis=SurfaceGeometry(
                distance=distance,
                state_init=state_init,
                logical_observable=log_obs,
                y_basis=True,
            ),
        )

        # Initiate Measurement Tracker
        self.tracker = MeasurementTracker(geometry=self.master_geometry.geometry_std)

        # Initialize Pairings for standard and y-basis
        self.master_pairings = MasterPairings(
            pairings_std=SurfacePairings(
                patch=self.master_geometry.geometry_std.coords,
                distance=distance,
            ),
            pairings_ybasis=SurfacePairings(
                patch=self.master_geometry.geometry_ybasis.coords,
                distance=distance,
                y_basis=True,
            ),
            pairings_ymemory=SurfacePairings(
                patch=self.master_geometry.geometry_ybasis.coords,
                distance=distance,
                y_basis=True,
                y_memory=True,
            ),
            pairings_yswitch=SurfacePairings(
                patch=self.master_geometry.geometry_ybasis.coords,
                distance=distance,
                y_basis=True,
                y_switch=True,
            ),
            pairings_log_h=SurfacePairings(
                patch=self.master_geometry.geometry_std.coords,
                distance=distance,
                is_flipped=True,
            ),
        )

    def build_circuit(self) -> stim.Circuit:
        # Initialize Empty Circuit
        self.return_circuit = stim.Circuit()

        # Adding Setup Resets
        self.return_circuit += self._adding_setup_resets()

        # Adding Initialiazion Circuit
        self.return_circuit += self._adding_initialization()

        # Addings Repetion Circuit
        self.return_circuit += self._adding_repetition()

        # Adding Conditional Circuits depending on Y Basis or Transversal H
        if self.y_sections_required:
            self.return_circuit += self._adding_y_basis_sections()
        elif self.logical_h is True:
            self.return_circuit += self._adding_logical_h_sections()

        # Adding Final Measurement Circuit -> Not for FT Y Basis
        if not self.y_sections_required:
            self.return_circuit += self._adding_final_measurement()

        # Adding Noise if specified
        self.return_circuit = self.apply_noise(
            input_circuit=self.return_circuit,
            geometry=self.master_geometry.geometry_std,
            noise=self.noise, 
            ft_init=FT_INIT, 
            ft_meas=FT_MEAS,
            num_tick_first_noise=self._get_noisy_tick(beginning=True),
            num_tick_last_noise=self._get_noisy_tick(beginning=False),
        )

        return self.return_circuit
    
    def _curr_type(self):
        # Determine type of circuit
        if not self.ft_init and self.state_init in {"Y+", "Y-"}:
            curr_type = "non_ft_init"
        elif self.ft_init and self.state_init in {"Y+", "Y-"}:
            curr_type = "y_basis"
        else:
            curr_type = "standard"

        return curr_type

    def _adding_setup_resets(self) -> stim.Circuit:

        # Adding Reset Circuit
        reset_builder = SurfaceReset(
            master_geometry=self.master_geometry,
            type=self._curr_type(),
        )

        reset_circ = reset_builder.build_circuit()

        # Adding TICK info
        self.tick_dict["reset"] = reset_circ.num_ticks

        return reset_circ

    def _adding_initialization(self) -> stim.Circuit:
        # Updating Measurement Tracker
        self.tracker.add_previous_measurements(
            count=self.return_circuit.num_measurements,
        )

        # Init Circuit depending on Y Basis or Standard
        init_builder = SurfaceInitialization(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type=self._curr_type(),
            tracker=self.tracker,
        )

        init_circ = init_builder.build_circuit()

        # Adding TICK info
        self.tick_dict["init"] = init_circ.num_ticks

        return init_circ

    def _adding_repetition(self) -> stim.Circuit:
        # Repetition Circuit depending on Y Basis or Standard
        repet_builder = SurfaceRepetitionCircuit(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type=self._curr_type(),
            tracker=self.tracker,
        )

        # Setting rec_list and getting circ
        repet_circ = repet_builder.build_circuit()
        self.rec_list = repet_builder.rec_list()

        # Adding TICK info
        self.tick_dict["repetition_per_round"] = repet_circ.num_ticks // (self.distance + 2)

        return repet_circ

    def _adding_y_basis_sections(self) -> stim.Circuit:
        # Y BASIS ONLY: Add Switch/Memory/Rev-Switch Circuits
        y_section_circuits = stim.Circuit()

        # Adding Y Switch Circuit
        y_switch_circ = YSwitchCircuit(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
        )
        y_section_circuits += y_switch_circ.build_circuit()

        # Adding Y Memory Circuit
        y_memory_circ = SurfaceRepetitionCircuit(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type="y_memory",
        )
        y_section_circuits += y_memory_circ.build_circuit()

        # Updating rec_list
        self.rec_list += y_memory_circ.rec_list()

        # If logical FLow is Y, add Rev Switch and get flows and repetition round (fault tolerance)
        # -> For X and Z measurements the measurement is not fault tolerant either way so
        #    no need to add extra rounds
        if self.log_obs == "Y":
            # Adding Y Reverse Switch Circuit
            y_rev_switch_circ = YRevSwitchCircuit(
                master_geometry=self.master_geometry,
                master_pairings=self.master_pairings,
            )
            y_section_circuits += y_rev_switch_circ.build_circuit()

            y_flow_getter = YBasisGetCircuitFlows(master_geometry=self.master_geometry)
            logical_creation_circ = y_flow_getter.get_flows(
                logical_creation_circ=y_switch_circ.build_circuit(),
                logical_contraction_circ=y_rev_switch_circ.build_circuit(),
                circuits_between=[y_memory_circ.build_circuit()],
            )
            y_section_circuits += logical_creation_circ

            # Adding another Repetition Circuit after Y Basis Memory for fault tolerance
            repet_circ_2 = SurfaceRepetitionCircuit(
                master_geometry=self.master_geometry,
                master_pairings=self.master_pairings,
                type="y_basis",
            )
            y_section_circuits += repet_circ_2.build_circuit()

        return y_section_circuits

    def _adding_logical_h_sections(self) -> stim.Circuit:
        # H GATE ONLY: Add Flipped Circuits: Flipped Init and Flipped Repetition

        h_section_circuits = stim.Circuit()

        # Adding Inital Flipped Circuit
        flip_init_circ = SurfaceInitialization(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type="log_h",
        )
        h_section_circuits += flip_init_circ.build_circuit()

        # Adding Flipped Repetition Circuit
        flip_repet_circ = SurfaceRepetitionCircuit(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type="log_h",
        )
        h_section_circuits += flip_repet_circ.build_circuit()

        return h_section_circuits

    def _adding_final_measurement(self) -> stim.Circuit:
        # Adding Final Measurement Circuit -> Not for Y Basis

        final_meas_circ = FinalMeasureCircuit(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type=("log_h" if self.logical_h is True else "standard"),
        )

        # Build measurement Circuit
        meas_circ = final_meas_circ.build_final_measurement_circuit()

        return meas_circ
    
    def _get_noisy_tick(self, beginning: bool) -> int:
        """
        Get the number of ticks until noise should be applied for the first time 
        (if beginning = True) or the last time (if beginning = False)
        """

        if beginning:
            num_tick = self.tick_dict.get("reset") \
                        + self.tick_dict.get("init")\
                        + self.tick_dict.get("repetition_per_round") * 1
            
        else:
            num_tick = self.tick_dict.get("reset") \
                    + self.tick_dict.get("init") \
                    + self.tick_dict.get("repetition_per_round") * (self.distance + 1) \
                    - 1
            
        return num_tick
