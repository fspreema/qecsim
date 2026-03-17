import stim
from tqecd import annotate_detectors_automatically

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

Coord = complex
Label = str
Index = int
Pair = tuple[Coord, Coord]

__all__ = ["SurfaceBuilder"]


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
        self.logical_h = logical_h
        self.noise = noise
        self.ft_init = ft_init
        self.ft_measurements = ft_measurements

        # Setting Up Builder Parameters
        self.y_sections_required = state_init in {"+i", "-i"}

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
        full_circuit = stim.Circuit()

        # Adding Setup Resets
        full_circuit += self._adding_setup_resets()

        # Adding Initialiazion Circuit
        full_circuit += self._adding_initilization()

        # Addings Repetion Circuit
        full_circuit += self._adding_repetition()

        # Adding Conditional Circuits depending on Y Basis or Transversal H
        if self.y_sections_required:
            full_circuit += self._adding_y_basis_sections()
        elif self.logical_h is True:
            full_circuit += self._adding_logical_h_sections()

        # Adding Final Measurement Circuit -> Not for Y Basis
        if not self.y_sections_required:
            full_circuit += self._adding_final_measurement()

        # Adding Detectors
        return_circuit = self._adding_detectors(input_circuit=full_circuit)

        # Adding Noise if specified
        return_circuit = self.apply_noise(input_circuit=return_circuit, 
                                        distance=self.distance,
                                        geometry= self.master_geometry.geometry_std,
                                        noise=self.noise, 
                                        ft_init=self.ft_init, 
                                        ft_meas=self.ft_measurements)

        return return_circuit

    def get_logical_meas_rec(self) -> list[int]:
        """
        Returns the list of measurement record positions that need to be xored together
        to get the final logical measurement

        -> This is used for Caluclation of the PTM
        """

        if self.y_sections_required:
            pass
        else:
            rec_list = self._get_x_z_basis_measurement_recs()

        return rec_list

    def _adding_setup_resets(self) -> stim.Circuit:
        # Adding Reset Circuit
        reset_circ = SurfaceReset(
            master_geometry=self.master_geometry,
            type="standard" if self.state_init not in {"+i", "-i"} else "y_basis",
        )

        return reset_circ.build_circuit()

    def _adding_initilization(self) -> stim.Circuit:
        # Init Circuit depending on Y Basis or Standard
        init_circ = SurfaceInitialization(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type="standard" if self.state_init not in {"+i", "-i"} else "y_basis",
        )

        return init_circ.build_circuit()

    def _adding_repetition(self) -> stim.Circuit:
        # Repetition Circuit depending on Y Basis or Standard
        repet_circ = SurfaceRepetitionCircuit(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type="standard" if self.state_init not in {"+i", "-i"} else "y_basis",
        )

        # Setting rec_list
        self.rec_list = repet_circ.rec_list()

        return repet_circ.build_circuit()

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

    def _get_x_z_basis_measurement_recs(self) -> list[int]:
        final_meas_circ = FinalMeasureCircuit(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type=("log_h" if self.logical_h is True else "standard"),
        )

        # Get Measurement records for everything except the y basis
        rec_list = final_meas_circ.build_observable_meas_rec()

        return rec_list

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

    @staticmethod
    def _adding_detectors(input_circuit: stim.Circuit) -> stim.Circuit:
        """
        Calculates all Detectors needed by using the tqecd package
        """

        # Annotate Detectors Automatically
        annotated_circuit = annotate_detectors_automatically(
            circuit=input_circuit,
        )

        return annotated_circuit
