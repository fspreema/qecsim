import stim
from tqecd import annotate_detectors_automatically

from src.codes.surface_code_rotated.circuits.final_measure import FinalMeasureCircuit
from src.codes.surface_code_rotated.circuits.initial import SurfaceInitialization
from src.codes.surface_code_rotated.circuits.repetition import SurfaceRepetitionCircuit
from src.codes.surface_code_rotated.circuits.reset import SurfaceReset
from src.codes.surface_code_rotated.circuits.y_rev_switch import YRevSwitchCircuit
from src.codes.surface_code_rotated.circuits.y_switch import YSwitchCircuit
from src.codes.surface_code_rotated.data_geometry import MasterGeometry, MasterPairings
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

        super().__init__(noise=noise)

    def build_circuit(self) -> stim.Circuit:
        # Initialize Empty Circuit
        self.full_circuit = stim.Circuit()

        ######################
        # Adding Reset Circuit
        ######################
        reset_circ = SurfaceReset(
            master_geometry=self.master_geometry,
            type="standard"
            if self.logical_h is False and self.state_init not in {"+i", "-i"}
            else "log_h"
            if self.logical_h
            else "ybasis",
        )
        self.full_circuit += reset_circ.build_circuit()

        ##############################
        # Adding Initilization Circuit
        ##############################

        # Init Circuit depending on Y Basis or Standard
        init_circ = SurfaceInitialization(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type="standard"
            if self.logical_h is False and self.state_init not in {"+i", "-i"}
            else "log_h_initial"
            if self.logical_h
            else "y_initial",
        )
        self.full_circuit += init_circ.build_circuit()

        ###########################
        # Adding Repetition Circuit
        ###########################

        # Repetition Circuit depending on Y Basis or Standard
        repet_circ = SurfaceRepetitionCircuit(
            master_geometry=self.master_geometry,
            master_pairings=self.master_pairings,
            type="y_repetition" if self.state_init in {"+i", "-i"} else "standard",
        )
        self.full_circuit += repet_circ.build_circuit()
        self.rec_list = repet_circ.rec_list()

        #####################################################
        # Y BASIS ONLY: Add Switch/Memory/Rev-Switch Circuits
        #####################################################

        if self.state_init in {"+i", "-i"}:
            # Adding Y Switch Circuit
            y_switch_circ = YSwitchCircuit(
                master_geometry=self.master_geometry,
                master_pairings=self.master_pairings,
            )
            self.full_circuit += y_switch_circ.build_circuit()

            # Adding Y Memory Circuit
            y_memory_circ = SurfaceRepetitionCircuit(
                master_geometry=self.master_geometry,
                master_pairings=self.master_pairings,
                type="y_memory",
            )
            self.full_circuit += y_memory_circ.build_circuit()
            self.rec_list += y_memory_circ.rec_list()

            # Adding Y Reverse Switch Circuit
            y_rev_switch_circ = YRevSwitchCircuit(
                master_geometry=self.master_geometry,
                master_pairings=self.master_pairings,
            )
            self.full_circuit += y_rev_switch_circ.build_circuit()

            # Adding another Repetition Circuit after Y Basis Memory for fault tolerance
            repet_circ_2 = SurfaceRepetitionCircuit(
                master_geometry=self.master_geometry,
                master_pairings=self.master_pairings,
                type="y_repetition",
            )
            self.full_circuit += repet_circ_2.build_circuit()

        ########################################################################
        # H GATE ONLY: Add Flipped Circuits: Flipped Init and Flipped Repetition
        ########################################################################

        if self.logical_h is True:
            # Adding Inital Flipped Circuit
            flip_init_circ = SurfaceInitialization(
                master_geometry=self.master_geometry,
                master_pairings=self.master_pairings,
                type="log_h_initial",
            )
            self.full_circuit += flip_init_circ.build_circuit()

            # Adding Flipped Repetition Circuit
            flip_repet_circ = SurfaceRepetitionCircuit(
                master_geometry=self.master_geometry,
                master_pairings=self.master_pairings,
                type="h_repetition",
            )
            self.full_circuit += flip_repet_circ.build_circuit()

        #####################################################
        # Adding Final Measurement Circuit -> Not for Y Basis
        #####################################################

        if self.state_init not in {"+i", "-i"}:
            final_meas_circ = FinalMeasureCircuit(
                master_geometry=self.master_geometry,
                master_pairings=self.master_pairings,
                type=("logical_h" if self.logical_h is True else "standard"),
            )
            meas_circ, self.rec_list = final_meas_circ.build_final_measurement_circuit()
            self.full_circuit += meas_circ

        ###########################
        # Adding Noise if specified
        ###########################
        self.full_circuit = self._apply_noise(self.full_circuit)

        return self.full_circuit, self.rec_list

    def _adding_detectors(self) -> stim.Circuit:
        # Annotate Detectors Automatically
        annotated_circuit = annotate_detectors_automatically(
            circuit=self.full_circuit,
        )

        return annotated_circuit, self.rec_list
