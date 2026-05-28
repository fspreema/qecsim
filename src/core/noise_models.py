from abc import ABC, abstractmethod

import numpy as np
import stim

from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.codes.surface_code_rotated.surface_geom import SurfaceGeometry
from src.codes.xzzx.xzzx_geom import XZZXGeometry

CLIFFORD_OPERATIONS = ["H", "CX", "S", "S_DAG", "CZ", "XCY", "SQRT_X_DAG"]
MEASUREMENT_OPERATIONS = ["M", "MX", "MY"]
RESET_OPERATIONS = ["R", "RX", "RY", "RZ"]

class NoiseModel(ABC):

    def __init__(self, 
                 circuit: stim.Circuit, 
                 noise: dict, 
                 geometry: SurfaceGeometry | SurgeryGeometry | XZZXGeometry,
                 num_tick_first_noise: int,
                 num_tick_last_noise: int,
                 ft_init: bool = True, 
                 ft_measurements: bool = True):
        
        """
        Abstract Base Class needed for the indivudal noise classes

        Args:
            Everything self explanatory except...
            -> custom_rounds (dict): If duration of one round changes during the circuit
                i.e., Lattice Surgery circuit (Merge and Split vs intial) pass dictionary
                including the individual rounds marked by the starting TICK number of that
                round
                -> CURRENTLY NOT USED DUE TO BEFORE ROUND DEPOL NOT IMPLEMENTED FOR PENOM NOISE
        """
        
        # Set up Preliminary Information
        self.circuit = circuit
        self.noise = noise
        self.ft_init = ft_init
        self.ft_measurements = ft_measurements
        self.curr_tick = 0
        self.num_tick_last_noise = num_tick_last_noise
        self.num_tick_first_noise = num_tick_first_noise
        self.data_qubits = geometry.data_idx
        self.bias = self.noise.get("bias", None)
        
        # Check type of Circuit
        self.is_lattice_surgery = isinstance(geometry, SurgeryGeometry)

        # Check what Type of Circuit is used
        if self.is_lattice_surgery:
            # Load Lattice Surgery Geometry
            self.shared_qubits = set(geometry.anc_x_bdy_b_stb_idx + geometry.anc_z_bdy_r_stb_idx)
            self.ancilla_qbts = set(geometry.anc_data_idx + geometry.anc_x_stb_idx + geometry.anc_z_stb_idx)
            self.unshared_ancilla = self.ancilla_qbts - self.shared_qubits
            self.all_qubits = geometry.data_idx + geometry.stab_idx

        elif isinstance(geometry, SurfaceGeometry) or isinstance(geometry, XZZXGeometry):
            # Load normal Surface Geometry
            self.all_qubits = geometry.data_idx + geometry.stab_idx

        else:
            raise ValueError(f"Current geometry of type {type(geometry)} is unsupported!")

    def apply(self) -> stim.Circuit:

        # Reset tick count
        self.curr_tick = 0
        
        # Create Noisy Circuit
        noisy_circuit = self._apply_recursive_noise_operations(circuit=self.circuit)

        return noisy_circuit
    
    def _get_idx_of_noisy_qbts(self, qbt_lst: list[int]) -> list[int]:
        """
        Returns for a given set of qubit idx that the current operation acts on
        the cleaned up version of all qbts that indeed need noise applied

        -> This is used in the case for lattice surgery as during the same TICK
            some qbt idx may already be in the next circuit round while others are 
            not (Therefore needed for adding noise at the beginning and end)
        
        ###Example###
        0 2 3 5 6 3 2 -> 0 2 6 2
        """

        # If not lattice surgery patch -> skip filtering!
        if not self.is_lattice_surgery:
            return qbt_lst

        noisy_qbt_idx: list[int] = []

        # If current TICK is at the start of noise application, only ANCILLA PATCH
        # can be applicable to noise while control and target are still on old round
        # ATTENTION: Shared qubit cannot be applicable to noise as it is used by control
        #           target and ancilla!

        # If ft_init is False, apply noise to 6 ticks earlier as num_tick_first_noise for ancilla
        # patch -> Begins new round wehn control target are still in old round
        if not self.ft_init and self.num_tick_first_noise <= self.curr_tick <= self.num_tick_first_noise + 5:
            for curr_qbt in qbt_lst:
                # Noise should only be applied for the non shared ancilla qubits during that time
                if curr_qbt not in self.shared_qubits and curr_qbt in self.ancilla_qbts:
                    noisy_qbt_idx.append(curr_qbt)
            return noisy_qbt_idx
        
        if not self.ft_measurements and self.num_tick_last_noise - 5 <= self.curr_tick <= self.num_tick_last_noise:
            for curr_qbt in qbt_lst:
                # Noise should be applied to everything except the ancilla qubits
                if curr_qbt not in self.unshared_ancilla:
                    noisy_qbt_idx.append(curr_qbt)
            return noisy_qbt_idx
        
        # If not in these specific Windows return all qbts as nosiy
        return qbt_lst

    def _prepare_noise_targets(self, instruction) -> tuple[stim.CircuitInstruction, list[int]]:
        """
        For any instruction this helper does the following:

        1) Retrieve all qubit idx of current instruction
        2) Filter out the idx which should get noise applied according to _get_idx_of_noisy_qubit
            (This is only really relevant for Lattice Surgery)
        3) Create the new instruction set to give to noise application helpers with only the filtered
            qbt idx list -> OUTPUT: new_instruction
        4) Get the list of all qubits not currently involved in any operations (Idling qubits)
        5) Filter the idling qubits according to _get_idx_of_noisy_qubit
            -> OUTPUT: untouched_noisy
        """

        # Retrieve all qubit indices which should get noise applied &
        # build new instruction set
        qubit_indices = [t.value for t in instruction.targets_copy()]
        noisy_qbts = self._get_idx_of_noisy_qbts(qubit_indices)
        new_instruction = self._modify_instruction_set(instruction, noisy_qbts)

        # Get untouched qubits and fiulter which should get noise applied
        noisy_set = set(noisy_qbts)
        untouched = [q for q in self.all_qubits if q not in noisy_set]
        untouched_noisy = self._get_idx_of_noisy_qbts(untouched)

        return new_instruction, untouched_noisy

    @staticmethod
    def _modify_instruction_set(old_instruction: stim.CircuitInstruction, new_qbt_lst: list[int]) -> stim.CircuitInstruction:
        """
        Static method to modify a given stim instruction set to remove qbt idx not present in new_qbt_lst
        
        I think I need to modify this for record based gates but these are not included so I do not care...
        """

        allowed_qubits = set(new_qbt_lst)
        new_targets = []
        
        for target in old_instruction.targets_copy():
            # Check integer value for standard gates
            if target.value in allowed_qubits:
                new_targets.append(target)

        # Reconstruct and return the modified instruction
        return stim.CircuitInstruction(
            old_instruction.name,
            new_targets,
            old_instruction.gate_args_copy()
        )

    def _apply_recursive_noise_operations(self, circuit: stim.Circuit) -> stim.Circuit:
        # Initilizing noisy circuit
        noisy_circuit = stim.Circuit()

        for instruction in circuit:
            # Filter Out Instructions
            if isinstance(instruction, stim.CircuitInstruction):
                # Check what type of Operation we have
                if instruction.name == "MR":
                    raise NotImplementedError("MR Operations are not supported")

                if instruction.name in CLIFFORD_OPERATIONS:
                    # Add Clifford gate
                    noisy_circuit.append(instruction)

                    # Add Noise After Gate and Check what Qubits need to be noisy
                    if self._should_apply_noise():
                        # Get Needed Noisy Infomration
                        new_instruction, untouched_noisy_qbts = self._prepare_noise_targets(instruction)

                        # Add Noise depending on NoiseModel
                        self._append_clifford_noise(noisy_circuit, new_instruction)
                        self._append_idling_error(noisy_circuit, untouched_noisy_qbts)

                elif instruction.name in MEASUREMENT_OPERATIONS:

                    # 1. Evaluate noise condition
                    apply_noise = self._should_apply_noise()

                    if apply_noise:
                        # Get Needed Noisy Infomration
                        new_instruction, untouched_noisy_qbts = self._prepare_noise_targets(instruction)

                        # Add Noise Before Measurement -> Measurement Flip
                        self._append_before_measurement_noise(noisy_circuit, new_instruction)


                    """
                    In theory one could also use M(prob) to simulate a faulty measurement but this
                    would require dismanetling the measurement operation and sometimes adding the clean
                    measurement but sometimes the noisy one...
                    -> As all measurements are done in the Z basis as H gates are used for transformation
                        X errors can easily simulate measurement flips and are therefore used here in this case...
                    """

                    # 3. Add Measurement (Always happens, sandwiched in the middle)
                    noisy_circuit.append(instruction)

                    if apply_noise:
                        # Add noise depending on Noise Model
                        # Add after measurement noise if needed
                        self._append_after_measurement_noise(noisy_circuit, new_instruction)
                        # Add Idle Noise
                        self._append_idling_error(noisy_circuit, untouched_noisy_qbts)
                        # Add additional Noise as these Qubits are currently not Measured
                        self._append_idling_error(noisy_circuit, untouched_noisy_qbts, waiting_for_r_m=True)

                elif instruction.name in RESET_OPERATIONS:
                    # Add Reset
                    noisy_circuit.append(instruction)

                    # Add Noise after Reset and Check what Qubits need to be noisy
                    if self._should_apply_noise():

                        # Get Needed Noisy Infomration
                        new_instruction, untouched_noisy_qbts = self._prepare_noise_targets(instruction)

                        # Add Noise depending on Noise Model
                        self._append_reset_noise(noisy_circuit, new_instruction)
                        self._append_before_round_depol(noisy_circuit, new_instruction)

                        # Add Idling error on all qubits not part of the current instruction:
                        self._append_idling_error(noisy_circuit, untouched_noisy_qbts)

                        # Add additional Noise as these Qubits are currently not Measured
                        self._append_idling_error(noisy_circuit, untouched_noisy_qbts, waiting_for_r_m = True)

                elif instruction.name == "TICK":
                    # Track Ticks
                    self.curr_tick += 1
                    noisy_circuit.append(instruction)

                else:
                    noisy_circuit.append(instruction)

            # Filter Out Repeat-Blocks
            elif isinstance(instruction, stim.CircuitRepeatBlock):
                body = instruction.body_copy()
                inner_noisy = self._apply_recursive_noise_operations(body)
                noisy_circuit += inner_noisy * instruction.repeat_count

                # Count ticks in inner noisy circuit
                # -> If nested blocks are used, then the if condition needs 
                # to include isinstance check
                ticks_inner_bdy = sum(1 for instr in body if instr.name == "TICK")
                
                # Update tick count as the TICK of inner_noiusy was only counted once
                self.curr_tick += (instruction.repeat_count - 1) * ticks_inner_bdy

            else:
                noisy_circuit.append(instruction)

        return noisy_circuit

    def _should_apply_noise(self) -> bool:
        """
        Noise gating rule:
        - If ft_measurements is False, do not add noise to the final small measurement cluster.
        - If ft_init is enabled, apply noise from the start.
        - Otherwise, only apply noise once the first "noisy reset" threshold is reached.
        """
        # Block everything before the transition window starts
        if not self.ft_init and self.curr_tick < self.num_tick_first_noise:
            return False
        
        # Block everything after the transition window ends
        if not self.ft_measurements and self.curr_tick > self.num_tick_last_noise:
            return False
    
        return True
    
    @staticmethod
    def _create_bias_noise_model(bias: list[float], bias_prob: float):
        """
        Creates the noise needed fot the Curstom Pauli Channels
        -> Given a bias list and full chance of a physical_error
        -> type of bias is [b_x, b_y, b_z] with b_x + b_y + b_z = 1
        -> Returns:
            - List of 3 probabilities for PAULI_CHANNEL_1
            - List of 15 probabilities for PAULI_CHANNEL_2

        -> These can then be used on the existing noise models to replace the
            depol1 and depol2 channels to incooporate biased noise!
        """

        # Check if Prob under 3/4 else BLoch sphere turned inside out
        if bias_prob > 3 / 4:
            raise ValueError("Probability for custom Pauli channel too high (> 3/4)")

        if np.any(bias):
            #######################################
            # Adding Noise for single Pauli Channel
            #######################################

            pauli_probs_single = [
                (bias_prob / sum(bias)) * bias[i] for i in range(3)
            ]

            ######################################
            # Adding Noise for multi Pauli Channel
            ######################################

            bx, by, bz = bias
            single_probs = np.array([1, bx, by, bz])

            pauli_probs_multi_unnorm: list = []

            # Probabilities for I{I,X,Y,Z}
            pauli_probs_multi_unnorm += list(single_probs)

            # Remove II prob.
            pauli_probs_multi_unnorm.pop(0)

            # Probabilities for X{I,X,Y,Z}
            pauli_probs_multi_unnorm += list(single_probs * bx)

            # Probabilities for Y{I,X,Y,Z}
            pauli_probs_multi_unnorm += list(single_probs * by)

            # Probabilities for Z{I,X,Y,Z}
            pauli_probs_multi_unnorm += list(single_probs * bz)

            # Normalize Weights
            total = sum(pauli_probs_multi_unnorm)

            pauli_probs_multi = [
                weights * (bias_prob / total)
                for weights in pauli_probs_multi_unnorm
            ]

        else:
            pauli_probs_single = [0] * 3
            pauli_probs_multi = [0] * 15

        return pauli_probs_single, pauli_probs_multi

    @abstractmethod
    def _append_clifford_noise(self, 
                               out: stim.Circuit, 
                               instruction: stim.CircuitInstruction,
                               ) -> None:
        pass

    @abstractmethod
    def _append_before_measurement_noise(self, 
                                  out: stim.Circuit, 
                                  instruction: stim.CircuitInstruction,
                                  ) -> None:
        pass

    @abstractmethod
    def _append_after_measurement_noise(self, 
                                  out: stim.Circuit, 
                                  instruction: stim.CircuitInstruction,
                                  ) -> None:
        pass

    @abstractmethod
    def _append_reset_noise(self, 
                            out: stim.Circuit, 
                            instruction: stim.CircuitInstruction,
                            ) -> None:
        pass

    @abstractmethod
    def _append_before_round_depol(self, 
                                   out: stim.Circuit, 
                                   instruction: stim.CircuitInstruction,
                                   ) -> None:
        pass

    @abstractmethod
    def _append_idling_error(self, 
                            out: stim.Circuit, 
                            qubits_idx: list[int],
                            waiting_for_r_m: bool = False
                            ) -> None:
        pass

class CircuitNoise(NoiseModel):

    def __init__(self, *args, **kwargs):

        """
        Implementation of the SI1000 noise model used by Gidney in arxiv:2302.07395
            -> 1 Gate Cliffords: have Depol after implementation (p/10)
            -> 2 Gate Cliffords: have Depol2 after miplementation (p)
            -> After Reset: X_error (2p)
            -> After measurement flip (5p) and Depol1(p)
            -> Qubits nots Measured or Reset during Rounds which include M or R on other
                qubits -> Depol1(2p)

        Args:
            circuit (stim.Circuit): The input quantum circuit to which noise will be added.
            noise (float): Specify Probability P
            distance: Distance of the code patch to determine which measurement are the last
            ft_init: Boolean to indicate if the noise should be applied at the beginning of the
                     circuit (If not ft init, noise is not added to the first round of ancilla
                     measurement)
            ft_measurements: Boolean to indicate if the noise should be applied at the end of the
                             circuit. If not ft measurements, noise is not added to the final
                             measurement readout of the data qubits
        """

        # Set up Preliminary Information
        super().__init__(*args, **kwargs)

        # Get Probabilities out of Dict
        prob_circ = self.noise.get("CircuitNoiseProbability")
        if prob_circ * 5 > 1:
            raise ValueError("Probability is to high, Measurement errors are 5 * p!")
        
        # Set Preliminary NoiseModel Information
        self.after_c_depol1_prob = prob_circ / 10
        self.after_c_depol2_prob = prob_circ
        self.after_r_flip_prob = 2 * prob_circ
        self.after_m_flip_prob = 5 * prob_circ
        self.after_m_depol_prob = prob_circ
        self.wait_m_r_prob = 2 * prob_circ
        self.skip_noise = True if prob_circ == 0 else False

        # Bias channels
        self.pauli_probs_multi      = None
        self.pauli_probs_single_p10 = None
        self.pauli_probs_single_p   = None
        self.pauli_probs_single_2p  = None

        if self.bias is not None:
            self.pauli_probs_single_p10, _ = self._create_bias_noise_model(self.bias, prob_circ / 10)
            self.pauli_probs_single_p, self.pauli_probs_multi = self._create_bias_noise_model(self.bias, prob_circ)
            self.pauli_probs_single_2p, _ = self._create_bias_noise_model(self.bias, 2 * prob_circ)

    def _append_clifford_noise(
        self,
        out: stim.Circuit,
        instruction: stim.CircuitInstruction,
    ) -> None:
        # Return without noise if no Clifford noise is present
        if self.skip_noise:
            return

        targets = instruction.targets_copy()
        is_two_qubit = instruction.name in {"CX", "CZ", "XCY"}

        if is_two_qubit:
            # Check if Multi-Qubit gate has record targets
            # -> Skip complelty as this needs to be handled as single qubit gate
            if any(t.is_measurement_record_target for t in targets):
                """
                It really makes no sense to have any noise here as there correction would be 
                implemented as a pauli fram correction and therefore classically tracked. 
                So there shouldn't be any noise here...!
                """
                pass

            else:
                qubits = [t.value for t in targets]
                # Check if Biased or normal Noise Channel needs to be used
                if self.pauli_probs_multi is None:
                    out.append("DEPOLARIZE2", qubits, self.after_c_depol2_prob)
                else:
                    out.append("PAULI_CHANNEL_2", qubits, self.pauli_probs_multi)
        else:
            # Single qubit gate
            qubits = [t.value for t in targets]
            # Check if Biased or normal Noise Channel needs to be used
            if self.pauli_probs_single_p10 is None:
                out.append("DEPOLARIZE1", qubits, self.after_c_depol1_prob)
            else:
                out.append("PAULI_CHANNEL_1", qubits, self.pauli_probs_single_p10)

    def _append_before_measurement_noise(
        self,
        out: stim.Circuit,
        instruction: stim.CircuitInstruction,
    ) -> None:
        # Return without noise if no measurement noise is present
        if self.skip_noise:
            return
        
        # Add Noise -> Adding X flip to simulate measruement flip
        for t in instruction.targets_copy():
            # Adding before Measurement Flip
            out.append("X_ERROR", [t.value], self.after_m_flip_prob)

    def _append_after_measurement_noise(
        self,
        out: stim.Circuit,
        instruction: stim.CircuitInstruction,
    ) -> None:
        # Return without noise if no measurement noise is present
        if self.skip_noise:
            return
        
        # Add Noise -> Depolarize after Measurement
        for t in instruction.targets_copy():
            # Check if Biased or normal Noise Channel needs to be used
            if self.pauli_probs_single_p is None:
                out.append("DEPOLARIZE1", [t.value], self.after_m_depol_prob)
            else:
                out.append("PAULI_CHANNEL_1", [t.value], self.pauli_probs_single_p)

    def _append_reset_noise(
        self,
        out: stim.Circuit,
        instruction: stim.CircuitInstruction,
    ) -> None:
        # Return without noise if no reset noise is present
        if self.skip_noise:
            return
        # Add Noise
        for t in instruction.targets_copy():
            # Adding after reset Flip
            out.append("X_ERROR", [t.value], self.after_r_flip_prob)
            
    def _append_idling_error(self, 
                             out: stim.Circuit, 
                             qubits_idx: list[int],
                             waiting_for_r_m: bool = False
                             ) -> None:
        
        # Return without noise if no before round depol noise is present
        if self.skip_noise:
            return
        
        # Check if normal idling error or additional wait error while measurement
        # or resets are perofrmed
        if waiting_for_r_m:
            # Check if Biased or normal Noise Channel needs to be used
            if self.pauli_probs_single_2p is None:
                out.append("DEPOLARIZE1", qubits_idx, self.wait_m_r_prob)
            else:
                out.append("PAULI_CHANNEL_1", qubits_idx, self.pauli_probs_single_2p)
        else:
            # Check if Biased or normal Noise Channel needs to be used
            if self.pauli_probs_single_p10 is None:
                out.append("DEPOLARIZE1", qubits_idx, self.after_c_depol1_prob)
            else:
                out.append("PAULI_CHANNEL_1", qubits_idx, self.pauli_probs_single_p10)

    def _append_before_round_depol(self, out, instruction):
        pass
                
class PenomenologicalNoise(NoiseModel):

    def __init__(self, *args, **kwargs):
        """
        Adding phenomological noise to a given Circuit
            -> Before Round Depolarization
            -> Before Measurement flip!

        Args:
            circuit (stim.Circuit): The input quantum circuit to which noise will be added.
            noise (dict): A dictionary specifying the noise parameters. Expected keys are:
                - phenom_prob: Probability of x error happening before measurement 
                & Probability of depolarizing noise before each round on data qubits
            distance: Distance of the code patch to determine which measurement are the last
            ft_init: Boolean to indicate if the noise should be applied at the beginning of the
                    circuit (If not ft init, noise is not added to the first round of ancilla
                    measurement)
            ft_measurements: Boolean to indicate if the noise should be applied at the end of the
                            circuit. If not ft measurements, noise is not added to the final
                            measurement readout of the data qubits
        """

        # Set up Preliminary Information
        super().__init__(*args, **kwargs)

        # Set Preliminary NoiseModel Information
        self.phenom_prob = self.noise.get("PhemoNoiseProbability", 0)
        self.skip_noise = True if self.phenom_prob == 0 else False

        # Bias Noise
        self.pauli_probs_single = None

        if self.bias is not None:
            self.pauli_probs_single, _ = self._create_bias_noise_model(
                self.bias, self.phenom_prob
            )

    def _append_before_measurement_noise(
        self,
        out: stim.Circuit,
        instruction: stim.CircuitInstruction,
    ) -> None:
        # Return without noise if no measurement noise is present
        if self.skip_noise:
            return
        # Add Noise -> Measurement Flip
        for t in instruction.targets_copy():
            # Adding before Measurement Flip
            out.append("X_ERROR", [t.value], self.phenom_prob)

    def _append_before_round_depol(
        self,
        out: stim.Circuit,
        instruction: stim.CircuitInstruction,
    ) -> None:
        
        """
        THIS DOES CURRENTLY NOT WORK AS ROUND TRACKER IS NOT IMPLEMENTED!
        -> Maybe add in future version but for thesis not needed!
        """

        # Return without noise if no before round depol noise is present
        if self.skip_noise:
            return
        # Add Noise -> Either Normal or Biased
        if self.pauli_probs_single is None:
            out.append("DEPOLARIZE1", [t.value for t in instruction.targets_copy() 
                                if t.value in self.data_qubits], self.phenom_prob)
        else:
            out.append("PAULI_CHANNEL_1", [t.value for t in instruction.targets_copy() 
                                if t.value in self.data_qubits], self.pauli_probs_single)
        
    def _append_clifford_noise(self, out, instruction):
        pass
    
    def _append_reset_noise(self, out, instruction):
        pass

    def _append_idling_error(self, out, qubits_idx, waiting_for_r_m):
        pass

    def _append_after_measurement_noise(self, out, instruction):
        pass