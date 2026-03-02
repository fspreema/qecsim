import numpy as np
import stim


class CircuitNoise:

    FIRST_NOISY_RESET = 2
    CLIFFORD_OPERATIONS = ["H", "CX", "S", "S_DAG", "CZ", "XCY", "SQRT_X_DAG"]
    MEASUREMENT_OPERATIONS = ["M", "MX", "MY"]
    RESET_OPERATIONS = ["R", "RX", "RY", "RZ"]

    def __init__(self,
                 circuit: stim.Circuit,
                 noise: dict,
                 distance: int,
                 ft_init: bool = True,
                 ft_measurements: bool = True):

        """
        Adding Circuit type noise to a given Circuit
            -> Clifford gates have Depol before implementation
            -> Before measurement flips
            -> After reset Flips

        Args:
            circuit (stim.Circuit): The input quantum circuit to which noise will be added.
            noise (dict): A dictionary specifying the noise parameters. Expected keys are:
                - before_m_flip_prob: Probability of x error happening before measurement
                - after_r_flip: Probability of x error after reset
                - after_c_depol_prob: Probability of depolarizing noise for Clifford gates
                - before_round_depol: Probability of depolarizing noise before each round
            distance: Distance of the code patch to determine which measurement are the last
            ft_init: Boolean to indicate if the noise should be applied at the beginning of the
                     circuit (If not ft init, noise is not added to the first round of ancilla
                     measurement)
            ft_measurements: Boolean to indicate if the noise should be applied at the end of the
                             circuit. If not ft measurements, noise is not added to the final
                             measurement readout of the data qubits
        """

        # Set up Preliminary Information
        self.circuit = circuit
        self.noise = noise
        self.ft_init = ft_init
        self.ft_measurements = ft_measurements
        self.distance = distance
        self.tracked_reset_idx = None
        self.curr_reset_num = 0

        # Get Probabilities out of Dict
        self.before_m_flip_prob = self.noise.get("before_m_flip_prob", 0)
        self.after_r_flip_prob = self.noise.get("after_r_flip", 0)
        self.after_c_depol_prob = self.noise.get("after_c_depol_prob", 0)
        self.before_round_depol = self.noise.get("before_round_depol", 0)

    def apply(self) -> stim.Circuit:

        # Create Noisy Circuit
        noisy_circuit = self._apply_noise_operations(circuit=self.circuit)

        return noisy_circuit

    def _apply_noise_operations(self, circuit: stim.Circuit) -> stim.Circuit:
        # Initilizing noisy circuit
        noisy_circuit = stim.Circuit()

        for instruction in circuit:
            # Filter Out Instructions
            if isinstance(instruction, stim.CircuitInstruction):
                # Check what type of Operation we have
                if isinstance(instruction, stim.CircuitInstruction):
                    if instruction.name == "MR":
                        raise NotImplementedError("MR Operations are not supported")

                    if instruction.name in self.CLIFFORD_OPERATIONS:
                        noisy_circuit.append(instruction)
                        # Add Noise
                        if self._should_apply_noise(instruction=instruction):
                            self._append_clifford_noise(noisy_circuit, instruction)

                    elif instruction.name in self.MEASUREMENT_OPERATIONS:
                        # Add Noise
                        if self._should_apply_noise(instruction=instruction):
                            self._append_measurement_noise(noisy_circuit, instruction)
                        noisy_circuit.append(instruction)

                    elif instruction.name in self.RESET_OPERATIONS:
                        noisy_circuit.append(instruction)
                        # Track Resets
                        self._track_resets(instruction)
                        # Add Noise
                        if self._should_apply_noise(instruction=instruction):
                            self._append_reset_noise(noisy_circuit, instruction)

                    else:
                        noisy_circuit.append(instruction)

            # Filter Out Repeat-Blocks
            elif isinstance(instruction, stim.CircuitRepeatBlock):
                inner_noisy = self._apply_noise_operations(instruction.body_copy())
                noisy_circuit += inner_noisy * instruction.repeat_count

            else:
                noisy_circuit.append(instruction)

        return noisy_circuit

    def _should_apply_noise(self, instruction: stim.CircuitInstruction) -> bool:
        """
        Noise gating rule:
        - If ft_measurements is False, do not add noise to the final small measurement cluster.
        - If ft_init is enabled, apply noise from the start.
        - Otherwise, only apply noise once the first "noisy reset" threshold is reached.
        """
        # If ft_measurements is False, do not add noise to the final small measurement cluster.
        if (
                not self.ft_measurements
                and instruction.name in self.MEASUREMENT_OPERATIONS
                and len(instruction.targets_copy()) <= self.distance
        ):
            return False

        # Add noise everywhere if ft_init is true
        if self.ft_init:
            return True

        # Add Noise if first noisy Reset is reached and ft_init is false
        if self.curr_reset_num >= self.FIRST_NOISY_RESET:
            return True

        return False

    def _track_resets(self, instruction: stim.CircuitInstruction) -> None:

        # Check if anything is already tracked:
        if self.tracked_reset_idx is None:

            # Add Reset to Tracker if needed
            if instruction.name in self.RESET_OPERATIONS :
                self.curr_reset_num += 1
                self.tracked_reset_idx = instruction.targets_copy()[0].value

        # If already tracked, check if the index is the same, if so add number of resets
        else:
            if instruction.targets_copy()[0].value == self.tracked_reset_idx:
                self.curr_reset_num += 1

    def _append_clifford_noise(
        self,
        out: stim.Circuit,
        instruction: stim.CircuitInstruction,
    ) -> None:
        # Return without noise if no Clifford noise is present
        if self.after_c_depol_prob <= 0:
            return

        targets = instruction.targets_copy()
        is_two_qubit = instruction.name in {"CX", "CZ", "XCY"}

        if is_two_qubit:
            # Check if Multi-Qubit gate has record targets
            # -> Skip complelty as this needs to be handled as single qubit gate
            if any(t.is_measurement_record_target for t in targets):
                # Single qubit Depolarize for rec dependent targets
                qubits = [t.value for t in targets if t.is_qubit_target]
                out.append("DEPOLARIZE1", qubits, self.after_c_depol_prob)
            else:
                qubits = [t.value for t in targets]
                out.append("DEPOLARIZE2", qubits, self.after_c_depol_prob)
        else:
            # Single qubit gate
            qubits = [t.value for t in targets]
            out.append("DEPOLARIZE1", qubits, self.after_c_depol_prob)

    def _append_measurement_noise(
        self,
        out: stim.Circuit,
        instruction: stim.CircuitInstruction,
    ) -> None:
        # Return without noise if no measurement noise is present
        if self.before_m_flip_prob <= 0:
            return
        # Add Noise
        for t in instruction.targets_copy():
            # Adding before Measurement Flip
            out.append("X_ERROR", [t.value], self.before_m_flip_prob)

    def _append_reset_noise(
        self,
        out: stim.Circuit,
        instruction: stim.CircuitInstruction,
    ) -> None:
        # Return without noise if no reset noise is present
        if self.after_r_flip_prob <= 0:
            return
        # Add Noise
        for t in instruction.targets_copy():
            # Adding after reset Flip
            out.append("X_ERROR", [t.value], self.after_r_flip_prob)

class CircuitNoiseYBasis:
    """
    Adds Circuit Noise analog to CircuitNoise but in Y-Basis
    -> Initialization etc. is not included in this noise model
    -> Only tick window inside the memory experiment is considered
    """

    def __init__(self, circuit: stim.Circuit, noise: dict):
        self.circuit = circuit
        self.noise = noise
        self.noise_before_operators: list[str] = ["H", "CX", "S", "S_DAG", "CZ"]


class BiasNoise:

    FIRST_NOISY_RESET = 2
    CLIFFORD_OPERATIONS = ["H", "CX", "S", "S_DAG", "CZ"]

    """
    Adds Bias Noise to a given Circuit

    Args:
        circuit (stim.Circuit): The input quantum circuit to which noise will be added.
        noise (dict): A dictionary specifying the noise parameters. Expected keys are:
            - bias: List of bias values [b_x, b_y, b_z] for Pauli channels
            - after_c_custom_noise: The complete Probability of An error happening after Clifford
                                    gates
                -> This is the probability which gets split up depending on the bias values

    Returns:
        stim.Circuit: The noisy quantum circuit with bias noise applied.
    """

    def __init__(self, circuit: stim.Circuit, noise: dict):
        self.after_c_p_xyz_multi = None
        self.after_c_p_xyz = None
        self.circuit = circuit
        self.noise = noise
        self.bias = noise.get("bias", [0, 0, 0])
        self.after_c_custom_noise = noise.get("after_c_custom_noise", 0)

    def _create_bias_noise_model(self):
        """
        Creates the noise needed fot the Curstom Pauli Channels
        -> Given a bias list and full chance of a physical_error
        -> type of bias is [b_x, b_y, b_z] with b_x + b_y + b_z = 1
        -> Returns:
            - List of 3 probabilities for PAULI_CHANNEL_1
            - List of 15 probabilities for PAULI_CHANNEL_2
        """

        # Check if Prob under 3/4 else BLoch sphere turned inside out
        if self.after_c_custom_noise > 3 / 4:
            raise ValueError("Probability for custom Pauli channel too high (> 3/4)")

        if np.any(self.bias):
            #######################################
            # Adding Noise for single Pauli Channel
            #######################################

            after_c_p_xyz = [
                (self.after_c_custom_noise / sum(self.bias)) * self.bias[i] for i in range(3)
            ]

            ######################################
            # Adding Noise for multi Pauli Channel
            ######################################

            bx, by, bz = self.bias
            single_probs = np.array([1, bx, by, bz])

            after_c_p_xyz_multi_unnorm: list = []

            # Probabilities for I{I,X,Y,Z}
            after_c_p_xyz_multi_unnorm += list(single_probs)

            # Remove II prob.
            after_c_p_xyz_multi_unnorm.pop(0)

            # Probabilities for X{I,X,Y,Z}
            after_c_p_xyz_multi_unnorm += list(single_probs * bx)

            # Probabilities for Y{I,X,Y,Z}
            after_c_p_xyz_multi_unnorm += list(single_probs * by)

            # Probabilities for Z{I,X,Y,Z}
            after_c_p_xyz_multi_unnorm += list(single_probs * bz)

            # Normalize Weights
            total = sum(after_c_p_xyz_multi_unnorm)

            after_c_p_xyz_multi = [
                weights * (self.after_c_custom_noise / total)
                for weights in after_c_p_xyz_multi_unnorm
            ]

        else:
            after_c_p_xyz = [0] * 3
            after_c_p_xyz_multi = [0] * 15

        return after_c_p_xyz, after_c_p_xyz_multi

    def apply(self) -> stim.Circuit:
        # Get Bias noise Values
        self.after_c_p_xyz, self.after_c_p_xyz_multi = self._create_bias_noise_model()

        # Create Noisy Circuit
        noisy_circuit = self._apply_noise_operations(circuit=self.circuit)

        return noisy_circuit

    def _apply_noise_operations(self, circuit: stim.Circuit) -> stim.Circuit:

        # Initilizing noisy circuit
        noisy_circuit = stim.Circuit()

        for instructions in circuit:
            # Filter Out Instructions
            if isinstance(instructions, stim.CircuitInstruction):
                # Check what type of Operation we have
                if instructions.name in self.CLIFFORD_OPERATIONS:
                    # Adding Clifford
                    noisy_circuit.append(instructions)

                    # Adding Depolarize Noise after Clifford
                    self._add_clifford_noise(noisy_circuit, instructions)

                else:
                    # For all other instructions, just append them without noise
                    noisy_circuit.append(instructions)

            # Filter Out Repeat-Blocks
            elif isinstance(instructions, stim.CircuitRepeatBlock):
                # Initlize noisy repeat Block
                noisy_repeat = stim.Circuit()

                # Deconstruct Repeat Body into Instructions
                inner_circuit = instructions.body_copy()
                repeat_count = instructions.repeat_count

                # Apply Noise to Inner Circuit
                noisy_repeat += self._apply_noise_operations(inner_circuit)

                # Append Repeat Block to Noisy Circuit
                noisy_circuit += noisy_repeat * repeat_count

            else:
                # For any other type of instruction, just append it
                noisy_circuit.append(instructions)

        return noisy_circuit

    def _add_clifford_noise(self,
                            out: stim.Circuit,
                            instruction: stim.CircuitInstruction) -> None:

        is_two_qubit = instruction.name in {"CX", "CZ"}

        # Adding Depolarize Noise after Clifford
        if is_two_qubit:
            # Check if 15 MPP is non zer0
            if np.any(self.after_c_p_xyz_multi):
                # Multi-qubit gate
                qubits = [targets.value for targets in instruction.targets_copy()]
                out.append(
                    "PAULI_CHANNEL_2",
                    qubits,
                    self.after_c_p_xyz_multi,
                )

        else:
            # Check if single Value is non zero
            if np.any(self.after_c_p_xyz):
                # Single qubit gate
                qubit = [targets.value for targets in instruction.targets_copy()]
                out.append(
                    "PAULI_CHANNEL_1",
                    qubit,
                    self.after_c_p_xyz,
                )
