import numpy as np
import stim


class CircuitNoise:
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
    """

    def __init__(self, circuit: stim.Circuit, noise: dict):
        self.circuit = circuit
        self.noise = noise
        self.noise_before_operators: list[str] = [
            "H",
            "CX",
            "S",
            "S_DAG",
            "CZ",
            "XCY",
            "SQRT_X_DAG",
        ]
        self.measurement_noise_operators: list[str] = ["M", "MX", "MY"]
        self.reset_noise_operators: list[str] = ["R", "RX", "RY", "RZ"]

    def apply(self) -> stim.Circuit:
        # Initilizing noisy circuit
        noisy_circuit = stim.Circuit()

        for instructions in self.circuit:
            # Filter Out Instructions
            if isinstance(instructions, stim.CircuitInstruction):
                # Check what type of Operation we have
                if instructions.name in self.noise_before_operators:
                    # Adding Clifford
                    noisy_circuit.append(instructions)

                    # Check if Dict has a non zero value
                    if self.noise.get("after_c_depol_prob", 0) > 0:
                        # Adding Depolarize Noise after Clifford
                        if instructions.name in {"CX", "CZ", "XCY"}:
                            # Check if Multi-Qubit gate has record targets
                            # -> Skip complelty as this needs to be handled as single qubit gate
                            if any(
                                target.is_measurement_record_target
                                for target in instructions.targets_copy()
                            ):
                                # Single qubit Depolarize for rec dependent targets
                                qubits = [
                                    targets.value
                                    for targets in instructions.targets_copy()
                                    if targets.is_qubit_target
                                ]

                                noisy_circuit.append(
                                    "DEPOLARIZE1",
                                    qubits,
                                    self.noise.get("after_c_depol_prob", 0),
                                )

                            else:
                                # Multi-qubit gate
                                qubits = [targets.value for targets in instructions.targets_copy()]
                                noisy_circuit.append(
                                    "DEPOLARIZE2",
                                    qubits,
                                    self.noise.get("after_c_depol_prob", 0),
                                )

                        else:
                            # Single qubit gate
                            qubit = [targets.value for targets in instructions.targets_copy()]
                            noisy_circuit.append(
                                "DEPOLARIZE1",
                                qubit,
                                self.noise.get("after_c_depol_prob", 0),
                            )

                elif instructions.name in self.measurement_noise_operators:
                    # Check if Dict has a non zero value
                    if self.noise.get("before_m_flip_prob", 0) > 0:
                        # Adding Before Measurement Flip
                        for targets in instructions.targets_copy():
                            qubit = targets.value
                            noisy_circuit.append(
                                "X_ERROR",
                                [qubit],
                                self.noise.get("before_m_flip_prob", 0),
                            )

                    # Adding Measurement
                    noisy_circuit.append(instructions)

                elif instructions.name in self.reset_noise_operators:
                    # Adding Reset
                    noisy_circuit.append(instructions)

                    # Check if Dict has a non zero value
                    if self.noise.get("after_r_flip", 0) > 0:
                        # Adding after reset noise
                        for targets in instructions.targets_copy():
                            qubit = targets.value
                            noisy_circuit.append(
                                "X_ERROR",
                                [qubit],
                                self.noise.get("after_r_flip", 0),
                            )

                else:
                    # For all other instructions, just append them without noise
                    noisy_circuit.append(instructions)

            # Check for unsupported MR operations
            elif isinstance(instructions, stim.CircuitInstruction) and instructions.name == "MR":
                raise NotImplementedError("MR Operations are not supported")

            # Filter Out Repeat-Blocks
            elif isinstance(instructions, stim.CircuitRepeatBlock):
                # Initlize noisy repeat Block
                noisy_repeat = stim.Circuit()

                # Deconstruct Repeat Body into Instructions
                inner_circuit = instructions.body_copy()
                repeat_count = instructions.repeat_count

                for inner_instructions in inner_circuit:
                    # Check for unsupported MR operations in repeat blocks
                    if inner_instructions.name == "MR":
                        raise NotImplementedError("MR Operations are not supported")

                    # Check what type of Operation we have
                    if inner_instructions.name in self.noise_before_operators:
                        # Adding Clifford
                        noisy_repeat.append(inner_instructions)

                        # Check if Dict has a non zero value
                        if self.noise.get("after_c_depol_prob", 0) > 0:
                            # Adding Depolarize Noise after Clifford
                            if inner_instructions.name in {"CX", "CZ", "XCY"}:
                                # Check if Multi-Qubit gate has record targets
                                # -> Skip complelty as this needs to be handled as single qubit gate
                                if any(
                                    target.is_measurement_record_target
                                    for target in inner_instructions.targets_copy()
                                ):
                                    # Single qubit Depolarize for rec dependent targets
                                    qubits = [
                                        targets.value
                                        for targets in inner_instructions.targets_copy()
                                        if targets.is_qubit_target
                                    ]

                                    noisy_circuit.append(
                                        "DEPOLARIZE1",
                                        qubits,
                                        self.noise.get("after_c_depol_prob", 0),
                                    )

                                else:
                                    # Multi-qubit gate
                                    qubits = [
                                        targets.value for targets in inner_instructions.targets_copy()
                                    ]
                                    noisy_repeat.append(
                                        "DEPOLARIZE2",
                                        qubits,
                                        self.noise.get("after_c_depol_prob", 0),
                                    )

                            else:
                                # Single qubit gate
                                qubit = [targets.value for targets in inner_instructions.targets_copy()]
                                noisy_repeat.append(
                                    "DEPOLARIZE1",
                                    qubit,
                                    self.noise.get("after_c_depol_prob", 0),
                                )

                    elif inner_instructions.name in self.measurement_noise_operators:
                        # Check if Dict has a non zero value
                        if self.noise.get("before_m_flip_prob", 0) > 0:
                            # Adding Before Measurement Flip
                            for targets in inner_instructions.targets_copy():
                                qubit = targets.value
                                noisy_repeat.append(
                                    "X_ERROR",
                                    qubit,
                                    self.noise.get("before_m_flip_prob", 0),
                                )

                        # Adding Measurement
                        noisy_repeat.append(inner_instructions)

                    elif inner_instructions.name in self.reset_noise_operators:
                        # Adding Reset
                        noisy_repeat.append(inner_instructions)

                        # Check if Dict has a non zero value
                        if self.noise.get("after_r_flip", 0) > 0:
                            # Adding after reset noise
                            for targets in inner_instructions.targets_copy():
                                qubit = targets.value
                                noisy_repeat.append(
                                    "X_ERROR",
                                    [qubit],
                                    self.noise.get("after_r_flip", 0),
                                )

                    else:
                        # For all other instructions in repeat block, just append them without noise
                        noisy_repeat.append(inner_instructions)

                # Append Repeat Block to Noisy Circuit
                noisy_circuit += noisy_repeat * repeat_count

            else:
                # For any other type of instruction, just append it
                noisy_circuit.append(instructions)

        return noisy_circuit


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
        self.circuit = circuit
        self.noise = noise
        self.bias = noise.get("bias", [0, 0, 0])
        self.after_c_custom_noise = noise.get("after_c_custom_noise", 0)
        self.noise_before_operators = ["H", "CX", "S", "S_DAG", "CZ"]

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
        after_c_p_xyz, after_c_p_xyz_multi = self._create_bias_noise_model()

        # Initilizing noisy circuit
        noisy_circuit = stim.Circuit()

        for instructions in self.circuit:
            # Filter Out Instructions
            if isinstance(instructions, stim.CircuitInstruction):
                # Check what type of Operation we have
                if instructions.name in self.noise_before_operators:
                    # Adding Clifford
                    noisy_circuit.append(instructions)

                    # Adding Depolarize Noise after Clifford
                    if instructions.name in {"CX", "CZ"}:
                        # Check if 15 MPP is non zer0
                        if np.any(after_c_p_xyz_multi):
                            # Multi-qubit gate
                            qubits = [targets.value for targets in instructions.targets_copy()]
                            noisy_circuit.append(
                                "PAULI_CHANNEL_2",
                                qubits,
                                after_c_p_xyz_multi,
                            )

                    elif instructions.name not in {"CX", "CZ"}:
                        # Check if single Value is non zero
                        if np.any(after_c_p_xyz):
                            # Single qubit gate
                            qubit = [targets.value for targets in instructions.targets_copy()]
                            noisy_circuit.append(
                                "PAULI_CHANNEL_1",
                                qubit,
                                after_c_p_xyz,
                            )

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

                for inner_instructions in inner_circuit:
                    # Check what type of Operation we have
                    if inner_instructions.name in self.noise_before_operators:
                        # Adding Clifford
                        noisy_repeat.append(instructions)

                        # Adding Depolarize Noise after Clifford
                        if inner_instructions.name in {"CX", "CZ"}:
                            # Check if 15 MPP is non zer0
                            if np.any(after_c_p_xyz_multi):
                                # Multi-qubit gate
                                qubits = [targets.value for targets in inner_instructions.targets_copy()]
                                noisy_repeat.append(
                                    "PAULI_CHANNEL_2",
                                    qubits,
                                    after_c_p_xyz_multi,
                                )

                        elif inner_instructions.name not in {"CX", "CZ"}:
                            # Check if single Value is non zero
                            if np.any(after_c_p_xyz):
                                # Single qubit gate
                                qubit = [targets.value for targets in inner_instructions.targets_copy()]
                                noisy_repeat.append(
                                    "PAULI_CHANNEL_1",
                                    qubit,
                                    after_c_p_xyz,
                                )

                    else:
                        # For all other instructions in repeat block, just append them without noise
                        noisy_repeat.append(inner_instructions)

                # Append Repeat Block to Noisy Circuit
                noisy_circuit += noisy_repeat * repeat_count

            else:
                # For any other type of instruction, just append it
                noisy_circuit.append(instructions)

        return noisy_circuit
