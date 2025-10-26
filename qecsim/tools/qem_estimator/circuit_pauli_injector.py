import numpy as np
import stim

__all__ = ["pauli_injector"]


def pauli_injector(
    circuit: stim.Circuit,
    noise_info: np.array,
    batch_size: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    This function uses stim.FlipSimulator in order to determine which measurements get flipped
    after an arbitary pauli error is included into the circuit

    Procedure:
    -> We take the infromation out of the noise_info array
    -> We run through each of the postitions and ticks and place either Identity or XYZ pauli
       depending on a probability of the individual channel
        -> Prob is created by the weights and is already provided in the array
    -> If the Circuit has run through and all neccessary changes have been made we give
       out the measurement flips
    -> Each of the measurement flips also include a sign switch which is also
       returned in a seperate array

    Arguments:
        circuit: Currently used stim.Circuit where the pauli should be inserted
        noise_info: Array giving information of what error happens at what Tick and Qubit Postition

    Returns:
    -> Measurement flips in an array
    -> Sign keeper where each index corresponds to the index of the qubit this
       sign needs to be applied to
    """

    #############################
    # Initlizie sign/gamma keeper
    #############################

    """
    We need this as we have also the sign as well as gamma which need to be multiplied by the final
    measurement outcome
    """

    sign_keeper = np.ones((circuit.num_qubits, batch_size), dtype=int)
    gamma_keeper = np.ones((circuit.num_qubits, batch_size), dtype=float)

    ###########################
    # Initlizing Flip Simulator
    ###########################

    flip = stim.FlipSimulator(
        batch_size=batch_size,
        num_qubits=circuit.num_qubits,
        disable_stabilizer_randomization=True,
    )

    ##################################################################
    # Getting Tick Information (All noisy ticks an current noisy tick)
    ##################################################################

    all_noisy_ticks = noise_info[:, 1]
    current_tick = -1

    ########################################################
    # Running through circuit tick by tick and adding paulis
    ########################################################

    for i in range(len(all_noisy_ticks)):
        """
        flip.do adds instruction to every batch instance -> therefore keep outside the batch loop
        """

        # Adding all instructions until Pauli insertion
        for index, operation in enumerate(circuit):
            if current_tick < index < all_noisy_ticks[i]:
                flip.do(operation)

        # Retireving Current Operation, qubit index, Prob
        qubit_index = list(noise_info[i, 2])
        current_err_prob = noise_info[i, 4]
        error_type = noise_info[i, 0]
        sgn_i = int(noise_info[i, 5])
        sgn_err = int(noise_info[i, 6])
        current_gamma = noise_info[i, 7]

        # Checking DEPOL2 -> No looping over the individual qubits, therefore outside the for loop
        if error_type == "DEPOL2":
            #############################################################################
            # Looping over Batch number and simulating batch_nuberm fo different outcomes
            #############################################################################

            for current_batch in range(batch_size):
                # Creating current random number
                current_rdm_nmbr = np.random.rand()

                # Checking if error is applied or not
                if current_rdm_nmbr < current_err_prob:
                    paulis = ["I", "X", "Y", "Z"]
                    depol2 = [(a, b) for a in paulis for b in paulis if not (a == "I" and b == "I")]

                    # Choose Random Multi Qubit Pauli for depol
                    rng = np.random.default_rng()
                    chosen_pauli1, chosen_pauli2 = depol2[rng.integers(15)]

                    # Add Paulis into noiseless circuit
                    flip.set_pauli_flip(
                        chosen_pauli1,
                        qubit_index=qubit_index[0],
                        instance_index=current_batch,
                    )
                    flip.set_pauli_flip(
                        chosen_pauli2,
                        qubit_index=qubit_index[1],
                        instance_index=current_batch,
                    )

                    # Updating sign keeper
                    """
                    We add the sign to one fo the qubits
                    -> As we multiply gamma together in the end it doesnt really matter what 
                    qubit, but as DEPOL2 is one operation we add the gamma only to one qubit!
                    """

                    sign_keeper[qubit_index[0], current_batch] *= sgn_err

                else:
                    # Update sign keeper
                    sign_keeper[qubit_index[0], current_batch] *= sgn_i

                # Add gamma keeper
                gamma_keeper[qubit_index[0], current_batch] *= current_gamma

            # Updating Boundaries for adding instructions to the flip sim
            current_tick = all_noisy_ticks[i] - 1
            continue

        # Looping over all qubit indexes
        for current_qubit in qubit_index:
            # Check what type of noise we have in order to determine what elements the inverse
            # channel has
            if error_type == "DEPOL":
                #############################################################################
                # Looping over Batch number and simulating batch_nuberm fo different outcomes
                #############################################################################

                for current_batch in range(batch_size):
                    # Creating current random number
                    current_rdm_nmbr = np.random.rand()

                    # Checking if error is applied or not
                    if current_rdm_nmbr < current_err_prob:
                        paulis = ["X", "Y", "Z"]

                        rng = np.random.default_rng()
                        chosen_pauli = paulis[rng.integers(3)]

                        # Add Pauli into noiseless circuit
                        flip.set_pauli_flip(
                            chosen_pauli,
                            qubit_index=current_qubit,
                            instance_index=current_batch,
                        )

                        # Updating sign keeper
                        sign_keeper[current_qubit, current_batch] *= sgn_err

                    else:
                        # Update sign keeper
                        sign_keeper[current_qubit, current_batch] *= sgn_i

                    # Add gamma keeper
                    gamma_keeper[current_qubit, current_batch] *= current_gamma

            elif error_type == "X_ERR":
                #############################################################################
                # Looping over Batch number and simulating batch_nuberm fo different outcomes
                #############################################################################

                for current_batch in range(batch_size):
                    # Creating current random number
                    current_rdm_nmbr = np.random.rand()

                    # Checking if error is applied or not
                    if current_rdm_nmbr < current_err_prob:
                        # Here we have only one option besides the identity for the inverse
                        chosen_pauli = "X"

                        # Add Pauli into nosieless Circ
                        flip.set_pauli_flip(
                            chosen_pauli,
                            qubit_index=current_qubit,
                            instance_index=current_batch,
                        )

                        # Updating sign keeper
                        sign_keeper[current_qubit, current_batch] *= sgn_err

                    else:
                        # Update sign keeper
                        sign_keeper[current_qubit, current_batch] *= sgn_i

                    # Add gamma keeper
                    gamma_keeper[current_qubit, current_batch] *= current_gamma

        # Updating Boundaries for adding instructions to the flip sim
        current_tick = all_noisy_ticks[i] - 1

    ########################
    # Adding rest of circuit
    ########################

    for index, operation in enumerate(circuit):
        if index > current_tick:
            flip.do(operation)

    ##################################################
    # Retrieving information about the changed results
    ##################################################

    meas_mask = flip.get_measurement_flips()
    det_mask = flip.get_detector_flips()
    obs_mask = flip.get_observable_flips()

    ########################################
    # Returning if emasurements were flipped
    ########################################

    return (
        meas_mask.astype(bool),
        det_mask.astype(bool),
        obs_mask.astype(bool),
        sign_keeper,
        gamma_keeper,
    )
