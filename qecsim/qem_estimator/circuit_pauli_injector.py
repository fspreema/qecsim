import stim
import numpy as np

__all__ = ["pauli_injector"]

def pauli_injector(circuit : stim.Circuit, 
                           noise_info : np.array) -> tuple[np.array, np.array, np.array, np.array]:

    '''
    This function uses stim.FlipSimualtor in order to determine which measurements get flipped after an arbitary pauli error is included into the circuit
    
    Procedure:
    -> We take the infromation out of the noise_info array
    -> We run through each of the postitions and ticks and place either Identity or XYZ pauli depending on a probability of the individual cahnnel
        -> Prob is created by the weights and is already provided in the array
    -> If the Circuit has run through and all neccessary changes have been made we give out the measurement flips
    -> Each of the measurement flips also include a sign switch which is also returned in a seperate array
    
    Arguments:
        circuit: Currently used stim.Circuit where the pauli should be inserted
        noise_info: Array giving information of what error happens at what Tick and Qubit Postition

    Returns:
    -> Measurement flips in an array
    -> Sign keeper where each index corresponds to the index of the qubit this sign needs to be applied to
    '''

    #############################
    # Initlizie sign/gamma keeper
    #############################

    """
    We need this as we have also the sign as well as gamma which need to be multiplied by the final measurement outcome
    """

    sign_keeper = np.ones(circuit.num_qubits, dtype=int)
    gamma_keeper = np.ones(circuit.num_qubits, dtype=float)

    ###########################
    # Initlizing Flip Simulator
    ###########################
    
    flip = stim.FlipSimulator(
        batch_size = 1,
        num_qubits = circuit.num_qubits,
        disable_stabilizer_randomization = True
        )
    
    ##################################################################
    # Getting Tick Information (All noisy ticks an current noisy tick)
    ##################################################################

    all_noisy_ticks = noise_info[:,1]
    current_tick = -1
    
    ########################################################
    # Running through circuit tick by tick and adding paulis
    ######################################################## 

    for i in range(len(all_noisy_ticks)):

        # Adding all instructions until Pauli insertion
        for index, operation in enumerate(circuit):
            if current_tick < index < all_noisy_ticks[i]:
                flip.do(operation)

        # Retireving Current Operation, qubit index, Prob
        qubit_index = list(noise_info[i,2])
        current_err_prob = noise_info[i,4]
        error_type = noise_info[i,0]
        sgn_I   = int(noise_info[i,5])
        sgn_err = int(noise_info[i,6])
        current_gamma = noise_info[i,7]

        # Looping over all qubit indexes
        for current_qubit in qubit_index:

            # Check what type of noise we have in order to determine what elements the inverse channel has
            if error_type == "DEPOL":

                # Creating current random number
                current_rdm_nmbr = np.random.rand()
            
                # Checking if error is applied or not
                if current_rdm_nmbr < current_err_prob:

                    # Choose Pauli for depol
                    pauli_random = np.random.rand()

                    if pauli_random < 1/3:
                        chosen_pauli = "X"

                    elif pauli_random < 2/3:
                        chosen_pauli = "Y"

                    else:
                        chosen_pauli = "Z"

                    # Add Pauli into noiseless circuit
                    flip.set_pauli_flip(chosen_pauli, qubit_index = current_qubit, instance_index = 0)

                    # Updating sign keeper
                    sign_keeper[current_qubit] *= sgn_err

                else:

                    # Update sign keeper
                    sign_keeper[current_qubit] *= sgn_I

                # Add gamma keeper
                gamma_keeper[current_qubit] *= current_gamma
                    
            
            elif error_type == "X_ERR":
                    
                # Creating current random number
                current_rdm_nmbr = np.random.rand()
            
                # Checking if error is applied or not
                if current_rdm_nmbr < current_err_prob:

                    # Here we have only one option besides the identity for the inverse
                    chosen_pauli = "X"

                    # Add Pauli into nosieless Circ
                    flip.set_pauli_flip(chosen_pauli, qubit_index = current_qubit, instance_index = 0)

                    # Updating sign keeper
                    sign_keeper[current_qubit] *= sgn_err

                else:

                    # Update sign keeper
                    sign_keeper[current_qubit] *= sgn_I

                # Add gamma keeper
                gamma_keeper[current_qubit] *= current_gamma

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
    det_mask  = flip.get_detector_flips()
    obs_mask  = flip.get_observable_flips()

    ########################################
    # Returning if emasurements were flipped
    ########################################

    return meas_mask.astype(bool), det_mask.astype(bool), obs_mask.astype(bool), sign_keeper, gamma_keeper


"""
                elif error_type == "DEPOL2":

                    # Choose Pauli for depol
                    pauli_random = np.random.rand()

                    if 
                    chosen_pauli = 0

                    # Iterate over both qubits
                    for i in 

                    flip.set_pauli_flip(chosen_pauli, qubit_index = current_qubit, instance_index = 0)
                """