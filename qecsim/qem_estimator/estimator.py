import stim
import numpy as np
import math

from .circuit_pauli_injector import pauli_injector
from .noise_array import return_noise_pos

__all__ = ["mc_estimator"]

def mc_estimator(circuit : stim.Circuit, circuit_noiseless : stim.Circuit, number_samples : int = 50, 
                  number_shots : int = 50 ) -> np.array:

    #############################################################
    # Init measurement result vec and sampler of original Circuit
    #############################################################

    meas_result = np.zeros(number_samples)
    sampler_I = circuit.compile_detector_sampler()

    ###################################################################
    # Calc the different possible error possibilities and their weights
    ###################################################################

    """
    Determine the postition and kind of error
    -> Also run through the list and calc every needed weight and gamma needed
    """

    noisy_info = return_noise_pos(circuit = circuit)

    ##############################################################################
    # Run Flip Simulator to get information on which measurements would be flipped
    ##############################################################################

    """
    -> What we should do here ist the following
        1. Run n shots
        2. In each shot go to each postion where there is the noisy gate and apply the inverse channel with the given weight of the postion (Or rather of this channel)
        3. After full implementing the Flip Simulator check what operators have flipped
        4. XOR these Operations onto the noisy run of the normal noisy sample
    """

    ########################################
    # Create Vector for mean result per shot
    ########################################

    mean_res = np.zeros(number_samples)

    for i in range(number_samples):

        ###########################################
        # Sample form the identity Circuit one shot
        ###########################################

        noisy_samples = sampler_I.sample(shots = 1, append_observables = True)

        ####################
        # Run pauli injector
        ####################

        """
        Maybe we runs x instance of the injector and we take the mean of the result?
        """

        pauli_frame_meas = pauli_injector(circuit = circuit_noiseless, noise_info = noisy_info)

        # Extract Frames
        meas_frame : np.array = pauli_frame_meas[0]
        log_frame : np.array = pauli_frame_meas[2]
        sgn_keeper : np.array = pauli_frame_meas[3]

        #########################
        # XOR Measurement results
        #########################

        # Extract gamma
        gamma_frame : np.array = pauli_frame_meas[4]
        
        # Multiplying all gammas/signs together as we look at the log operator
        full_gmsgn = 1

        for value in gamma_frame:
            full_gmsgn *= value

        for value in sgn_keeper:
            full_gmsgn *= value

        # XOR logical Operator and multiply all gammas to it
        obs_bits = noisy_samples[0, -circuit.num_observables:]
        xored_vec = obs_bits.astype(bool) ^ log_frame     
        xored_res = bool(xored_vec[0])       

        # Convert into non boolian result
        non_bool_res = 0

        if xored_res == True:
            non_bool_res = -1 * full_gmsgn
        else:
            non_bool_res = 1 * full_gmsgn

        ##########################
        # Return Final Measurement
        ##########################

        mean_res[i] = non_bool_res

    return mean_res