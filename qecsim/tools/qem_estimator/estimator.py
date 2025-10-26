import numpy as np
import stim

from .circuit_pauli_injector import pauli_injector
from .noise_array import return_noise_pos

__all__ = ["mc_estimator"]


def mc_estimator(
    *,
    circuit: stim.Circuit,
    circuit_noiseless: stim.Circuit,
    real_measurement: bool = False,
    number_samples: int = 50,
    number_batches: int = 1,
) -> np.array:
    #############################################################
    # Init measurement result vec and sampler of original Circuit
    #############################################################

    # Build the sampler which is needed (i.e. real measurement or measurements flipped?)
    if real_measurement is False:
        sampler_i = circuit.compile_detector_sampler()
    else:
        sampler_i = circuit.compile_sampler()

    #############################################
    # Determine qubit indexes of logical operator
    #############################################

    # init measurement index
    ms_idx: list[list] = []

    # Determine all emasurement refrences
    for _tick, inst in enumerate(circuit.flattened()):
        if inst.name in {"OBSERVABLE_INCLUDE"}:
            ms_idx.append([i.value for i in inst.targets_copy()])

    # init qubit index
    _qbt_index: list = []

    # going through the circuit until one has gotten to the i-th measurement
    for _tick, inst in enumerate(reversed(circuit.flattened())):
        if inst.name in {"M", "MR", "MX", "MZ"}:
            # Currently not used
            _curr_qubit_meas = [i.value for i in inst.targets_copy()]

            """
            if check current measurement if we have one which is included in ms_idx
            """

        """
        Else shorten qubit index with -current + len(curr_qubit_meas) and then continue 
        looking in the next measurement operator
        """

    ###################################################################
    # Calc the different possible error possibilities and their weights
    ###################################################################

    """
    Determine the postition and kind of error
    -> Also run through the list and calc every needed weight and gamma needed

    Info: .flattend() just means that stim does not simplify the circuit by adding repeat blocks
    -> This is needed here as we would need additional if conditions for reap blocks otherwise
    """

    noisy_info = return_noise_pos(circuit=circuit.flattened())

    ##############################################################################
    # Run Flip Simulator to get information on which measurements would be flipped
    ##############################################################################

    """
    -> What we should do here ist the following
        1. Run n shots
        2. In each shot go to each postion where there is the noisy gate and apply 
        the inverse channel with the given weight of the postion (Or rather of this channel)
        3. After full implementing the Flip Simulator check what operators have flipped
        4. XOR these Operations onto the noisy run of the normal noisy sample
    """

    # Create Vector for mean result per shot
    mean_res = np.zeros(number_samples)

    for i in range(number_samples):
        block_sum = 0

        ###########################################
        # Sample form the identity Circuit one shot
        ###########################################

        noisy_samples = sampler_i.sample(shots=number_batches, append_observables=True)

        ####################
        # Run pauli injector
        ####################

        pauli_frame_meas = pauli_injector(
            circuit=circuit_noiseless.flattened(),
            noise_info=noisy_info,
            batch_size=number_batches,
        )

        for current_batch in range(number_batches):
            # Extract Frames
            # meas_frame is not used; omit to reduce lint noise
            log_frame: np.array = pauli_frame_meas[2][:, current_batch]
            sgn_keeper: np.array = pauli_frame_meas[3][:, current_batch]
            gamma_frame: np.array = pauli_frame_meas[4][:, current_batch]

            #########################
            # XOR Measurement results
            #########################

            if real_measurement is False:
                # Multiplying all gammas/signs together as we look at the log operator
                """
                NEEDS REWORK
                -> BRAKES if logical Operator is not on all qubits!
                """

                full_gmsgn = 1

                for value in gamma_frame:
                    full_gmsgn *= value

                for value in sgn_keeper:
                    full_gmsgn *= value

                # XOR logical Operator and multiply all gammas to it
                obs_bits = noisy_samples[current_batch, -circuit.num_observables :]
                xored_vec = obs_bits.astype(bool) ^ log_frame
                xored_res = bool(xored_vec[0])

                # Convert into non boolian result
                non_bool_res = 0

                if xored_res:
                    non_bool_res = -1 * full_gmsgn
                else:
                    non_bool_res = 1 * full_gmsgn

            else:
                raise ValueError("Currently not supported")

            ##########################
            # Return Final Measurement
            ##########################

            block_sum += non_bool_res

        mean_res[i] = block_sum / float(number_batches)

    return mean_res
