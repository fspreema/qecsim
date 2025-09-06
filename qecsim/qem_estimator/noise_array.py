import stim
import numpy as np
import math

__all__ = ["return_noise_pos"]

def return_noise_pos(circuit : stim.Circuit) -> np.array:

    """
    Returns the kind of error, qubit-index, time slice postion, prob_i, prob_err
    -> In Monte Carlo we go to each of those locations and choose the channel depending on 
        the weight (or rather prob) of each channel
    
    IMPORTANT
    -> P_err is the prob for any error happening! -> Choose random between XYZ or IX, IY etc. in the different channels!
    """

    noise_overview : list = []

    q_already_calc_depol = False
    q_already_calc_x = False
    q_already_calc_depol2 = False

    """
    We can't use the enumeration as a TICK method as the noisy and noiseless circuit have different lengths
    -> gate_num: How many gate where before this error?
    """

    gate_num = 0

    for tick, inst in enumerate(circuit):

        # Filter Different Error Channels
        if inst.name in {"DEPOLARIZE1"}:


            if q_already_calc_depol is False:
                # Get Prob and calc q
                prob_depol = inst.gate_args_copy()[0]
                weight_depol = -1 * prob_depol / (1 - (4/3 * prob_depol))

                # Calc overhead and physical prob
                gamma_depol = abs(weight_depol) + abs(1 - weight_depol)
                prob_I_depol = abs(1  -weight_depol)/gamma_depol
                prob_err_depol = abs(weight_depol)/(gamma_depol)
                sgn_I_depol = math.copysign(1, 1 - weight_depol)
                sgn_err_depol = math.copysign(1, weight_depol)

                # Get qubit index
                q_number = tuple(t.value for t in inst.targets_copy() if t.is_qubit_target)

                # Append into List
                noise_overview.append(["DEPOL", gate_num, q_number, prob_I_depol, prob_err_depol, sgn_I_depol, sgn_err_depol, gamma_depol])

                # Setting Bool to True
                q_already_calc_depol = True

            else:
                # Get qubit index
                q_number = tuple(t.value for t in inst.targets_copy() if t.is_qubit_target)

                # Append into List
                noise_overview.append(["DEPOL", gate_num, q_number, prob_I_depol, prob_err_depol, sgn_I_depol, sgn_err_depol, gamma_depol])


        elif inst.name in {"X_ERROR"}:
            
            if q_already_calc_x is False:
                # Get Prob and calc q
                prob_x = inst.gate_args_copy()[0]
                weight_x = -1 * prob_x / (1 - (2 * prob_x))

                # Calc overhead and physical prob
                gamma_x = abs(weight_x) + abs(1 - weight_x)
                prob_I_x = abs(1 - weight_x)/gamma_x
                prob_err_x = abs(weight_x)/(gamma_x)
                sgn_I_x = math.copysign(1, 1 - weight_x)
                sgn_err_x = math.copysign(1, weight_x)

                # Get qubit index
                q_number = tuple(t.value for t in inst.targets_copy() if t.is_qubit_target)

                # Append into List
                noise_overview.append(["X_ERR", gate_num, q_number, prob_I_x, prob_err_x, sgn_I_x, sgn_err_x, gamma_x])

                # Setting Bool to True
                q_already_calc_x = True

            else:
                # Get qubit index
                q_number = tuple(t.value for t in inst.targets_copy() if t.is_qubit_target)

                # Append into List
                noise_overview.append(["X_ERR", gate_num, q_number, prob_I_x, prob_err_x, sgn_I_x, sgn_err_x, gamma_x])

        elif inst.name in {"DEPOLARIZE2"}:

            if q_already_calc_depol2 is False:
                # Get Prob and calc q
                prob_depol2 = inst.gate_args_copy()[0]
                weight_depol2 = -1 * prob_depol2 / (1 - (16/15 * prob_depol2))

                # Calc overhead and physical prob
                gamma_depol2 = abs(weight_depol2) + abs(1 - weight_depol2)
                prob_I_depol2 = abs(1 - weight_depol2)/gamma_depol2
                prob_err_depol2 = abs(weight_depol2)/(gamma_depol2)
                sgn_I_depol2 = math.copysign(1, 1 - weight_depol2)
                sgn_err_depol2 = math.copysign(1, weight_depol2)

                # Get qubit index
                q_number = tuple(t.value for t in inst.targets_copy() if t.is_qubit_target)

                # Append into List
                noise_overview.append(["DEPOL2", gate_num, q_number, prob_I_depol2, prob_err_depol2, sgn_I_depol2, sgn_err_depol2, gamma_depol2])

                # Setting Bool to True
                q_already_calc_depol2 = True

            else:
                # Get qubit index
                q_number = tuple(t.value for t in inst.targets_copy() if t.is_qubit_target)

                # Append into List
                noise_overview.append(["DEPOL2", gate_num, q_number, prob_I_depol2, prob_err_depol2, sgn_I_depol2, sgn_err_depol2, gamma_depol2])

        else:
            gate_num += 1

    ret_noise = np.array(noise_overview, dtype=object)

    return ret_noise