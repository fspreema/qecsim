import numpy as np
import stim
import pymatching

__all__ = ["logical_estimator"]

class CircuitAssets:
    def __init__(self, circuit, shots):

        # Build DEM and Matcher
        dem = circuit.detector_error_model(decompose_errors=True)
        self.matcher = pymatching.Matching.from_detector_error_model(dem)

        # Compile DEM sampler and sample
        self.sampler = circuit.compile_detector_sampler()
        dets, obs = self.sampler.sample(shots, separate_observables=True)

        # Match and get prediction flip on observable
        self.pred  = self.matcher.decode_batch(dets)
        self.obs   = obs

def _only_diag(ptm: np.ndarray) -> tuple[float, dict, dict]:

    """
    Builds the pauli basis (I,X,Y,Z), solves equation (By hand implemented) to get eta
    
    Returns:
        * Gamma: Overhead required
        * Probs: Probabilities of each Clifford Basis
        * Signs: Each sign of the corresponding Paulis that go into the weight of the sampling
    """

    # Create Inverse of Diagonal ptm Matrix
    inv_x = 1/ptm[0,0]
    inv_y = 1/ptm[1,1]
    inv_z = 1/ptm[2,2]

    # Calc linear combination of elementary matricees X Y Z I are needed to build ptm 
    # -> Quasi Probabilities
    eta_I = 1/2 * (inv_x + inv_y)
    eta_X = 1/2 * (inv_x - inv_z)
    eta_Y = 1/2 * (inv_y - inv_z)
    eta_Z = 0

    eta_dict = {'I': eta_I, 'X': eta_X, 'Y': eta_Y, 'Z': eta_Z}

    # Compute overhead and probs for sampling
    gamma = np.sum(abs(v) for v in eta_dict.values())
    probs = {k: abs(v)/gamma for k,v in eta_dict.items()}
    signs = {k: 1 if v == 0 else np.sign(v) for k,v in eta_dict.items()}

    return gamma, probs, signs

# Sample Pauli
def sample_pauli(probs: dict):

    order = ['I','X','Y','Z']
    pvec  = np.array([probs[k] for k in order])

    rng = np.random.default_rng()

    return rng.choice(order, p=pvec)  

def _create_clifford_basis() -> dict:

    """
    Creates the 24 needed Cliffords inside a dictionary with all permutations included
    """



    return None

def _full_mtrx(ptm: np.ndarray) -> tuple[float, dict, dict]:

    """
    Builds the 24 needed Cliffords Basis, solves equation with minimal overhead required(L1)
    
    Returns:
        * Gamma: Overhead required
        * Probs: Probabilities of each Clifford Basis
        * Signs: Each sign of the corresponding Cliffords that go into the weight of the sampling
    """
    #########################
    # Invert Full ptm Matrix:
    #########################

    ptm_inv = np.linalg.inv(ptm)

    ##########################
    # Build Equation and solve
    ##########################
    """
    We solve here ptm_inv = sum(eta_clifford * ptm)
    -> Many solutions exist, choose such that eta is minimal (Minimal Overhead)
    """

    # Build equation

    # Solve -> eta is vec of 24 entries
    eta_values, *_ = np.linalg.lstsq()

    # Compute overhead and probs for sampling
    gamma = np.sum(abs(v) for v in eta_values)
    probs = {}
    signs = {}

    return gamma, probs, signs

def logical_estimator(*, 
                      ptm: np.ndarray, 
                      ptm_circuit_x: stim.Circuit,
                      ptm_circuit_y: stim.Circuit,
                      ptm_circuit_z: stim.Circuit,
                      logical_obs: str,
                      logical_frame: bool,
                      shots: int) -> tuple:
    
    """
    Estimates the choosen logical observable on the given circuit

    Arguemnts:
        * ptm: Pauli Transfer Matrix which got build by smapling the circuit ptm_circuit
        * ptm_circuit: Circuit which gets sampled -> F.ex memory in x/z/y basis
        * logical_obs: Logical Observable which is measured at the end and which expectation value is returned at the end
        * logical_frame: What logical value does the state have i.e. logical gate applied or not?
        * shots: Number of shots this is repeated in order to ahve a big sampling pool

    Returns:
        * Expectation-Value on a given observable which has been given by logical_obs
        * Error estimate
    """

    ################################
    # Check if ptm has off-Diagonals
    ################################

    use_cliffords = False
    eps = 1e-3

    #for i in range(len(ptm)):
    #    for j in range(len(ptm)):
            #if ptm[i,j] >= eps and i != j:
                #use_cliffords = True
    
    #----------------------------------------------------
    # Do only Diagonals (quasi probs eta only of I,X,Y,Z)
    #----------------------------------------------------

    if use_cliffords == False:

        gamma, probs, signs = _only_diag(ptm= ptm)

        ##########################
        # Precompute Frame updater
        ##########################
        
        """
        In Later runs with Lattice Surgery this needs to be updated as we measure mutli qubit Paulis
        -> i.e. X -> XX or ZZ -> Z etc.
        """

        conj_sign = {
            ('I','X'): +1, ('I','Y'): +1, ('I','Z'): +1,
            ('X','X'): +1, ('X','Y'): -1, ('X','Z'): -1,
            ('Y','X'): -1, ('Y','Y'): +1, ('Y','Z'): -1,
            ('Z','X'): -1, ('Z','Y'): -1, ('Z','Z'): +1,
        }

    #-----------------------
    # Full 24 Clifford Basis
    #-----------------------

    else:

        gamma, probs, signs = _full_mtrx(ptm = ptm)


        # Gather inforamtion about every memory basis
        # -> All are needed due to clifford mapping (i.e. H x -> z etc.)

    #####################################################
    # Adding Dictionary with all Circuits and needed info
    #####################################################

    circuit_basis = {
            'X': CircuitAssets(ptm_circuit_x, shots),
            'Y': CircuitAssets(ptm_circuit_y, shots),
            'Z': CircuitAssets(ptm_circuit_y, shots),
        }

    ###########################
    # Sample shots from circuit
    ###########################

    # In what state was it initilized?
    clean_meas = logical_frame
    all_contributions : list = []

    #--------------------------------------------
    # Sample from the same circuit as obs measured
    #--------------------------------------------

    if use_cliffords == False:

        # Current measurement circuit stays the same
        current_circuit = circuit_basis[logical_obs]

        for curr_shot in range(shots):

            # Choose what logical frame update given by upper probs
            curr_pauli = sample_pauli(probs= probs)

            # Determine current noisy Operator states (i.e. xor obs from det sample with noiseless Measurement outcome)
            noisy_meas = current_circuit.obs[curr_shot,0] ^ clean_meas

            # XOR flip with noiseless Measurement
            # -> Taking first entry for first logical observable
            final_meas = 1 - 2 * (current_circuit.pred[curr_shot,0].astype(np.int8) ^ np.int8(noisy_meas))

            # Calculate Current Weight and with that shot based result
            weight = signs[curr_pauli] * gamma * conj_sign[curr_pauli, logical_obs]

            all_contributions.append(weight * final_meas)

    #-----------------------------------------------
    # Sample from different Circuits (Non Diagonals)
    #-----------------------------------------------

    ###############################
    # Calculate Estimate and Error:
    ###############################

    logical_estimate = np.sum(all_contributions) / shots

    if shots > 1:

        inner_term = 0
        for curr_shot in range(shots):
            inner_term += (all_contributions[curr_shot] - logical_estimate) ** 2

        sampling_err = np.sqrt(1/(shots * (shots - 1)) * inner_term)

        return logical_estimate, sampling_err
    
    else:
        return logical_estimate

