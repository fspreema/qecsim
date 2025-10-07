import numpy as np
import cvxpy as cp
import itertools
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
def _sample_pauli_or_cliff(probs: dict) -> str:
    """
    Sample any Pauli or Clifford with given probability

    Returns:
        * Str: Selected Pauli for current sampling
    """

    possible_gates = list(probs.keys())
    pvec  = np.array([probs[k] for k in possible_gates])

    rng = np.random.default_rng()

    return rng.choice(possible_gates, p=pvec)

def _even_permutation_check(sigma: tuple[int, int, int]) -> int:
    """
    Returns a +1 if given permutiation is even coming from (0,1,2)
    """

    inversions = 0
    for i in range(3):
        for j in range(i+1, 3):
            if sigma[i] > sigma[j]:
                inversions +=1

    if inversions % 2 == 0:
        return +1
    else:
        return -1
    
def _create_matrix(sigma: tuple[int, int, int],
                  signs: tuple[int, int, int]) -> np.ndarray:
    
    """
    Returns the 3x3 Matrix used for the solver
    -> Describes how the Operator acts on the Paulis i.e. how do each pauli get transformed?

    Example:
    -> sigma = (2,0,1) i.e. x -> z y -> x z -> x
    -> signs = (-1, +1, +1)

    -> We receive: np.ndarray([0,1,0],
                              [0,0,-1],
                              [-1,0,0])
    """

    mat = np.zeros((3,3), dtype= int)
    for j in range(3):
        i = sigma[j]
        curr_sign = signs[j]
        mat[i,j] = curr_sign

    return mat

def _perm_sign_dict_conv(sigma: tuple[int, int, int],
                        signs: tuple[int, int, int]) -> tuple[dict, dict]:
    
    """
    Converts the sigma and sign convention into proper dictionaries

    Returns:
        * perm["X"] = Z
        * sign["Z"] = -1
    """

    axes = ["X", "Y", "Z"]

    # Intilize dict
    perm_dict = {}
    sgn_dict = {}

    # Calc perm dict
    for curr_idx in range(len(sigma)):
        transf_idx = sigma[curr_idx]
        perm_dict[axes[curr_idx]] = axes[transf_idx]

    # Calc sgn dict
    for curr_idx in range(len(sigma)):
        curr_sign = signs[curr_idx]
        sgn_dict[axes[curr_idx]] = curr_sign

    return perm_dict, sgn_dict

def _create_clifford_basis() -> list[dict]:

    """
    Creates the 24 needed Cliffords inside a dictionary with all permutations included
    """

    # Initilize Basis and current Cliff transf
    basis: list[dict]= []
    curr_idx = 0

    # Create all possible permutations
    all_perm = itertools.permutations(range(3))

    for curr_sigma in all_perm:
        # Calc parity
        par = _even_permutation_check(curr_sigma)

        # Try every possible sign combination
        for curr_sgns in itertools.product([-1, 1], repeat= 3):

            #Check validitiy of Transform (Only keep det = +1 -> Unitary Operations)
            sgn_x, sgn_y, sgn_z = curr_sgns

            if par * sgn_x * sgn_y * sgn_z == 1:

                # Calc perm and sgn dict:
                perm_dict, sgn_dict = _perm_sign_dict_conv(curr_sigma, curr_sgns)

                # Build tranf Matrix
                mat = _create_matrix(curr_sigma, curr_sgns)
                
                basis.append({
                    "name": f"C{curr_idx}",
                    "mat": mat,
                    "perm": perm_dict,
                    "sgn": sgn_dict
                })

                curr_idx += 1

    return basis

def _has_off_diagonals(ptm : np.ndarray) -> bool:

    off_diags = False
    eps = 1e-10

    for i in range(len(ptm)):
        for j in range(len(ptm)):
            if abs(ptm[i,j]) >= eps and i != j:
                off_diags = True

    return off_diags

def _stack_vectors_of_mtrx(mtx: np.ndarray) -> np.ndarray:
    """
    Equation which gets solved stacks the mtrx on top so 3x3 goes to 9 element vector
    """

    reshaped_mtx = mtx.reshape(9, order="F")

    return reshaped_mtx

def _full_mtrx_slv(ptm: np.ndarray) -> tuple[float, dict, dict, list]:

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
    basis = _create_clifford_basis()
    A_cols = [_stack_vectors_of_mtrx(c['mat']) for c in basis]
    A = np.column_stack(A_cols)
    b = _stack_vectors_of_mtrx(ptm_inv)

    eta = cp.Variable(24)
    eps_val: float = 0.0

    # Solve -> eta is vec of 24 entries
    constraints = [cp.norm2(A @ eta - b) <= eps_val]
    problem = cp.Problem(cp.Minimize(cp.norm1(eta)), constraints)
    problem.solve(solver=cp.ECOS, abstol=1e-9, reltol=1e-9, feastol=1e-9)
    eta_values = np.array(eta.value, dtype=float).reshape(-1)

    # Compute overhead and probs for sampling
    abs_eta = np.abs(eta_values)
    gamma = np.sum(abs(v) for v in eta_values)
    probs = {basis[k]['name']: float(abs_eta[k] / gamma) for k in range(24)}
    signs = {basis[k]['name']: (1 if eta_values[k] >= 0 else -1) for k in range(24)}

    return gamma, probs, signs, basis

def logical_estimator(*, 
                      ptm: np.ndarray, 
                      ptm_circuit_x: stim.Circuit,
                      ptm_circuit_y: stim.Circuit,
                      ptm_circuit_z: stim.Circuit,
                      logical_obs: str,
                      frame_flip: bool,
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

    use_cliffords = _has_off_diagonals(ptm)
    
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

        gamma, probs, signs, basis = _full_mtrx_slv(ptm = ptm)

    #####################################################
    # Adding Dictionary with all Circuits and needed info
    #####################################################

    circuit_basis = {
            'X': CircuitAssets(ptm_circuit_x, shots),
            'Y': CircuitAssets(ptm_circuit_y, shots),
            'Z': CircuitAssets(ptm_circuit_z, shots),
        }

    ###########################
    # Sample shots from circuit
    ###########################

    # In what state was it initilized?
    clean_meas = frame_flip
    all_contributions : list = []

    #--------------------------------------------
    # Sample from the same circuit as obs measured
    #--------------------------------------------

    if use_cliffords == False:

        # Current measurement circuit stays the same
        current_circuit = circuit_basis[logical_obs]

        for curr_shot in range(shots):

            # Choose what logical frame update given by upper probs
            curr_pauli = _sample_pauli_or_cliff(probs= probs)

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

    else: 

        name_to_elem = {c["name"]: c for c in basis}

        for curr_shot in range(shots):

            # Choose what logical frame update given by upper probs
            curr_cliff = _sample_pauli_or_cliff(probs= probs)
            curr_cliff_dict = name_to_elem[curr_cliff]

            # Get infromation what circuit needs to be measured and what sign flip we have
            which_meas_basis = curr_cliff_dict['perm'][logical_obs]
            curr_sign = curr_cliff_dict['sgn'][logical_obs]

            # Current measurement dependent on the current pauli sampled
            current_circuit = circuit_basis[which_meas_basis]

            # Determine current noisy Operator states (i.e. xor obs from det sample with noiseless Measurement outcome)
            noisy_meas = current_circuit.obs[curr_shot,0] ^ clean_meas

            # XOR flip with noiseless Measurement
            # -> Taking first entry for first logical observable
            final_meas = 1 - 2 * (current_circuit.pred[curr_shot,0].astype(np.int8) ^ np.int8(noisy_meas))

            # Calculate Current Weight and with that shot based result
            weight = signs[curr_cliff] * gamma * curr_sign

            all_contributions.append(weight * final_meas)

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

