import numpy as np
import pymatching
import stim

__all__ = ["calc_ptm"]


def _define_off_diag_log_rec(samples_from_sampler: np.ndarray, curr_shot: int, rec_pos: list[int]) -> np.ndarray:
    """
    Returns per shot full logical measurement by xoring the needed measurements
    -> All emasurements provided by samples_from_sampler
    -> Index of measurements needed to build the logical measurement outcome
    """

    final_meas = -99

    for curr_rec in rec_pos:
        # Check if anything is added right now
        if final_meas == -99:
            final_meas = samples_from_sampler[curr_shot, curr_rec]

        else:
            # XOR
            final_meas ^= samples_from_sampler[curr_shot, curr_rec]

    return final_meas


def _xor_meas_and_decoder(decoder_prediction: np.ndarray, logical_state_meas: list, curr_shot: int) -> int:
    """
    Returns the final measurement for the current shot
    -> XOR logical measurement with the decoder prediction iof the observable was flipped
    -> i.e. correct measured observable through decoder prediction
    """

    final_log_meas = 1 - 2 * (decoder_prediction[curr_shot, 0].astype(np.int8) ^ np.int8(logical_state_meas[curr_shot]))
    return final_log_meas


def _build_dem(circuit: stim.Circuit, dem_needed: bool = False):
    """
    Creates for a given stim Circuit the detector error model as well as the pymachting matcher
    """

    dem = circuit.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(dem)

    if dem_needed:
        return matcher, dem

    else:
        return matcher


def calc_ptm(
    *,
    circuit_0_z: stim.Circuit,
    circuit_p_z: stim.Circuit,
    circuit_m_z: stim.Circuit,
    circuit_pi_z: stim.Circuit,
    circuit_mi_z: stim.Circuit,
    circuit_0_x: stim.Circuit,
    circuit_1_x: stim.Circuit,
    circuit_p_x: stim.Circuit,
    circuit_pi_x: stim.Circuit,
    circuit_mi_x: stim.Circuit,
    circuit_pi_y: stim.Circuit,
    circuit_0_y: stim.Circuit,
    circuit_1_y: stim.Circuit,
    circuit_p_y: stim.Circuit,
    circuit_m_y: stim.Circuit,
    rec_pos_log_x: list[int],
    rec_pos_log_z: list[int],
    rec_pos_log_y01: list[int],
    rec_pos_log_ypm: list[int],
    rec_pos_log_z_pmi: list[int],
    rec_pos_log_x_pmi: list[int],
    samples: int = 1_000,
) -> np.ndarray:
    """
    Building PTM by building each possible memory circuit and measurement combination
    -> + logical -> Measuring in x/z/y basis
    -> - logical -> Measuring in x/z/y basis
    -> 0 logical -> Measuring in x/z/y basis
    -> 1 logical -> Measuring in x/z/y basis
    -> +i logical -> Measuring in x/z/y basis

    Since 0 and 1 logical or + and - etc. only differ by logical_operator
    -> Only one circuit is required (Logical is inside Pauli Frame)

    What we do:
        * Determine the noise clean expectation of the measurement
        * Compute the noisy current observable status (xor noiseless initlial state by obs state)
        * Run QEC and check whether the noisy observable should be flipped or not
        * Flip the observable i.e. the clean expectation value
        * Build up PTM out of the expectaition values

    Arguemnts:
        *assume_pauli: Only the 3 diagonal Values are calculated

    """

    #####################################
    # Builing Diagonal Entries of the PTM
    #####################################

    # Correct nosieless Measurements
    """
    Does this even make sense to define so specifically?
        -> We have only a pauli correction which could be done classically
        -> Initilizing in +x or -x or +z or -z makes no difference (Can be calssically correct, no?)
    """
    clean_meas_0_z = False
    clean_meas_1_z = True
    clean_meas_p_x = False
    clean_meas_m_x = True
    clean_meas_pi_y = False
    clean_meas_mi_y = True

    # Build a graphlike DEM and a matcher
    m_0_z = _build_dem(circuit_0_z)
    m_p_x = _build_dem(circuit_p_x)
    m_pi_y = _build_dem(circuit_pi_y)

    # Take m shots of detection events
    sampler_0_z = circuit_0_z.compile_detector_sampler()
    dets_0_z, obs_0_z = sampler_0_z.sample(samples, separate_observables=True)

    sampler_p_x = circuit_p_x.compile_detector_sampler()
    dets_p_x, obs_p_x = sampler_p_x.sample(samples, separate_observables=True)

    sampler_pi_y = circuit_pi_y.compile_detector_sampler()
    dets_pi_y, obs_pi_y = sampler_pi_y.sample(samples, separate_observables=True)

    # Determine current noisy Operator states (i.e. xor obs from det sample with noiseless Measurement outcome)
    """
    THIS IS TOO MUCH WORK AND CAN BE SIMPLIFIED
    -> Just return freom each circuit the log rec tragets not only for the non determinstic emasurements
    -> Do analogue method than for the off diagonals
    -> This works but is more complicated and results in two methods beeing used in one function which can be mitigated
    """

    noisy_meas_0_z = []
    noisy_meas_1_z = []
    noisy_meas_p_x = []
    noisy_meas_m_x = []
    noisy_meas_pi_y = []
    noisy_meas_mi_y = []

    for curr_shot in range(samples):
        """
        Here we now use the same observable flip for 0 or 1 etc. and only flip with pauli frame
        """
        noisy_meas_0_z.append(obs_0_z[curr_shot, 0] ^ clean_meas_0_z)
        noisy_meas_1_z.append(obs_0_z[curr_shot, 0] ^ clean_meas_1_z)
        noisy_meas_p_x.append(obs_p_x[curr_shot, 0] ^ clean_meas_p_x)
        noisy_meas_m_x.append(obs_p_x[curr_shot, 0] ^ clean_meas_m_x)
        noisy_meas_pi_y.append(obs_pi_y[curr_shot, 0] ^ clean_meas_pi_y)
        noisy_meas_mi_y.append(obs_pi_y[curr_shot, 0] ^ clean_meas_mi_y)

    # Decode all sample round in one Batch decode
    pred_0_z = m_0_z.decode_batch(dets_0_z)
    pred_p_x = m_p_x.decode_batch(dets_p_x)
    pred_pi_y = m_pi_y.decode_batch(dets_pi_y)

    # XOR flip with noiseless Measurement
    # -> Taking first entry for first logical observable
    final_meas_0_z = []
    final_meas_1_z = []
    final_meas_p_x = []
    final_meas_m_x = []
    final_meas_pi_y = []
    final_meas_mi_y = []

    for curr_shot in range(samples):
        """
        Here we now use the same decoder predictions for 0 or 1 etc. and only flip with pauli frame
        """
        final_meas_0_z.append(_xor_meas_and_decoder(pred_0_z, noisy_meas_0_z, curr_shot))
        final_meas_1_z.append(_xor_meas_and_decoder(pred_0_z, noisy_meas_1_z, curr_shot))
        final_meas_p_x.append(_xor_meas_and_decoder(pred_p_x, noisy_meas_p_x, curr_shot))
        final_meas_m_x.append(_xor_meas_and_decoder(pred_p_x, noisy_meas_m_x, curr_shot))
        final_meas_pi_y.append(_xor_meas_and_decoder(pred_pi_y, noisy_meas_pi_y, curr_shot))
        final_meas_mi_y.append(_xor_meas_and_decoder(pred_pi_y, noisy_meas_mi_y, curr_shot))

    # Take the Average of each of them
    mu_z_pz = np.average(final_meas_0_z)
    mu_z_mz = np.average(final_meas_1_z)
    mu_x_px = np.average(final_meas_p_x)
    mu_x_mx = np.average(final_meas_m_x)
    mu_y_py = np.average(final_meas_pi_y)
    mu_y_my = np.average(final_meas_mi_y)

    # Calc PTM entries
    r_zz = 1 / 2 * (mu_z_pz - mu_z_mz)
    r_xx = 1 / 2 * (mu_x_px - mu_x_mx)
    r_yy = 1 / 2 * (mu_y_py - mu_y_my)

    ######################################
    # Building Off-Diagonals of the Matrix
    ######################################

    """
    For the Off-Diagonals we can't directly use the DEM as we need the 
    raw measurements to infer what logical state we have

    -> We use compile sampler to infer the logical state
    -> Convert into a dem to run the matching
    """

    # Build the normal measurement smaples and sample n shots
    smpl_0_x = circuit_0_x.compile_sampler()
    rstls_smpls_0_x = smpl_0_x.sample(shots=samples)
    smpl_0_y = circuit_0_y.compile_sampler()
    rstls_smpls_0_y = smpl_0_y.sample(shots=samples)
    smpl_1_x = circuit_1_x.compile_sampler()
    rstls_smpls_1_x = smpl_1_x.sample(shots=samples)
    smpl_1_y = circuit_1_y.compile_sampler()
    rstls_smpls_1_y = smpl_1_y.sample(shots=samples)
    smpl_p_z = circuit_p_z.compile_sampler()
    rstls_smpls_p_z = smpl_p_z.sample(shots=samples)
    smpl_p_y = circuit_p_y.compile_sampler()
    rstls_smpls_p_y = smpl_p_y.sample(shots=samples)
    smpl_m_y = circuit_m_y.compile_sampler()
    rstls_smpls_m_y = smpl_m_y.sample(shots=samples)
    smpl_m_z = circuit_m_z.compile_sampler()
    rstls_smpls_m_z = smpl_m_z.sample(shots=samples)
    smpl_pi_z = circuit_pi_z.compile_sampler()
    rstls_smpls_pi_z = smpl_pi_z.sample(shots=samples)
    smpl_pi_x = circuit_pi_x.compile_sampler()
    rstls_smpls_pi_x = smpl_pi_x.sample(shots=samples)
    smpl_mi_z = circuit_mi_z.compile_sampler()
    rstls_smpls_mi_z = smpl_mi_z.sample(shots=samples)
    smpl_mi_x = circuit_mi_x.compile_sampler()
    rstls_smpls_mi_x = smpl_mi_x.sample(shots=samples)

    # Getting the individual logical states after measurement
    log_state_0_x = []
    log_state_0_y = []
    log_state_1_x = []
    log_state_1_y = []
    log_state_p_y = []
    log_state_m_y = []
    log_state_p_z = []
    log_state_m_z = []
    log_state_pi_x = []
    log_state_mi_x = []
    log_state_pi_z = []
    log_state_mi_z = []

    for curr_shot in range(samples):
        measurement_res_0_x = _define_off_diag_log_rec(rstls_smpls_0_x, curr_shot, rec_pos_log_x)
        measurement_res_1_x = _define_off_diag_log_rec(rstls_smpls_1_x, curr_shot, rec_pos_log_x)
        measurement_res_p_z = _define_off_diag_log_rec(rstls_smpls_p_z, curr_shot, rec_pos_log_z)
        measurement_res_m_z = _define_off_diag_log_rec(rstls_smpls_m_z, curr_shot, rec_pos_log_z)
        measurement_res_0_y = _define_off_diag_log_rec(rstls_smpls_0_y, curr_shot, rec_pos_log_y01)
        measurement_res_1_y = _define_off_diag_log_rec(rstls_smpls_1_y, curr_shot, rec_pos_log_y01)
        measurement_res_p_y = _define_off_diag_log_rec(rstls_smpls_p_y, curr_shot, rec_pos_log_ypm)
        measurement_res_m_y = _define_off_diag_log_rec(rstls_smpls_m_y, curr_shot, rec_pos_log_ypm)
        measurement_res_pi_x = _define_off_diag_log_rec(rstls_smpls_pi_x, curr_shot, rec_pos_log_x_pmi)
        measurement_res_mi_x = _define_off_diag_log_rec(rstls_smpls_mi_x, curr_shot, rec_pos_log_x_pmi)
        measurement_res_pi_z = _define_off_diag_log_rec(rstls_smpls_pi_z, curr_shot, rec_pos_log_z_pmi)
        measurement_res_mi_z = _define_off_diag_log_rec(rstls_smpls_mi_z, curr_shot, rec_pos_log_z_pmi)

        # Appending the current measurements to the list
        log_state_0_x.append(measurement_res_0_x)
        log_state_1_x.append(measurement_res_1_x)
        log_state_p_z.append(measurement_res_p_z)
        log_state_m_z.append(measurement_res_m_z)
        log_state_0_y.append(measurement_res_0_y)
        log_state_1_y.append(measurement_res_1_y)
        log_state_p_y.append(measurement_res_p_y)
        log_state_m_y.append(measurement_res_m_y)
        log_state_pi_x.append(measurement_res_pi_x)
        log_state_mi_x.append(measurement_res_mi_x)
        log_state_pi_z.append(measurement_res_pi_z)
        log_state_mi_z.append(measurement_res_mi_z)

    # Build a graphlike DEM and a matcher
    m_0_x, dem_0_x = _build_dem(circuit_0_x, dem_needed=True)
    m_1_x, dem_1_x = _build_dem(circuit_1_x, dem_needed=True)
    m_p_z, dem_p_z = _build_dem(circuit_p_z, dem_needed=True)
    m_m_z, dem_m_z = _build_dem(circuit_m_z, dem_needed=True)
    m_0_y, dem_0_y = _build_dem(circuit_0_y, dem_needed=True)
    m_1_y, dem_1_y = _build_dem(circuit_1_y, dem_needed=True)
    m_p_y, dem_p_y = _build_dem(circuit_p_y, dem_needed=True)
    m_m_y, dem_m_y = _build_dem(circuit_m_y, dem_needed=True)
    m_pi_x, dem_pi_x = _build_dem(circuit_pi_x, dem_needed=True)
    m_mi_x, dem_mi_x = _build_dem(circuit_mi_x, dem_needed=True)
    m_pi_z, dem_pi_z = _build_dem(circuit_pi_z, dem_needed=True)
    m_mi_z, dem_mi_z = _build_dem(circuit_mi_z, dem_needed=True)

    # Converting the measurement sample into DEM sample and continue as usual with decoding
    cvrtr_0_x = circuit_0_x.compile_m2d_converter()
    dets_ops_0_x = cvrtr_0_x.convert(measurements=rstls_smpls_0_x, append_observables=True)
    num_dets = dem_0_x.num_detectors
    dets_0_x = dets_ops_0_x[:, :num_dets]

    cvrtr_1_x = circuit_1_x.compile_m2d_converter()
    dets_ops_1_x = cvrtr_1_x.convert(measurements=rstls_smpls_1_x, append_observables=True)
    num_dets = dem_1_x.num_detectors
    dets_1_x = dets_ops_1_x[:, :num_dets]

    cvrtr_p_z = circuit_p_z.compile_m2d_converter()
    dets_ops_p_z = cvrtr_p_z.convert(measurements=rstls_smpls_p_z, append_observables=True)
    num_dets = dem_p_z.num_detectors
    dets_p_z = dets_ops_p_z[:, :num_dets]

    cvrtr_m_z = circuit_m_z.compile_m2d_converter()
    dets_ops_m_z = cvrtr_m_z.convert(measurements=rstls_smpls_m_z, append_observables=True)
    num_dets = dem_m_z.num_detectors
    dets_m_z = dets_ops_m_z[:, :num_dets]

    cvrtr_0_y = circuit_0_y.compile_m2d_converter()
    dets_ops_0_y = cvrtr_0_y.convert(measurements=rstls_smpls_0_y, append_observables=True)
    num_dets = dem_0_y.num_detectors
    dets_0_y = dets_ops_0_y[:, :num_dets]

    cvrtr_1_y = circuit_1_y.compile_m2d_converter()
    dets_ops_1_y = cvrtr_1_y.convert(measurements=rstls_smpls_1_y, append_observables=True)
    num_dets = dem_1_y.num_detectors
    dets_1_y = dets_ops_1_y[:, :num_dets]

    cvrtr_p_y = circuit_p_y.compile_m2d_converter()
    dets_ops_p_y = cvrtr_p_y.convert(measurements=rstls_smpls_p_y, append_observables=True)
    num_dets = dem_p_y.num_detectors
    dets_p_y = dets_ops_p_y[:, :num_dets]

    cvrtr_m_y = circuit_m_y.compile_m2d_converter()
    dets_ops_m_y = cvrtr_m_y.convert(measurements=rstls_smpls_m_y, append_observables=True)
    num_dets = dem_m_y.num_detectors
    dets_m_y = dets_ops_m_y[:, :num_dets]

    cvrtr_pi_x = circuit_pi_x.compile_m2d_converter()
    dets_ops_pi_x = cvrtr_pi_x.convert(measurements=rstls_smpls_pi_x, append_observables=True)
    num_dets = dem_pi_x.num_detectors
    dets_pi_x = dets_ops_pi_x[:, :num_dets]

    cvrtr_mi_x = circuit_mi_x.compile_m2d_converter()
    dets_ops_mi_x = cvrtr_mi_x.convert(measurements=rstls_smpls_mi_x, append_observables=True)
    num_dets = dem_mi_x.num_detectors
    dets_mi_x = dets_ops_mi_x[:, :num_dets]

    cvrtr_pi_z = circuit_pi_z.compile_m2d_converter()
    dets_ops_pi_z = cvrtr_pi_z.convert(measurements=rstls_smpls_pi_z, append_observables=True)
    num_dets = dem_pi_z.num_detectors
    dets_pi_z = dets_ops_pi_z[:, :num_dets]

    cvrtr_mi_z = circuit_mi_z.compile_m2d_converter()
    dets_ops_mi_z = cvrtr_mi_z.convert(measurements=rstls_smpls_mi_z, append_observables=True)
    num_dets = dem_mi_z.num_detectors
    dets_mi_z = dets_ops_mi_z[:, :num_dets]

    # Decode all sample round in one Batch decode
    pred_0_x = m_0_x.decode_batch(dets_0_x)
    pred_1_x = m_1_x.decode_batch(dets_1_x)

    pred_p_z = m_p_z.decode_batch(dets_p_z)
    pred_m_z = m_m_z.decode_batch(dets_m_z)

    pred_0_y = m_0_y.decode_batch(dets_0_y)
    pred_1_y = m_1_y.decode_batch(dets_1_y)

    pred_p_y = m_p_y.decode_batch(dets_p_y)
    pred_m_y = m_m_y.decode_batch(dets_m_y)

    pred_pi_x = m_pi_x.decode_batch(dets_pi_x)
    pred_mi_x = m_mi_x.decode_batch(dets_mi_x)

    pred_pi_z = m_pi_z.decode_batch(dets_pi_z)
    pred_mi_z = m_mi_z.decode_batch(dets_mi_z)

    # Init final state List after EC
    final_meas_0_x: list = []
    final_meas_1_x: list = []

    final_meas_p_z: list = []
    final_meas_m_z: list = []

    final_meas_0_y: list = []
    final_meas_1_y: list = []

    final_meas_p_y: list = []
    final_meas_m_y: list = []

    final_meas_pi_x: list = []
    final_meas_mi_x: list = []

    final_meas_pi_z: list = []
    final_meas_mi_z: list = []

    # XOR flip with noiseless Measurement
    # (looping over all samples and compoaring to the current logical outcome in the list)
    # -> Taking first entry for first logical observable
    for curr_shot in range(samples):
        # Adding final XORed measurement into the list
        final_meas_0_x.append(_xor_meas_and_decoder(pred_0_x, log_state_0_x, curr_shot))
        final_meas_1_x.append(_xor_meas_and_decoder(pred_1_x, log_state_1_x, curr_shot))

        final_meas_p_z.append(_xor_meas_and_decoder(pred_p_z, log_state_p_z, curr_shot))
        final_meas_m_z.append(_xor_meas_and_decoder(pred_m_z, log_state_m_z, curr_shot))

        final_meas_0_y.append(_xor_meas_and_decoder(pred_0_y, log_state_0_y, curr_shot))
        final_meas_1_y.append(_xor_meas_and_decoder(pred_1_y, log_state_1_y, curr_shot))

        final_meas_p_y.append(_xor_meas_and_decoder(pred_p_y, log_state_p_y, curr_shot))
        final_meas_m_y.append(_xor_meas_and_decoder(pred_m_y, log_state_m_y, curr_shot))

        final_meas_pi_x.append(_xor_meas_and_decoder(pred_pi_x, log_state_pi_x, curr_shot))
        final_meas_mi_x.append(_xor_meas_and_decoder(pred_mi_x, log_state_mi_x, curr_shot))

        final_meas_pi_z.append(_xor_meas_and_decoder(pred_pi_z, log_state_pi_z, curr_shot))
        final_meas_mi_z.append(_xor_meas_and_decoder(pred_mi_z, log_state_mi_z, curr_shot))

    # Take the Average of each of them
    mu_x_pz = np.average(final_meas_0_x)
    mu_x_mz = np.average(final_meas_1_x)
    mu_z_px = np.average(final_meas_p_z)
    mu_z_mx = np.average(final_meas_m_z)
    mu_y_pz = np.average(final_meas_0_y)
    mu_y_mz = np.average(final_meas_1_y)
    mu_y_px = np.average(final_meas_p_y)
    mu_y_mx = np.average(final_meas_m_y)
    mu_x_py = np.average(final_meas_pi_x)
    mu_x_my = np.average(final_meas_mi_x)
    mu_z_py = np.average(final_meas_pi_z)
    mu_z_my = np.average(final_meas_mi_z)

    # Calc PTM entries
    r_zx = 1 / 2 * (mu_z_px - mu_z_mx)
    r_xz = 1 / 2 * (mu_x_pz - mu_x_mz)
    r_xy = 1 / 2 * (mu_y_px - mu_y_mx)
    r_zy = 1 / 2 * (mu_y_pz - mu_y_mz)
    r_yx = 1 / 2 * (mu_x_py - mu_x_my)
    r_yz = 1 / 2 * (mu_z_py - mu_z_my)

    # Create PTM Matrix
    ptm = [[r_xx, r_xy, r_xz], [r_yx, r_yy, r_yz], [r_zx, r_zy, r_zz]]

    return np.array(ptm)
