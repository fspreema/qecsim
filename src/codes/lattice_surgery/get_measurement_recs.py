import stim

__all__ = ["get_measurement_recs"]


def get_measurement_recs(circuit: stim.Circuit, observable_index: int) -> list[int]:
    """
    Returns the list of measurement record positions that need to be xored together
    to get the final logical measurement.
    """

    # Building final Measure Class
    if circuit is None:
        raise ValueError("Circuit has not been built yet. Please build the circuit first.")

    # Init Measurement Records list
    measurement_records = []
    check_for_measurements = False
    measurement_in_last_tick = False
    curr_meas_offset = 0

    # Loop through all instructions in the circuit and check for OBSERVABLE_INCLUDE instructions
    # We need to check if after the first observable arguments where captured no new
    # measurements are done -> If that is the case, add the measurements to calc new offsets
    # for next observable_include command

    for instruction in circuit:
        # Check Instruction name
        if instruction.name == "OBSERVABLE_INCLUDE":
            # If previously measurements detected
            if measurement_in_last_tick:
                # Shift current measurement records by applying offset and add to final
                # measurement records
                # Build new shifted list and clear current measurement records
                shifted_measurement_records = [
                    rec - curr_meas_offset for rec in measurement_records
                ]
                # Update measurement records with shifted values
                measurement_records = []
                measurement_records = shifted_measurement_records
                # Reset current measurement records and lookout for measurements
                check_for_measurements = False
                curr_meas_offset = 0

            # Check if observable index matches with the one we want to track
            if int(instruction.gate_args_copy()[0]) == observable_index:
                # Append all measurement rec values
                for target in instruction.targets_copy():
                    measurement_records.append(target.value)

            # Enable lookout for measurements
            check_for_measurements = True

        # Check for Measurements
        if instruction.name in {"M", "MX", "MZ", "MY"} and check_for_measurements:
            # Add to offset
            curr_meas_offset += len(instruction.targets_copy())
            measurement_in_last_tick = True

    return measurement_records
