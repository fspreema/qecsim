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
    total_measurements = 0

    # Loop through all instructions in the circuit and check for OBSERVABLE_INCLUDE instructions
    for instruction in circuit:
        if instruction.name == "OBSERVABLE_INCLUDE":
            # Check if observable index matches with the one we want to track
            if int(instruction.gate_args_copy()[0]) == observable_index:
                # Get absolute measurement indices
                for target in instruction.targets_copy():
                    # Check if observable tracks records and not indices -> Else skip
                    if not target.is_measurement_record_target:
                        continue

                    # absolute_index = current_total_measurements + target_value
                    measurement_records.append(total_measurements + target.value)

        # Check for Measurements to update total_measurements
        if instruction.name in {"M", "MX", "MZ", "MY"}:
            total_measurements += len(instruction.targets_copy())
        elif instruction.name in {"DETECTOR", "OBSERVABLE_INCLUDE"}:
            # These don't increase measurement count
            pass

    # If circuit completly ran through, convert absolute indices back to relative
    # relative_index = absolute_index - total_measurements
    final_recs = [idx - total_measurements for idx in measurement_records]

    return final_recs
