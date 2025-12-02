import stim


def normalize_ticks(circ: stim.Circuit) -> stim.Circuit:
    """
    Return a new `stim.Circuit` instance with consecutive TICK instructions in `circ` collapsed
    to a single TICK. This is a cleanup that removes empty time slices.

    Notes:
    - Only checks outside of repeat blocks!
    - The input circuit is not modified; a new circuit is returned.
    """

    out_circ = stim.Circuit()
    prev_tick = False

    for inst in circ:
        # Check for Ticks and only add single TICKs back to output
        if inst.name == "TICK":
            if prev_tick:
                continue
            prev_tick = True
            out_circ.append(inst)
        # Reset Marker and add arbitrary instructions
        else:
            prev_tick = False
            out_circ.append(inst)

    return out_circ


def has_consecutive_ticks(circ: stim.Circuit) -> bool:
    """
    Utility predicate to detect whether a circuit contains back-to-back TICKs.
    """

    last_tick = False

    for inst in circ:
        # Check for double ticks and set marker accordingly
        if inst.name == "TICK":
            if last_tick:
                return True
            last_tick = True

        else:
            # If Operation found, set Tick marker back to False
            last_tick = False
    return False
