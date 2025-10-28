__all__ = ["get_logical_strings"]


def get_logical_strings(
    q2i: dict[complex, int],
    distance: int,
    *,
    shift_cx_for_y: bool = False,
    shift_tz_for_y: bool = False,
    shift_tx_for_y: bool = False,
    shift_cz_for_y: bool = False,
) -> dict[str, list[int]]:
    """
    Args:
        q2i: Mapping from complex coordinates to qubit indices
        distance: Code distance

    Returns:
        Dictionary with keys:
            - 'a_z': Ancilla Z-logical indices (horizontal string at y=1)
            - 't_x': Target X-logical indices (vertical string at x=distance*2+1)
            - 't_z': Target Z-logical indices (horizontal string at y=1, offset region)
            - 't_y': Target Y-logical indices (Z-string + Y-corner + X-string)
            - 'c_x': Control X-logical indices (vertical string at x=1, offset region)
            - 'c_z': Control Z-logical indices (horizontal string at y=distance*2+1)
            - 'c_y': Control Y-logical indices (Z-string + Y-corner + X-string)
    """
    # Control Y logical observable components
    # Y is at bottom right corner of control region: (distance*2-1, distance*4-1)
    c_y_z_string = [q2i[real + (distance * 4 - 1) * 1j] for real in range(1, distance * 2 - 1, 2)]
    c_y_corner = [q2i[distance * 2 - 1 + (distance * 4 - 1) * 1j]]
    c_y_x_string = [
        q2i[(distance * 2 - 1) + imag * 1j] for imag in range(distance * 2 + 1, distance * 4 - 1, 2)
    ]

    # Target Y logical observable components
    # Y is at bottom right corner of target region: (distance*4-1, distance*2-1)
    t_y_z_string = [
        q2i[real + (distance * 2 - 1) * 1j] for real in range(distance * 2 + 1, distance * 4 - 1, 2)
    ]
    t_y_corner = [q2i[distance * 4 - 1 + (distance * 2 - 1) * 1j]]
    t_y_x_string = [q2i[(distance * 4 - 1) + imag * 1j] for imag in range(1, distance * 2 - 1, 2)]

    # Offsets for conditional shifts in Y-including flows (symmetric for both directions)
    cx_shift_real = (distance * 2 - 1) if shift_cx_for_y else 1
    tz_shift_imag = (distance * 2 - 1) if shift_tz_for_y else 1
    tx_shift_real = (distance * 4 - 1) if shift_tx_for_y else (distance * 2 + 1)
    cz_shift_imag = (distance * 4 - 1) if shift_cz_for_y else (distance * 2 + 1)

    return {
        "a_z": [q2i[real + 1j] for real in range(1, distance * 2, 2)],
        "t_x": [q2i[tx_shift_real + imag * 1j] for imag in range(1, distance * 2, 2)],
        "t_z": [
            q2i[real + tz_shift_imag * 1j] for real in range(distance * 2 + 1, distance * 4, 2)
        ],
        "t_y": {
            "z_string": t_y_z_string,
            "y_corner": t_y_corner,
            "x_string": t_y_x_string,
        },
        "c_x": [
            q2i[cx_shift_real + imag * 1j] for imag in range(distance * 2 + 1, distance * 4, 2)
        ],
        "c_z": [q2i[real + cz_shift_imag * 1j] for real in range(1, distance * 2, 2)],
        "c_y": {
            "z_string": c_y_z_string,
            "y_corner": c_y_corner,
            "x_string": c_y_x_string,
        },
    }
