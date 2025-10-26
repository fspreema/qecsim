import stim

from qecsim.core.data_models import CircuitResult, ConfigSurface as Config, Context, Patch

Coord = complex

__all__ = ["reset"]


def reset(
    *,
    lct: Context,
    patches: dict[str, Patch],
    cfg: Config,
    logical_h: bool = False,
    skip_coords: bool = False,
    offset: complex = 0 + 0j,
) -> CircuitResult:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # -Loading in Patches
    patch = patches["patch"]

    # -Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    distance = cfg.distance
    init_state = cfg.state_init
    log_obs = cfg.obs

    # -Retrieving Data Coords
    data = patch.data

    # -Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab

    ##################
    # Helper Fuctions
    ##################

    def _logical_x_indices() -> list[int]:
        # Vertical string at x=1 (odd grid), along imag axis
        return [q2i[1 + imag * 1j] for imag in range(1, 2 * distance, 2)]

    def _logical_z_indices() -> list[int]:
        # Horizontal string at y=1, along real axis
        return [q2i[real + 1j] for real in range(1, 2 * distance, 2)]

    def _logical_y_indices() -> list[int]:
        # Both
        z_string = [q2i[real + 1j] for real in range(3, 2 * distance, 2)]
        x_string = [q2i[1 + imag * 1j] for imag in range(3, 2 * distance, 2)]
        y_string = [q2i[1 + 1j]]

        return (x_string, y_string, z_string)

    ########################
    # Define Initial Circuit
    ########################

    reset_circuit = stim.Circuit()

    # -----BUILDING-INITILIZATION-CIRCUIT----

    # Appending Coords
    if not skip_coords:
        for q, i in q2i.items():
            reset_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

    #########################
    # Appending Resets for 0/1
    #########################

    if init_state in {"0", "1"}:
        reset_circuit.append("RZ", data + x_stab_index + z_stab_index)

        # Create logical X Data String:
        data_log = []

        for imag in range(1, (distance * 2), 2):
            data_log.append(q2i[3 + imag * 1j])

        if init_state == "1":
            reset_circuit.append("X", data_log)

        if log_obs == "X":
            """
            We need to remove the added Pauli measurement from the end of the circuit
            -> Else the X paulis tring would anticommute with the RZ reset of the data
            """

            if not logical_h:
                # Getting corresponding logical string and rec
                log_x = _logical_x_indices()

                # XORing the observable away
                reset_circuit.append("OBSERVABLE_INCLUDE", [f"X{index}" for index in log_x], 0)

        if log_obs == "Y":
            # XORing the observable away
            reset_circuit.append(
                "OBSERVABLE_INCLUDE",
                [f"X{index}" for index in _logical_y_indices()[0]]
                + [f"Y{index}" for index in _logical_y_indices()[1]]
                + [f"Z{index}" for index in _logical_y_indices()[2]],
                0,
            )

    #########################
    # Appending Resets for +/-
    #########################

    elif init_state in {"+", "-"}:
        reset_circuit.append("RX", data)

        # ancilla are still init in 0
        reset_circuit.append("RZ", x_stab_index + z_stab_index)

        # Create logical X Data String:
        data_log = []

        for real in range(1, (distance * 2), 2):
            data_log.append(q2i[real + 3j])

        if init_state == "-":
            reset_circuit.append("Z", data_log)

        if log_obs == "Z":
            """
            We need to remove the added Pauli measurement from the end of the circuit 
            -> Else the Z paulis tring would anticommute with the RZ reset of the data
            """

            if not logical_h:
                # Getting corresponding logical string and rec
                log_z = _logical_z_indices()

                # XORing the observable away
                reset_circuit.append("OBSERVABLE_INCLUDE", [f"Z{index}" for index in log_z], 0)

        if log_obs == "Y":
            # XORing the observable away
            reset_circuit.append(
                "OBSERVABLE_INCLUDE",
                [f"X{index}" for index in _logical_y_indices()[0]]
                + [f"Y{index}" for index in _logical_y_indices()[1]]
                + [f"Z{index}" for index in _logical_y_indices()[2]],
                0,
            )

    ###########################
    # Appending Resets for +i/-i
    ###########################

    elif init_state in {"+i", "-i"}:
        """
        Look at Crumble circuit for a better understanding
        -> Half Half initlization of x and z basis
        """

        data_rx = []
        data_rz = []

        xs = [i2q[i].real - offset.real for i in data]
        ys = [i2q[i].imag - offset.imag for i in data]

        # Calc threshold for diagonal cut
        s0 = (min(xs) + max(xs)) / 2 + (min(ys) + max(ys)) / 2

        skip_coord = 1 + 1j + offset

        for data_index in data:
            c = i2q[data_index]
            if c == skip_coord:
                continue

            # Diagonal Cut
            if (c.real - offset.real + c.imag - offset.imag) >= s0:
                data_rz.append(q2i[c])
            else:
                data_rx.append(q2i[c])

        reset_circuit.append("RX", data_rx)
        reset_circuit.append("RZ", data_rz)

        # ancilla are still init in 0
        reset_circuit.append("RZ", x_stab_index + z_stab_index)

    else:
        ValueError("Not a valid init Basis")

    return CircuitResult(circuit=reset_circuit)
