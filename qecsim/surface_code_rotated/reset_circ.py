import stim
import numpy as np
from .data_models import Config, Patch, Context, CircuitResult

Coord = complex

__all__ = ["reset"]

def reset(*, 
        lct: Context, 
        patches: dict[str, Patch], 
        cfg: Config,
        logical_h: bool) -> CircuitResult:
    
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Loading in Patches
    patch = patches["patch"]

    #-Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    distance = cfg.distance
    init_state = cfg.state_init
    log_obs = cfg.obs

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab

    ########################
    # Define Initial Circuit
    ########################

    reset_circuit = stim.Circuit()

    #-----BUILDING-INITILIZATION-CIRCUIT----

    #Appending Coords
    for q, i in q2i.items():
        reset_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

    #Appending Reset
    if init_state in {"0", "1"}:
        reset_circuit.append("RZ", data + x_stab_index + z_stab_index)

        #Create logical X Data String:
        data_log = []

        for imag in range(1,(distance * 2),2):
            data_log.append(q2i[3 + imag * 1j])

        if init_state == "1":
            reset_circuit.append("X", data_log)

        if log_obs == "X":

            if not logical_h:

                """
                We need to remove the added Pauli measurement from the end of the circuit
                -> Else the X paulis tring would anticommute with the RZ reset of the data
                """

                # Getting corresponding logical string and rec
                log_x = []

                for imag in range(1, (distance * 2), 2):
                    log_x.append(q2i[1 + imag*1j])

                # XORing the observable away
                reset_circuit.append("OBSERVABLE_INCLUDE", [f"X{index}" for index in log_x], 0)

    elif init_state in {"+", "-"}:
        reset_circuit.append("RX", data)

        # ancilla are still init in 0
        reset_circuit.append("RZ", x_stab_index + z_stab_index)

        #Create logical X Data String:
        data_log = []

        for real in range(1,(distance * 2),2):
            data_log.append(q2i[real + 3j])

        if init_state == "-":
            reset_circuit.append("Z", data_log)

        if log_obs == "Z":

            if not logical_h:

                """
                We need to remove the added Pauli measurement from the end of the circuit
                -> Else the Z paulis tring would anticommute with the RZ reset of the data
                """

                # Getting corresponding logical string and rec
                log_z = []

                for real in range(1, (distance * 2), 2):
                    log_z.append(q2i[real + 1j])

                # XORing the observable away
                reset_circuit.append("OBSERVABLE_INCLUDE", [f"Z{index}" for index in log_z], 0)

    elif init_state in {"+i", "-i"}:

        """
        Look at Crumble circuit for a better understanding
        -> Half Half initlization of x and z basis
        """

        data_rx = []
        data_rz = []

        xs = [i2q[i].real for i in data]
        ys = [i2q[i].imag for i in data]

        # Calc threshold for diagonal cut        
        s0 = (min(xs)+max(xs))/2 + (min(ys)+max(ys))/2

        skip_coord = 1 + 1j

        for data_index in data:
            c = i2q[data_index]
            if c == skip_coord:
                continue
            
            # Diagonal Cut
            if (c.real + c.imag) >= s0:
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