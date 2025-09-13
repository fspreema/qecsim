from typing import Dict, Tuple, Mapping
import stim
from dataclasses import dataclass
from .dataclasses import Config, Patch, Context

Coord = complex

__all__ = ["reset"]

def reset(*, lct : Context, patches: dict[str, Patch], 
            cfg : Config) -> stim.Circuit:
    
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
    stab_to_data = lct.stab_to_data

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

    elif init_state in {"+", "-"}:
        reset_circuit.append("RX", data)

        # ancilla are still init in 0
        reset_circuit.append("R", x_stab_index + z_stab_index)

        #Create logical X Data String:
        data_log = []

        for real in range(1,(distance * 2),2):
            data_log.append(q2i[real + 3j])

        if init_state == "-":
            reset_circuit.append("Z", data_log)

    else:
        ValueError("Not a valid init Basis")

    return reset_circuit