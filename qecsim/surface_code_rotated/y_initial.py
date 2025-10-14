import stim
from .data_models import Config, Patch, Context, NoiseModel, CircuitResult

Coord = complex

__all__ = ["y_initial"]

def y_initial(*, 
              lct: Context, 
              patch: dict[str, Patch], 
              cfg: Config,
              noise: NoiseModel,
              offset: complex = 0 + 0j) -> CircuitResult:
    
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    #-Retrieving Global Infomration
    q2i = lct.q2i
    i2q = lct.i2q
    stab_to_data = lct.stab_to_data

    #-Retrieving Data Coords
    data = patch.data

    #-Retrieving Index from Stabilizers of the Lattices
    x_stab_index = patch.x_stab
    z_stab_index = patch.z_stab

    """
    For the Y basis initilization we exclude the upper right data qubit
    -> One weight 4 stabilizer reduced to weight 3 and boundary stabilizer missing
    """

    # Finding Upper right qubit index -> need to look in 2-CX
    y_index = q2i[1 + 1j + offset]

    ########################
    # Define Initial Circuit
    ########################

    initial_circuit = stim.Circuit()

    initial_circuit.append("TICK")

    #1) Reset/ Basis
    initial_circuit.append("H", x_stab_index)
    initial_circuit.append("TICK")

    #2) CX Operations

    def _pairs_for(order: str) -> list[list[int]]:
        # Return the Pairs needed at the current order
        return [[q2i[cp[1]], q2i[cp[0]]] for cp, o in stab_to_data.items() if o == order]

    def _append_by_order(op: str, order: str, noise: float = 0.0) -> None:
        # Getting pair info
        for pair in _pairs_for(order):
            # Checking for upper corner CX and leave it out!
            if y_index not in pair:
                #Adding Pair on Operation
                if op == "CX":
                    initial_circuit.append(op, pair)

    # Adding all the CX gates
    for order in ("1-CX", "2-CX", "3-CX", "4-CX"):
        _append_by_order("CX", order)
        initial_circuit.append("TICK")

    #-------Continue-Circuit------------

    #3) Basis/ Measurement
    initial_circuit.append("H", x_stab_index)

    #-------Continue-Circuit------------

    initial_circuit.append("TICK")

    #-------Continue-Circuit------------
    
    initial_circuit.append("M", x_stab_index + z_stab_index)

    return CircuitResult(circuit=initial_circuit)
