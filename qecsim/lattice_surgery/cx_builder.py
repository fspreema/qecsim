import stim
from typing import Iterable
from .dataclasses import NoiseModel

__all__ = ["cx_builder"]

def cx_builder(*,
               q2i: dict,
               stab_to_data: dict,
               circuit: stim.Circuit,
               orders: Iterable[str] = ("1-CX", "2-CX", "3-CX", "4-CX"),
               noise: NoiseModel) -> None:

    """
    Used to build the CX Gates into the Circuit
    """

    def _pairs_for(order: str) -> list[list[int]]:
            # Return the Pairs needed at the current order
            return [[q2i[cp[1]], q2i[cp[0]]] for cp, o in stab_to_data.items() if o == order]

    def _append_by_order(op: str, order: str, noise: float = 0.0) -> None:
        # Getting pair info
        for pair in _pairs_for(order):
            #Adding Pair on Operation
            if op == "CX":
                circuit.append(op, pair)
            elif op == "DEPOLARIZE2":
                circuit.append(op, pair, noise)
            else:
                 raise ValueError("Unsupported Operator")
            
    # Adding all the CX gates
    for order in orders:
        _append_by_order("CX", order)
        if noise.after_c_depol_prob > 0:
            _append_by_order("DEPOLARIZE2", order, noise.after_c_depol_prob)
        circuit.append("TICK")

