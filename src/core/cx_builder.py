from collections.abc import Iterable

import stim

__all__ = ["cx_builder"]


def cx_builder(
    *,
    q2i: dict,
    stab_to_data: dict,
    circuit: stim.Circuit,
    orders: Iterable[str] = ("1-CX", "2-CX", "3-CX", "4-CX"),
    excluded_index: int | None = None,
    add_operator_before: dict[str, list[tuple[str, list[int]]]] = None,
    add_operator_after: dict[str, list[tuple[str, list[int]]]] = None,
) -> None:
    """
    Used to build the CX Gates into the Circuit

    Args:
        q2i: A mapping from Qubit to Index
        stab_to_data: A mapping from Stabilizer to Data Qubit and Order
        circuit: The Circuit to build into
        orders: What strings in stab_to_data to apply in what order
        excluded_index: An index to exclude from CX application
        noise_overwrite: Whether to apply noise after CX gates
        noise: The Noise Model to use
        add_operator_before: Operators to add before the CX gates
        -> Dict with order as key and list of (operator, indices) tuples as value
        -> Example: {"2-CX": [("H", [1,2,3]), ("S_DAG", [1,2,3])]}
        add_operator_after: Operators to add after the CX gates
        -> Dict with order as key and list of (operator, indices) tuples as value
        -> Example: {"2-CX": [("S", [1,2,3]), ("H", [1,2,3])]}

    Returns:
        None
    """

    def _pairs_for(order: str) -> list[list[int]]:
        # Return the Pairs needed at the current order
        return [[q2i[cp[1]], q2i[cp[0]]] for cp, o in stab_to_data.items() if o == order]

    def _append_by_order(
        op: str,
        order: str,
    ) -> None:
        # Getting pair info
        for pair in _pairs_for(order):
            # Check for excluded index
            if excluded_index not in pair:
                # Adding Pair on Operation
                if op == "CX":
                    circuit.append(op, pair)
                elif op == "CZ":
                    circuit.append(op, pair)
                else:
                    raise ValueError("Unsupported Operator")

    def _append_by_single_index(
        op: str,
        index: int,
    ) -> None:
        # Adding Single Qubit Operation
        if op in ["X", "Y", "Z", "H", "S", "S_DAG", "T", "T_DAG"]:
            circuit.append(op, index)
        else:
            raise ValueError("Unsupported Single Qubit Operator")

    # Adding all the CX gates
    for order in orders:
        # Add before CX operations
        if add_operator_before is not None and order in add_operator_before:
            for op, indices in add_operator_before[order]:
                _append_by_single_index(op, indices)
            circuit.append("TICK")

        # Check if CX or CZ
        if order.endswith("CX"):
            # Add CX Operations
            _append_by_order("CX", order)
        elif order.endswith("CZ"):
            # Add CZ Operations
            _append_by_order("CZ", order)
        else:
            raise ValueError("Unsupported Multi Qubit Operator inside Order")

        # Add after CX operations
        if add_operator_after is not None and order in add_operator_after:
            circuit.append("TICK")
            for op, indices in add_operator_after[order]:
                _append_by_single_index(op, indices)

        circuit.append("TICK")

    return None
