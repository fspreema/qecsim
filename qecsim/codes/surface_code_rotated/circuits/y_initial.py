import stim

from qecsim.core.cx_builder import cx_builder
from qecsim.core.data_models import (
    CircuitResult,
    Context,
    Patch,
)

Coord = complex

__all__ = ["y_initial"]


def y_initial(
    *,
    lct: Context,
    patch: dict[str, Patch],
    offset: complex = 0 + 0j,
) -> CircuitResult:
    #################################################
    # Exporting all necessary values from Dataclasses
    #################################################

    # -Retrieving Global Infomration
    q2i = lct.q2i
    stab_to_data = lct.stab_to_data

    # -Retrieving Index from Stabilizers of the Lattices
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

    # 1) Reset/ Basis
    initial_circuit.append("H", x_stab_index)
    initial_circuit.append("TICK")

    # 2) CX Operations
    cx_builder(
        q2i=q2i,
        stab_to_data=stab_to_data,
        circuit=initial_circuit,
        excluded_index=y_index,
    )

    # 3) Basis/ Measurement
    initial_circuit.append("H", x_stab_index)
    initial_circuit.append("TICK")
    initial_circuit.append("M", x_stab_index + z_stab_index)

    return CircuitResult(circuit=initial_circuit)
