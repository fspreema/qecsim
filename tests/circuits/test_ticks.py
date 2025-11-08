from qecsim.codes.lattice_surgery.builder import surgery_circuit
from qecsim.core.circuit_utils import has_consecutive_ticks


# X Flows alalogue to Z flow
def test_no_double_ticks_basic_x_flow():
    c = surgery_circuit(
        distance=3,
        target_state_init="X+",
        control_state_init="X+",
        flow_observable="X -> X",
    )

    assert not has_consecutive_ticks(c)


# Y flow initializes Y Basis
def test_no_double_ticks_with_y_flow():
    c = surgery_circuit(
        distance=3,
        target_state_init="Z0",
        control_state_init="Y+",
        flow_observable="YZ -> XY",
    )

    assert not has_consecutive_ticks(c)
