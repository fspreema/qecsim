import stim

from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry


class YBasisGetCircuitFlows:
    def __init__(
        self,
        geometry: SurgeryGeometry,
    ) -> stim.Circuit:
        """
        For a given circuit without final measurement and reset, the needed Measurements can be
        deterermined by the flow type.

        -> These are given as a list of qubit indicies to be measured.
        -> These are then applied to the circuit as additional measurement
        records to the same observable

        Returns:
            stim.Circuit: Corrected Circuit with the correct logical observable included
        """

        # Determine if the flow exists and if so what measurements are needed
        self.return_circuit = stim.Circuit()
        self.geometry = geometry

    def get_flows(
        self,
        flow_circ: stim.Circuit,
        circuits_between: list[stim.Circuit],
    ) -> stim.Circuit:
        # Getting logical y string
        x_idx, y_idx, z_idx = self.geometry.get_logical_observables(
            "Y",
            fixed_coord=self.geometry.distance * 2 - 1,
        )

        # Adding logical z string
        logical_xyz_string = "*".join(
            [f"Z{idz}" for idz in z_idx] + [f"Y{y_idx[0]}"] + [f"X{idx}" for idx in x_idx],
        )

        logical_creation = f"{1} -> {logical_xyz_string}"
        logical_contraction = f"{logical_xyz_string} -> {1}"

        (logical_creation_rec,) = logical_creation_circ.solve_flow_measurements(
            [stim.Flow(logical_creation)],
        )
        (logical_contraction_rec,) = logical_contraction_circ.solve_flow_measurements(
            [stim.Flow(logical_contraction)],
        )

        # Calculating target rec pos
        rec_pos = []

        contraction_records = logical_contraction_circ.num_measurements
        creation_records = (
            logical_contraction_circ.num_measurements
            + logical_creation_circ.num_measurements
            + sum(circuit.num_measurements for circuit in circuits_between)
        )

        # Adding the Observable
        for index_creation in logical_creation_rec:
            current_rec_crea = creation_records - index_creation
            rec_pos.append(-current_rec_crea)

        for index_contraction in logical_contraction_rec:
            current_rec_cont = contraction_records - index_contraction
            rec_pos.append(-current_rec_cont)

        # Adding measurements to the logical observable
        self.return_circuit.append("OBSERVABLE_INCLUDE", [stim.target_rec(k) for k in rec_pos], 0)
        self.return_circuit.append("TICK")

        return self.return_circuit
