import stim

from src.codes.surface_code_rotated.data_geometry import MasterGeometry

Coord = complex

__all__ = ["SurfaceReset"]


class SurfaceReset:
    def __init__(
        self,
        master_geometry: MasterGeometry,
        type: str,
    ):
        """
        Initializes the Surface Reset Circuit.

        Args:
            master_geometry (MasterGeometry)
            type (str): Type of reset circuit to build 
            -> If non_ft Y obs added by reset XX Y ZZ for d=3 f.ex.
        """

        if type not in {"standard", "y_basis", "non_ft_init"}:
            raise ValueError(
                f"Invalid type '{type}' for SurfaceReset. "
                f"Must be 'standard', 'y_basis', or 'non_ft_init'.",
            )

        # Initialize Geometry depending on the type
        self.type = type

        if self.type in {"standard", "non_ft_init"}:
            self.geometry = master_geometry.geometry_std
        elif self.type == "y_basis":
            self.geometry = master_geometry.geometry_ybasis

    def build_circuit(self) -> stim.Circuit:
        circuit = stim.Circuit()

        # Adding Qubit Coords
        circuit += self._adding_qubit_coords()

        # Adding Resets
        circuit += self._adding_resets()

        # Adding Logical Observables if non deterministic basis/ observable
        circuit += self._adding_logical_observables()

        return circuit

    def _adding_qubit_coords(self):
        # Init Coords
        coords_circuit = stim.Circuit()
        # Appending Coords
        for q, i in self.geometry.q2i.items():
            coords_circuit.append("QUBIT_COORDS", [i], [q.real, q.imag])

        return coords_circuit

    def _adding_resets(self):
        # Init Coords
        reset_circuit = stim.Circuit()
        # Depending on init state add resets and logical operator
        if self.geometry.state_init in {"0", "1"}:
            reset_circuit.append("RZ", self.geometry.data_idx)

            if self.geometry.state_init == "1":
                # Create logical X Data String:
                reset_circuit.append("X", self.geometry.get_logical_observables("X"))

        elif self.geometry.state_init in {"+", "-"}:
            reset_circuit.append("RX", self.geometry.data_idx)

            if self.geometry.state_init == "-":
                # Create logical Z Data String:
                reset_circuit.append("Z", self.geometry.get_logical_observables("Z"))

        elif self.geometry.state_init in {"+i", "-i"} and self.type == "y_basis":
            
            reset_circuit.append("RX", self.geometry.data_rx_idx)
            reset_circuit.append("RZ", self.geometry.data_rz_idx)

        elif self.geometry.state_init in {"+i", "-i"} and self.type == "non_ft_init":
            mx_idx, my_idx, mz_idx = self.geometry.get_logical_observables(
                    "Y",
                    fixed_coord=(self.geometry.distance * 2 - 1),
                )
            
            reset_circuit.append("RX", mx_idx)
            reset_circuit.append("RZ", mz_idx)
            reset_circuit.append("RY", my_idx)

            if self.geometry.state_init == "-i":
                # Create logical Y Data String:
                reset_circuit.append("X", self.geometry.get_logical_observables("X"))

        return reset_circuit

    def _adding_logical_observables(self):
        """
        These Observables are added for measurement which in the given
        basis would be non deterministic
        -> We need to remove them in order to simulate the 50/50 outcome

        Inside the Y basis the observables are removed after the switch as the observable
        "lives" only inside the memory rounds of the y circuit
        """
        log_circuit = stim.Circuit()

        if self.geometry.state_init in {"0", "1"}:
            if self.geometry.obs == "X" and self.type == "standard":
                # Getting corresponding logical string and rec
                log_x = self.geometry.get_logical_observables("X")

                # XORing the observable away
                log_circuit.append("OBSERVABLE_INCLUDE", [f"X{index}" for index in log_x], 0)

            if self.geometry.obs == "Y":
                mx_idx, my_idx, mz_idx = self.geometry.get_logical_observables(
                    "Y",
                    fixed_coord=(self.geometry.distance * 2 - 1),
                )

                # XORing the observable away
                log_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [f"X{index}" for index in mx_idx]
                    + [f"Y{index}" for index in my_idx]
                    + [f"Z{index}" for index in mz_idx],
                    0,
                )

        elif self.geometry.state_init in {"+", "-"}:
            if self.geometry.obs == "Z" and self.type == "standard":
                # Getting corresponding logical string and rec
                log_z = self.geometry.get_logical_observables("Z")

                # XORing the observable away
                log_circuit.append("OBSERVABLE_INCLUDE", [f"Z{index}" for index in log_z], 0)

            if self.geometry.obs == "Y":
                mx_idx, my_idx, mz_idx = self.geometry.get_logical_observables(
                    "Y",
                    fixed_coord=(self.geometry.distance * 2 - 1),
                )

                # XORing the observable away
                log_circuit.append(
                    "OBSERVABLE_INCLUDE",
                    [f"X{index}" for index in mx_idx]
                    + [f"Y{index}" for index in my_idx]
                    + [f"Z{index}" for index in mz_idx],
                    0,
                )

        # If NON-FT we do not have an observable defined between the y siwtch and y reverse siwtch,
        # as we do not have these subsections -> Need to remove directly before the reset

        elif self.geometry.state_init in {"+i", "-i"} and self.type == "non_ft_init":
            if self.geometry.obs == "X":
                # Getting corresponding logical string and rec
                log_x = self.geometry.get_logical_observables("X")

                # XORing the observable away
                log_circuit.append("OBSERVABLE_INCLUDE", [f"X{index}" for index in log_x], 0)
            
            elif self.geometry.obs == "Z":
                # Getting corresponding logical string and rec
                log_z = self.geometry.get_logical_observables("Z")

                # XORing the observable away
                log_circuit.append("OBSERVABLE_INCLUDE", [f"Z{index}" for index in log_z], 0)

        return log_circuit
