__all__ = ["MeasurementTracker"]


class MeasurementTracker:
    def __init__(self):
        """
        Initializes the MeasurementTracker with no measurements and an empty tag dictionary.
        """

        self.total_measurements = 0
        # Tags store: {"tag_name": [indices of measurements associated with this tag]}
        self.tags: dict[str, list[int]] = {}

    def add_previous_measurements(self, count: int):
        """
        Adds a specified number of previous measurements to the total count.

        Args:
            count (int): Number of previous measurements to add.
        """

        self.total_measurements += count

    def add_measurements(
        self,
        measured_qubits: list[int],
        specific_qubits: list[int] = None,
        tag: str = None,
    ):
        """
        Adds measurements for the given qubits to the tracker and optionally tags them.

        Args:
            measured_qubits (list): List of qubit indices to be measured.
            specific_qubits (list, optional): Specific qubits out of the measured qubits which
            should be tracked.
                -> If None, all measured qubits are tracked.
            specific_round (int, optional): Specific round in which the specific_qubits selected
            or the whole meeasurement block should be tracked.
                -> Used for measurements inside a stim repeat block.
            repeats (int, optional): Number of times the measurements are repeated.
                -> Needed if Measurement is inside a stim repeat block.
            tag (str, optional): Tag to associate with these measurements.
        """

        if tag:
            if specific_qubits is None:
                # Track all measured qubits: their positions are simply the full range.
                tracked_positions = [
                    self.total_measurements + i + 1 for i in range(len(measured_qubits))
                ]
            else:
                # Map the specified qubit IDs to their positions within measured_qubits
                tracked_positions = []
                for pos, qubit in enumerate(measured_qubits):
                    if qubit in specific_qubits:
                        tracked_positions.append(self.total_measurements + pos + 1)

            # Apply offset to all tracked positions
            final_positions = [pos for pos in tracked_positions]

            self.tags[tag] = final_positions

        self.total_measurements += len(measured_qubits)

    def get_tagged_measurements(self, tag: str) -> list[int]:
        """
        Retrieves the measurement indices associated with the given tag.

        Args:
            tag (str): The tag for which to retrieve measurement indices.

        Returns:
            list[int]: List of measurement indices associated with the tag.
        """

        if tag not in self.tags:
            raise ValueError(f"Tag '{tag}' not found in MeasurementTracker.")

        indices = self.tags[tag]

        # Generate list of measurement indices for the tag
        return [-self.total_measurements - 1 + i for i in indices]
