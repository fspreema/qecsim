import stim

__all__ = ["MeasurementTracker"]


class MeasurementTracker:
    def __init__(self):
        """
        Initializes the MeasurementTracker with no measurements and an empty tag dictionary.
        """

        self.total_measurements = 0
        # Tags store: {"tag_name": [indices of measurements associated with this tag]}
        self.tags: dict[str, list[int]] = {}
        # Dictionary to store the Detector information for each qubit index
        self.detector_dict: dict[str, list[tuple[int, int]]] = {
            "Ancilla": [],
            "Control_&_Target": [],
            "Merge_AC": [],
            "Split_AC": [],
            "Merge_AT": [],
            "Split_AT": [],
            "Merge_AC_untouched": [],
            "Split_AC_untouched": [],
            "Merge_AT_untouched": [],
            "Split_AT_untouched": [],
        }
        # Dict to store the complete measurement history
        self.full_measurement_history: dict[str, list[tuple[int, int]]] = {
            "Ancilla": [],
            "Control_&_Target": [],
            "Merge_AC": [],
            "Split_AC": [],
            "Merge_AT": [],
            "Split_AT": [],
            "Merge_AC_untouched": [],
            "Split_AC_untouched": [],
            "Merge_AT_untouched": [],
            "Split_AT_untouched": [],
        }
        # Track the most recent measurement batch to avoid double-counting
        self._last_measured_qubits: list[int] | None = None
        self._last_batch_start: int | None = None
        self._last_batch_counted: bool = False

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

        batch_start = self.total_measurements

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
        self._last_measured_qubits = list(measured_qubits)
        self._last_batch_start = batch_start
        self._last_batch_counted = True

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

    def add_measurements_to_detector_dict(
        self,
        patch_type: str,
        measured_qubits: list[int],
        selected_qubits: list[int] | None = None,
    ):
        """
        Adds all the measurements in "selected qubits" into the dictionary

        Dict Functionality:
            This Dict has the following information.
            -> Key: Patch_Type (i.e. Ancilla, Target, Control)
            -> Entries: Qubit Index, Absolute Index of the Measurement

        Function:
            -> Creates Dict (Detector_dict and Full Measurement History Dict)
            -> Detector Dict only has the most recent measurement which was selected with
                "selected_qubits" inside this method
                -> Gets cleared every time this method is called again for the same patch type
            -> Full Measurement History Dict has all measurements of the patch type, even if
                they were not selected with "selected_qubits"
                -> This is needed to check for previous measurements of the same qubit index in case
                there is only one entry in the Detector Dict for a specific qubit index
                -> This Dict includes only the last 2 recent measurements for each qubit index,
                meaning if there are more than 2 measurements of the same qubit index,
                the oldest one gets deleted
        """

        # Check Validity of Method Arguemnts
        if patch_type not in {
            "Ancilla",
            "Control_&_Target",
            "Merge_AC",
            "Split_AC",
            "Merge_AT",
            "Split_AT",
            "Merge_AC_untouched",
            "Split_AC_untouched",
            "Merge_AT_untouched",
            "Split_AT_untouched",
        }:
            raise ValueError(
                f"Invalid patch type '{patch_type}'. Expected 'Ancilla', "
                "'Control_&_Target', 'Merge_AC', 'Split_AC', "
                "'Merge_AT', 'Split_AT', 'Merge_AC_untouched', 'Split_AC_untouched', "
                "'Merge_AT_untouched', or 'Split_AT_untouched'.",
            )

        if selected_qubits is None:
            selected_qubits = measured_qubits

        num_new_measurements = len(measured_qubits)

        # If the same batch was already counted via add_measurements, reuse its start index.
        already_counted = (
            self._last_batch_counted
            and self._last_measured_qubits == measured_qubits
            and self._last_batch_start is not None
            and self._last_batch_start + num_new_measurements == self.total_measurements
        )
        batch_start = self._last_batch_start if already_counted else self.total_measurements

        # Clear the current round's detectors for the patch type
        self.detector_dict[patch_type] = []

        ###################
        # REGULAR DETECTORS
        ###################

        # Enumerating current measured qubit batch and their positions
        for offset, qubit in enumerate(measured_qubits):
            # Calc absolute index of the current measurement
            abs_idx = batch_start + offset

            # If Full Measurement History Dict already has 2 entries on the same qubit index,
            # remove the oldest one (the one with the smallest absolute index)
            existing_entries_full = [
                entry for entry in self.full_measurement_history[patch_type] if entry[0] == qubit
            ]
            if len(existing_entries_full) >= 2:
                # Find the entry with the smallest absolute index (oldest measurement)
                oldest_entry = min(existing_entries_full, key=lambda x: x[1])
                self.full_measurement_history[patch_type].remove(oldest_entry)

            # Add Entry to Full Dict if qubit is in selected_qubits
            self.full_measurement_history[patch_type].append((qubit, abs_idx))

            # Add Entry to Detector Dict if qubit is in selected_qubits
            if qubit in selected_qubits:
                self.detector_dict[patch_type].append((qubit, abs_idx))

        ###################
        # MERGING DETECTORS
        ###################

        #####################
        # SPLITTING DETECTORS
        #####################

        # Update total measurements only if this batch was not already counted
        if not already_counted:
            self.total_measurements += num_new_measurements
            self._last_measured_qubits = list(measured_qubits)
            self._last_batch_start = batch_start
            self._last_batch_counted = True

    def get_records_for_detectors(self, patch_type: str) -> list[list[stim.GateTarget]]:
        """
        *Return:
            Returns List of Lists with measurements of the same qubit index on different
            time steps
                -> I.e. [[-1,-20], [-2,-65], ...]

            Be aware: These Int Records are already in the stim.target_rec format!

        *Functionality:
            For the current Patch get Info of Detector Dictionary and FUll Measurement Histroy Dict
            -> If the Dcit Entrie only has one entry, double check for previous measurements of the
            same qubit index in the full history dict and add this to the records list
            -> If the Dict Entry has two entries, simply add these to the records list
        """

        # Check Validity of Method Arguemnts
        if patch_type not in {
            "Ancilla",
            "Control_&_Target",
            "Merge_AC",
            "Split_AC",
            "Merge_AT",
            "Split_AT",
            "Merge_AC_untouched",
            "Split_AC_untouched",
            "Merge_AT_untouched",
            "Split_AT_untouched",
        }:
            raise ValueError(
                f"Invalid patch type '{patch_type}'. Expected 'Ancilla', "
                "'Control_&_Target', 'Merge_AC', 'Split_AC', "
                "'Merge_AT', 'Split_AT', 'Merge_AC_untouched', 'Split_AC_untouched', "
                "'Merge_AT_untouched', or 'Split_AT_untouched'.",
            )

        records: list[list[stim.GateTarget]] = []

        for qubit, curr_abs_idx in self.detector_dict[patch_type]:
            # Find all entries in the full history dict for the same qubit index
            full_history_entries = [
                entry for entry in self.full_measurement_history[patch_type] if entry[0] == qubit
            ]

            # Calculate the Relative Index for the current measurement (for stim)
            current_rel_idx = curr_abs_idx - self.total_measurements

            # If there is only one entry in the full history dict, add the current rec just once
            if len(full_history_entries) == 1:
                records.append([stim.target_rec(current_rel_idx)])
            # If there are two entries, add both to the records list
            elif len(full_history_entries) == 2:
                # Calculate the Relative Index for the previous measurement (for stim)
                previous_rel_idx = full_history_entries[0][1] - self.total_measurements
                records.append(
                    [stim.target_rec(current_rel_idx), stim.target_rec(previous_rel_idx)],
                )
            else:
                raise ValueError(
                    f"Unexpected number of entries for qubit {qubit} in full measurement history.",
                )

        return records
