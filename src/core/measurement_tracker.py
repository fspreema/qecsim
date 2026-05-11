import stim

from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.codes.surface_code_rotated.surface_geom import SurfaceGeometry

__all__ = ["MeasurementTracker"]


class MeasurementTracker:
    VALID_SURGERY_DETECTOR_TYPES = {
        "Ancilla",
        "Control_&_Target",
        "Ancilla_Split_AC",
        "Ancilla_Split_AT",
        "Control_&_Target_Split_AC",
        "Control_&_Target_Split_AT",
        "Merge_AC",
        "Merge_AT",
        "Merge_AC_untouched",
        "Merge_AT_untouched",
        "None",
        "Final_Measurement_Anc",
        "Final_Measurement_C",
        "Final_Measurement_T",
    }

    VALID_SURFACE_DETECTOR_TYPES = {
        "STD_PATCH",
        "None",
    }


    def __init__(self, geometry: SurgeryGeometry | SurfaceGeometry):
        """
        Initializes the MeasurementTracker with no measurements and an empty tag dictionary.
        """
        self.geometry = geometry

        if isinstance(geometry, SurgeryGeometry):
            self.VALID_DETECTOR_TYPES = self.VALID_SURGERY_DETECTOR_TYPES
        elif isinstance(geometry, SurfaceGeometry):
            self.VALID_DETECTOR_TYPES = self.VALID_SURFACE_DETECTOR_TYPES
        else:
            raise ValueError("Invalid geometry type provided to MeasurementTracker.")

        self.total_measurements = 0
        # Tags store: {"tag_name": [indices of measurements associated with this tag]}
        self.tags: dict[str, list[int]] = {}
        # Dictionary to store the Detector information for each qubit index
        self.detector_dict: dict[str, list[tuple[int, int]]] = {
            pt: [] for pt in self.VALID_DETECTOR_TYPES
        }
        # Dictionary to store split measurements that await their partner from the other patch
        self.storage_for_split_matching: dict[int, list[int, int, int]] = {}
        # Dict to store the complete measurement history
        self.full_measurement_history: dict[str, list[tuple[int, int]]] = {
            pt: [] for pt in self.VALID_DETECTOR_TYPES
        }

    def add_previous_measurements(self, count: int):
        """
        Adds a specified number of previous measurements to the total count.

        Args:
            count (int): Number of previous measurements to add.
        """

        self.total_measurements += count

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

    def add_measurements_to_tracker(
        self,
        patch_type: str,
        measured_qubits: list[int],
        qubits_for_detectors: list[int] | None = None,
        tagged_qubits: list[int] | None = None,
        tag: str | None = None,
    ):
        """
        Adds all the measurements in "qubits_for_detectors" into the dictionary as well as
        adds taggs for selected measurements if needed.

        Dict Functionality:
            This Dict has the following information.
            -> Key: Detector Type (i.e. Ancilla, Target, Control)
            -> Entries: Qubit Index, Absolute Index of the Measurement

        Function:
            -> Creates Dict (Detector_dict and Full Measurement History Dict)
            -> Detector Dict only has the most recent measurement which was selected with
                "selected_qubits" inside this method
                -> Gets cleared every time this method is called again for the same detector type
            -> Full Measurement History Dict has all measurements of the detector type, even if
                they were not selected with "selected_qubits"
                -> This is needed to check for previous measurements of the same qubit index in case
                there is only one entry in the Detector Dict for a specific qubit index
                -> This Dict includes only the last 2 recent measurements for each qubit index,
                meaning if there are more than 2 measurements of the same qubit index,
                the oldest one gets deleted

        Args:
            patch_type (str): The type of patch for which the measurements are being added.
                -> If None is selected no detectors will be build based on these measurements!
                -> Soley for tracking the measurements in the tracker and for tagging if needed
            measured_qubits (list[int]): List of qubit indices that were measured in the
                current round.
            qubits_for_detectors (list[int] | None): List of qubit indices from measured_qubits that
                should be included in the detector dictionary. If None, all measured_qubits
                are included.
            tagged_qubits (list[int] | None): List of qubit indices from measured_qubits that should
                be tagged with the provided tag. If None, no qubits will be tagged.
            tag (str | None): The tag to associate with the tagged_qubits. If None, no tagging will
                occur.
        """

        ######################################
        # Check Validity of Method Arguemnts #
        ######################################

        if patch_type not in self.VALID_DETECTOR_TYPES:
            raise ValueError(
                f"Invalid patch type '{patch_type}'. Expected"
                f" one of {self.VALID_DETECTOR_TYPES}.",
            )

        ########################
        # Add Tag if specified #
        ########################

        if tag:
            if tagged_qubits is None:
                raise ValueError("Tagged qubits must be provided if a tag is specified.")
            else:
                # Map the specified qubit IDs to their positions within measured_qubits
                tracked_positions = []
                for pos, qubit in enumerate(measured_qubits):
                    if qubit in tagged_qubits:
                        tracked_positions.append(self.total_measurements + pos + 1)

            # Apply offset to all tracked positions
            final_positions = [pos for pos in tracked_positions]

            self.tags[tag] = final_positions

        ############################################################################
        # Adding Measurements into Detector Dict and Full Measurement History Dict #
        ############################################################################

        if qubits_for_detectors is None:
            qubits_for_detectors = measured_qubits

        # Clear the current round's detectors for the patch type
        self.detector_dict[patch_type] = []

        # Enumerating current measured qubit batch and their positions
        for offset, qubit in enumerate(measured_qubits):
            # Calc absolute index of the current measurement
            abs_idx = self.total_measurements + offset

            # If Full Measurement History Dict already has 2 entries on the same qubit index,
            # remove the oldest one (the one with the smallest absolute index)
            existing_entries_full = [
                entry for entry in self.full_measurement_history[patch_type] if entry[0] == qubit
            ]
            if len(existing_entries_full) >= 2:
                # Find the entry with the smallest absolute index (oldest measurement)
                oldest_entry = min(existing_entries_full, key=lambda x: x[1])
                self.full_measurement_history[patch_type].remove(oldest_entry)

            # Add Entry to Full Dict
            self.full_measurement_history[patch_type].append((qubit, abs_idx))

            # Add Entry to Detector Dict if qubit is in qubits_for_detectors
            if qubit in qubits_for_detectors:
                self.detector_dict[patch_type].append((qubit, abs_idx))

        ################################################################
        # Updating total measurements count after processing the batch #
        ################################################################

        self.total_measurements += len(measured_qubits)

    def get_records_for_detectors(
        self,
        patch_type: str,
    ) -> list[list[stim.GateTarget]]:
        """
        *Return:
            Returns List of Lists with measurements of the same qubit index on different
            time steps
                -> I.e. [[-1,-20], [-2,-65], ...]

            Be aware: These Int Records are already in the stim.target_rec format!

        *Functionality:
            For the current Patch get Info of Detector Dictionary and Full Measurement History Dict
            -> If the Dict Entry only has one entry, double check for previous measurements of the
            same qubit index in the full history dict and add this to the records list
            -> If the Dict Entry has two entries, simply add these to the records list

            For Merging Regions, we build up the detector dict differently
            -> Create Full Meas dict out of the Old still existing measurements
                from Ancilla, Control And Target
                -> Check indices measured in the Merg and compare with indices inside
                the other dict
        """

        ######################################
        # Check Validity of Method Arguemnts #
        ######################################

        if patch_type not in self.VALID_DETECTOR_TYPES:
            raise ValueError(
                f"Invalid patch type '{patch_type}'. Expected"
                f" one of {self.VALID_DETECTOR_TYPES}.",
            )

        #####################
        # REGULAR DETECTORS #
        #####################

        if patch_type in {
            "Ancilla",
            "Control_&_Target",
            "STD_PATCH",
        }:
            records = self._get_regular_records(patch_type)

        #######################
        # SPLITTING DETECTORS #
        #######################

        elif patch_type in {
            "Ancilla_Split_AC",
            "Ancilla_Split_AT",
            "Control_&_Target_Split_AC",
            "Control_&_Target_Split_AT",
        }:
            records = self._get_splitting_records(patch_type)

        #####################
        # MERGING DETECTORS #
        #####################

        elif patch_type in {
            "Merge_AC_untouched",
            "Merge_AT_untouched",
        }:
            records = self._get_untouched_merging_records(patch_type)

        elif patch_type in {
            "Merge_AC",
            "Merge_AT",
        }:
            records = self._get_merging_records(patch_type)

        elif patch_type in {
            "Final_Measurement_Anc",
            "Final_Measurement_C",
            "Final_Measurement_T",
        }:
            records = self._get_final_measurement_records(patch_type)

        elif patch_type == "None":
            records = []

        else:
            raise ValueError(f"Unknown patch type encountered: {patch_type}")

        return records

    def _get_regular_records(self, patch_type: str) -> list[list[stim.GateTarget]]:
        # Initialize empty records list
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

    def _get_untouched_merging_records(self, patch_type: str) -> list[list[stim.GateTarget]]:
        # Initialize empty records list
        records: list[list[stim.GateTarget]] = []

        # Determine the previous patch type based on the current untouched merging patch type
        if patch_type == "Merge_AC_untouched":
            previous_patch_type = "Control_&_Target"
        elif patch_type == "Merge_AT_untouched":
            previous_patch_type = "Control_&_Target_Split_AC"

        for qubit, curr_abs_idx in self.detector_dict[patch_type]:
            # Calculate the Relative Index for the current measurement (for stim)
            current_rel_idx = curr_abs_idx - self.total_measurements

            # First check if this qubit was previously measured in this merge patch type
            full_history_entries = [
                entry for entry in self.full_measurement_history[patch_type] if entry[0] == qubit
            ]

            # If there are 2 entries in full history, treat like regular detector
            if len(full_history_entries) == 2:
                previous_rel_idx = full_history_entries[0][1] - self.total_measurements
                records.append(
                    [stim.target_rec(current_rel_idx), stim.target_rec(previous_rel_idx)],
                )

            # If only 1 entry in full history, this is the first time - check special cases
            elif len(full_history_entries) == 1:
                for qubit_ct, curr_abs_idx_ct in self.detector_dict[previous_patch_type]:
                    # Filtering For Regular Stabilizers: These are the ones which have a
                    # measurement entry in the Control & Target
                    if qubit_ct == qubit:
                        previous_rel_idx_ct = curr_abs_idx_ct - self.total_measurements
                        records.append(
                            [
                                stim.target_rec(current_rel_idx),
                                stim.target_rec(previous_rel_idx_ct),
                            ],
                        )

                        # Found the latest measurement entry
                        break

        return records

    def _get_merging_records(self, patch_type: str) -> list[list[stim.GateTarget]]:
        # Initialize empty records list
        records: list[list[stim.GateTarget]] = []

        # Determine the previous patch type based on the current merging patch type
        if patch_type == "Merge_AC":
            previous_patch_type_ct = "Control_&_Target"
            previous_patch_type_anc = "Ancilla"
            married_stab_indices = self.geometry.anc_x_bdy_b_stb_idx

        elif patch_type == "Merge_AT":
            previous_patch_type_ct = "Control_&_Target_Split_AC"
            previous_patch_type_anc = "Ancilla_Split_AC"
            married_stab_indices = self.geometry.anc_z_bdy_r_stb_idx

        for qubit, curr_abs_idx in self.detector_dict[patch_type]:
            # Initialize boolean to track if a partner record has been found
            partner_found = False

            # Calculate the Relative Index for the current measurement (for stim)
            current_rel_idx = curr_abs_idx - self.total_measurements

            # First check if this qubit was previously measured in this merge patch type
            full_history_entries = [
                entry for entry in self.full_measurement_history[patch_type] if entry[0] == qubit
            ]

            # If there are 2 entries in full history, treat like regular detector
            if len(full_history_entries) == 2:
                previous_rel_idx = full_history_entries[0][1] - self.total_measurements
                records.append(
                    [stim.target_rec(current_rel_idx), stim.target_rec(previous_rel_idx)],
                )

            # If only 1 entry in full history, this is the first time - check special cases
            elif len(full_history_entries) == 1:
                # Init Married Stab list for current qubit
                stored_rel_married_idx: int = 0

                for qubit_anc, curr_abs_idx_anc in self.detector_dict[previous_patch_type_anc]:
                    # Filtering For Regular Stabilizers: These are the ones which have a
                    # measurement entry in either the Ancilla Dict or the Control &
                    # Target Dict but not in both on the same qubit index
                    if qubit_anc == qubit and qubit_anc not in married_stab_indices:
                        previous_rel_idx_anc = curr_abs_idx_anc - self.total_measurements

                        records.append(
                            [
                                stim.target_rec(current_rel_idx),
                                stim.target_rec(previous_rel_idx_anc),
                            ],
                        )

                        # Found the latest measurement entry
                        partner_found = True
                        break

                    # Filtering For Married Stabilizers: These are the ones which have a
                    # measurement entry in both the Ancilla Dict and the Control & Target
                    if qubit_anc == qubit and qubit_anc in married_stab_indices:
                        previous_rel_idx_anc = curr_abs_idx_anc - self.total_measurements

                        # Store the Ancilla record and wait for the
                        # Control & Target record to appear in the loop below
                        stored_rel_married_idx = previous_rel_idx_anc

                        break

                if not partner_found:
                    for qubit_ct, curr_abs_idx_ct in self.detector_dict[previous_patch_type_ct]:
                        if qubit_ct == qubit and qubit_ct not in married_stab_indices:
                            previous_rel_idx_ct = curr_abs_idx_ct - self.total_measurements

                            records.append(
                                [
                                    stim.target_rec(current_rel_idx),
                                    stim.target_rec(previous_rel_idx_ct),
                                ],
                            )
                            # Found the latest measurement entry
                            break

                        # Filtering For Married Stabilizers: These are the ones which have a
                        # measurement entry in both the Ancilla Dict and the Control & Target
                        if qubit_ct == qubit and qubit_ct in married_stab_indices:
                            previous_rel_idx_ct = curr_abs_idx_ct - self.total_measurements

                            # Get Info of the corresponding Ancilla record stored in
                            records.append(
                                [
                                    stim.target_rec(current_rel_idx),
                                    stim.target_rec(previous_rel_idx_ct),
                                    stim.target_rec(stored_rel_married_idx),
                                ],
                            )
                            # Found the latest measurement entry
                            break

            else:
                raise ValueError(
                    f"Unexpected number of entries for qubit {qubit} in detector dict.",
                )

        return records

    def _get_splitting_records(self, patch_type: str) -> list[list[stim.GateTarget]]:
        # Initialize empty records list and boolean
        records: list[list[stim.GateTarget]] = []

        # Define seperated stab indices
        all_seperated_stab_indices = (
            self.geometry.anc_x_bdy_b_stb_idx + self.geometry.anc_z_bdy_r_stb_idx
        )

        # Determine the previous patch type based on the current splitting patch type
        if patch_type in {"Ancilla_Split_AC", "Control_&_Target_Split_AC"}:
            previous_patch_type_untouched = "Merge_AC_untouched"
            previous_patch_type_merged = "Merge_AC"
            specific_seperated_stab_indices = self.geometry.anc_x_bdy_b_stb_idx

        elif patch_type in {"Ancilla_Split_AT", "Control_&_Target_Split_AT"}:
            previous_patch_type_untouched = "Merge_AT_untouched"
            previous_patch_type_merged = "Merge_AT"
            specific_seperated_stab_indices = self.geometry.anc_z_bdy_r_stb_idx

        for qubit, curr_abs_idx in self.detector_dict[patch_type]:
            # Initialize boolean to track if a partner record has been found
            partner_found = False

            # Calculate the Relative Index for the current measurement (for stim)
            current_rel_idx = curr_abs_idx - self.total_measurements

            # Find all entries in the full history dict for the same qubit index
            full_history_entries = [
                entry for entry in self.full_measurement_history[patch_type] if entry[0] == qubit
            ]

            # If there are two entries, add both to the records list
            if len(full_history_entries) == 2:
                # Calculate the Relative Index for the previous measurement (for stim)
                previous_rel_idx = full_history_entries[0][1] - self.total_measurements
                records.append(
                    [stim.target_rec(current_rel_idx), stim.target_rec(previous_rel_idx)],
                )

            # If only one record exists, find older measurements inside Merge Dict
            # and add to records list
            elif len(full_history_entries) == 1:
                # Look for regular stabilizers in both Merge Dicts (merged and untouched)
                # and add to the combined dict
                for qubit_merge, abs_idx_merge in self.detector_dict[previous_patch_type_untouched]:
                    if qubit_merge == qubit and qubit not in all_seperated_stab_indices:
                        previous_rel_idx_merge = abs_idx_merge - self.total_measurements
                        records.append(
                            [
                                stim.target_rec(current_rel_idx),
                                stim.target_rec(previous_rel_idx_merge),
                            ],
                        )
                        partner_found = True
                        break

                # Continue looking in other Patch
                if not partner_found:
                    for qubit_merge, abs_idx_merge in self.detector_dict[
                        previous_patch_type_merged
                    ]:
                        if qubit_merge == qubit and qubit not in all_seperated_stab_indices:
                            previous_rel_idx_merge = abs_idx_merge - self.total_measurements
                            records.append(
                                [
                                    stim.target_rec(current_rel_idx),
                                    stim.target_rec(previous_rel_idx_merge),
                                ],
                            )
                            break

                # Look for divorced stabilizers in the Merge Dict and add to the combined dict
                # -> Records for divorced stabilizers are old record and both new records from
                #    Ancilla and Control & Target Split
                for qubit_merge, abs_idx_merge in self.detector_dict[previous_patch_type_merged]:
                    if qubit_merge == qubit and qubit in specific_seperated_stab_indices:
                        stored_list = self.storage_for_split_matching.get(qubit)

                        if stored_list is None:
                            self.storage_for_split_matching[qubit] = [
                                curr_abs_idx,
                                abs_idx_merge,
                                None,
                            ]
                            break
                        elif stored_list is not None:
                            records.append(
                                [
                                    stim.target_rec(current_rel_idx),
                                    stim.target_rec(stored_list[1] - self.total_measurements),
                                    stim.target_rec(stored_list[0] - self.total_measurements),
                                ],
                            )
                            # After finding both records for the divorced stabilizer, we can
                            # remove the entry from the storage dict
                            self.storage_for_split_matching.pop(qubit)
                            break
                        else:
                            raise ValueError(
                                f"Unexpected state of split matching storage for qubit {qubit}.",
                            )

            else:
                raise ValueError(
                    f"Unexpected number of entries for qubit {qubit} in full measurement history.",
                )

        return records

    def _get_final_measurement_records(self, patch_type: str) -> list[list[stim.GateTarget]]:
        """
        Helper Function to create Target Recod Pairs for Detectors out of the final Measuerement
        created by data measruements and the last measurement of real ancillary qubits

        Currently only used for the final measurement of the ancilla patch,
        but can be easily adapted
        """

        # Initialize empty records list
        records: list[list[stim.GateTarget]] = []

        # Create dict stab_idx to data_idx
        if patch_type == "Final_Measurement_Anc":
            measurement_basis = "Z"
            stab_idx_to_data_idx = self.geometry.get_patch_stabilizer_to_data_mapping(
                patch_coords=self.geometry.coords_ancilla,
                type=measurement_basis,
            )
            patch_before = "Ancilla_Split_AT"
        elif patch_type == "Final_Measurement_C":
            measurement_basis = None
            stab_idx_to_data_idx = self.geometry.get_patch_stabilizer_to_data_mapping(
                patch_coords=self.geometry.coords_ancilla,
                type=measurement_basis,
            )
            patch_before = "Control_&_Target_Split_AT"
        elif patch_type == "Final_Measurement_T":
            measurement_basis = None
            stab_idx_to_data_idx = self.geometry.get_patch_stabilizer_to_data_mapping(
                patch_coords=self.geometry.coords_ancilla,
                type=measurement_basis,
            )
            patch_before = "Control_&_Target_Split_AT"

        # Last Measurement -> Therefore only one Entry
        # Loop over Stab Measurements from past
        for qubit, curr_abs_idx in self.detector_dict[patch_before]:
            # Get the corresponding data qubit index for the current stabilizer qubit index
            data_qubit_idx = stab_idx_to_data_idx.get(qubit)

            # Skip different Basis Stabs
            if data_qubit_idx is None:
                continue

            # Create List for storage of all abs idx of data meas
            stored_list_of_abs_idx = [curr_abs_idx]

            # Loop over current data Measurements in present time step
            for qubit_data, curr_abs_idx_data in self.detector_dict[patch_type]:
                if qubit_data in data_qubit_idx:
                    # Store Abs Idx
                    stored_list_of_abs_idx.append(curr_abs_idx_data)

            # If loop is finished add get relative idx and add to records list
            curr_record_list = [
                stim.target_rec(curr_abs_idx - self.total_measurements)
                for curr_abs_idx in stored_list_of_abs_idx
            ]

            records.append(curr_record_list)

        return records
