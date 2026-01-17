import stim

from src.codes.lattice_surgery.data_geometry import MasterPairings
from src.codes.lattice_surgery.surgery_geom import SurgeryGeometry
from src.core.cx_builder import cx_builder

Coord = complex

__all__ = ["SurgerySplit"]


class SurgerySplit:
    def __init__(
        self,
        geometry: SurgeryGeometry,
        master_pairings: MasterPairings,
        split_type: str,
    ):
        # Preliminary Setup
        self.geometry = geometry
        self.master_pairings = master_pairings

        # Getting Specific Splitting Type Info
        self.split_type = split_type
        if split_type == "AT":
            self.x_stab_index_untouched_circ = self.geometry.control_x_stb_idx
            self.z_stab_index_untouched_circ = self.geometry.control_z_stb_idx
            self.combined_x_stab_merging_lattices, self.combined_z_stab_merging_lattices = (
                self.geometry.get_combined_xz_stabs_merging_lattice(split_type="AT")
            )
        elif split_type == "AC":
            self.x_stab_index_untouched_circ = self.geometry.target_x_stb_idx
            self.z_stab_index_untouched_circ = self.geometry.target_z_stb_idx
            self.combined_x_stab_merging_lattices, self.combined_z_stab_merging_lattices = (
                self.geometry.get_combined_xz_stabs_merging_lattice(split_type="AC")
            )
        else:
            raise ValueError("No valid splitting Type in Function selected!")

    def build_circuit(self) -> stim.Circuit:
        # Init return Circuit
        return_circuit = stim.Circuit()

        # Add Inital Split Round
        return_circuit += self._init_split_circuit()

        # Add Repetition Rounds
        return_circuit += self._repeat_split_circuit()

        # Add Final Split Round
        return_circuit += self._final_split_circuit()

        return return_circuit

    def _init_split_circuit(self):
        # Define initial split Circuit
        split_init_circuit = stim.Circuit()

        # Adding resets from merge ac & at
        split_init_circuit.append("TICK")
        split_init_circuit.append(
            "R",
            self.x_stab_index_untouched_circ + self.z_stab_index_untouched_circ,
        )

        split_init_circuit.append("TICK")
        split_init_circuit.append("H", self.geometry.combined_x_stab_idx_filtered)

        # CX Operations for Ancilla
        split_init_circuit.append("TICK")
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_init_circuit,
        )

        # Retreive Boundary + Normal Stabilizers Ancilla
        # I.e. Basis switch and measurement of ancillas
        split_init_circuit.append("H", self.geometry.anc_x_stb_idx)
        split_init_circuit.append("TICK")

        split_init_circuit.append("M", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_init_circuit.append("TICK")

        split_init_circuit.append("R", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_init_circuit.append("TICK")

        split_init_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
        split_init_circuit.append("TICK")

        # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_init_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Retreive Boundary + Normal Stabilizers from Target and Control
        # (Basis Change + Measurement):
        split_init_circuit.append(
            "H",
            self.geometry.control_x_stb_idx + self.geometry.target_x_stb_idx,
        )
        split_init_circuit.append("TICK")

        split_init_circuit.append("M", self.geometry.control_target_all_stab_idx)

        return split_init_circuit

    def _repeat_split_circuit(self):
        # Implementing Repeat Block
        split_repeat_circuit = stim.Circuit()

        # Reset Ancilla and prepare measurement basis
        split_repeat_circuit.append("TICK")
        split_repeat_circuit.append("R", self.geometry.control_target_all_stab_idx)

        split_repeat_circuit.append("TICK")
        split_repeat_circuit.append("H", self.geometry.combined_x_stab_idx_filtered)

        split_repeat_circuit.append("TICK")

        # CX Operations for Ancilla
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_repeat_circuit,
        )

        # Retreive Boundary + Normal Stabilizers Ancilla:
        # I.e. Basis switch and measurement of ancillas
        split_repeat_circuit.append("H", self.geometry.anc_x_stb_idx)
        split_repeat_circuit.append("TICK")

        split_repeat_circuit.append("M", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_repeat_circuit.append("TICK")

        split_repeat_circuit.append("R", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_repeat_circuit.append("TICK")

        split_repeat_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
        split_repeat_circuit.append("TICK")

        # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_repeat_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Retreive Boundary + Normal Stabilizers from Target and Control
        # (Basis Change + Measurement):
        split_repeat_circuit.append(
            "H",
            self.geometry.control_x_stb_idx + self.geometry.target_x_stb_idx,
        )
        split_repeat_circuit.append("TICK")

        split_repeat_circuit.append("M", self.geometry.control_target_all_stab_idx)

        return split_repeat_circuit * (self.geometry.distance - 2)

    def _final_split_circuit(self):
        # Adding Final Circ
        """
        In this Section we add the Conditional X_L and Z_L depending on the XX and ZZ Measurements
        """

        split_final_circuit = stim.Circuit()

        # Adding Reset and basis preparation
        split_final_circuit.append("TICK")
        split_final_circuit.append("R", self.geometry.control_target_all_stab_idx)

        split_final_circuit.append("TICK")
        split_final_circuit.append("H", self.geometry.combined_x_stab_idx_filtered)

        split_final_circuit.append("TICK")

        # CX Operations for Ancilla
        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_final_circuit,
        )

        # Retreive Boundary + Normal Stabilizers Ancilla:
        # I.e. Basis switch and measurement of ancillas
        split_final_circuit.append("H", self.geometry.anc_x_stb_idx)
        split_final_circuit.append("TICK")

        split_final_circuit.append("M", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_final_circuit.append("TICK")

        split_final_circuit.append("R", self.geometry.anc_x_stb_idx + self.geometry.anc_z_stb_idx)
        split_final_circuit.append("TICK")

        split_final_circuit.append("H", self.geometry.anc_x_bdy_b_stb_idx)
        split_final_circuit.append("TICK")

        # -------------------------------------------------------------------------------
        """
        As these are dependent on the measurement outcome and also on the type of merge
        -> Seperated into two Parts conditioned on the Merge/Split type

        The logical ZZ Observable is already defined by the newly implemented Z stabilizers 
        on the merge
        -> Product of the stabilizers give measurement result
        """

        # Determining Position in the measurement Run of only the Ancilla
        combined_x_stab_merging_lattices: list = []
        combined_z_stab_merging_lattices: list = []

        """
        Only Comparing the inner stabilizers and the boundarys who are 
        not involed in the merging with one another
        """

        if self.split_type == "AT":
            # Adding h gate for X stabilizers only on merging lattices
            # -> Filtering out double coords in big lattice
            combined_x_stab_merging_lattices, combined_z_stab_merging_lattices = (
                self.geometry.get_combined_xz_stabs_merging_lattice(split_type="AT")
            )

        # Finding newly generated Stabilizer postion in the Measurement-Rec of the Merge
        pos_to_index_newly_gen_stabs: list = []

        if self.split_type == "AC":
            for pos, index in enumerate(
                combined_z_stab_merging_lattices + combined_x_stab_merging_lattices,
            ):
                if index in (self.geometry.surgery_z_m_stb_idx + self.geometry.surgery_z_l_stb_idx):
                    pos_to_index_newly_gen_stabs.append([pos, index])

        elif self.split_type == "AT":
            for pos, index in enumerate(
                combined_z_stab_merging_lattices + combined_x_stab_merging_lattices,
            ):
                if index in (self.geometry.surgery_x_m_stb_idx + self.geometry.surgery_x_b_stb_idx):
                    pos_to_index_newly_gen_stabs.append([pos, index])

        logical_obs_rec_tar = []

        inner_record = len(combined_z_stab_merging_lattices + combined_x_stab_merging_lattices)
        skipped_records = len(
            self.x_stab_index_untouched_circ
            + self.z_stab_index_untouched_circ
            + self.geometry.distance
            * (
                self.geometry.control_target_all_stab_idx
                + self.geometry.anc_x_stb_idx
                + self.geometry.anc_z_stb_idx,
            ),
        )

        for index_pos_merge in pos_to_index_newly_gen_stabs:
            logical_obs_rec_tar.append(index_pos_merge[0] - inner_record - skipped_records)

        #####################################################
        # Defining Logical Data Qubit string for all lattices
        #####################################################

        log_strings = self.geometry._get_logical_strings()
        a_log_obs_z_index = log_strings["a_z"]
        t_log_obs_z_index = log_strings["t_z"]
        c_log_obs_x_index = log_strings["c_x"]

        ####################################
        # Adding Conditional CZ/CX-Operators
        ####################################

        if self.split_type == "AC":
            for records in logical_obs_rec_tar:
                for data in c_log_obs_x_index:
                    split_final_circuit.append("CX", [stim.target_rec(records), data])

        elif self.split_type == "AT":
            for records in logical_obs_rec_tar:
                for data in t_log_obs_z_index:
                    split_final_circuit.append("CZ", [stim.target_rec(records), data])

            ###################################
            # Meassuring Ancilla in the Z Basis
            ###################################

            # Meassuring Data
            split_final_circuit.append("TICK")
            split_final_circuit.append("MZ", self.geometry.anc_data_idx)
            split_final_circuit.append("TICK")

            # Adding the conditional Gate on Control
            for rec_tar, index in enumerate(self.geometry.anc_data_idx):
                if index in a_log_obs_z_index:
                    for data in c_log_obs_x_index:
                        split_final_circuit.append(
                            "CX",
                            [stim.target_rec(-len(self.geometry.anc_data_idx) + rec_tar), data],
                        )

        # Continue CX-Implementation for Target and Control (As Ancilla already has a full run)
        split_final_circuit.append("TICK")

        cx_builder(
            q2i=self.geometry.q2i,
            stab_to_data=self.master_pairings.std_pairings.get_schedule(),
            circuit=split_final_circuit,
            orders=("5-CX", "6-CX"),
        )

        # Retreive Boundary + Normal Stabilizers from Target and Control
        # (Basis Change + Measurement):
        split_final_circuit.append(
            "H",
            self.geometry.control_x_stb_idx + self.geometry.target_x_stb_idx,
        )
        split_final_circuit.append("TICK")

        split_final_circuit.append("M", self.geometry.control_target_all_stab_idx)

        return split_final_circuit
