class Mapping:
    # Setting up Class Items
    def __init__(self, surface_dict, color_dict):
        self.surface_dict = surface_dict
        self.color_dict = color_dict

        # Get Joint Coordinates
        self.joint_coords: dict = {}
        for key in self.surface_dict.keys():
            self.joint_coords[key] = self.surface_dict[key]

        for key_2 in self.color_dict.keys():
            self.joint_coords[key_2] = self.color_dict[key_2]

        # Get Index Mappings
        self.q2i: dict[complex, int] = {
            q: i
            for i, q in enumerate(
                sorted(self.joint_coords, key=lambda v: (v.real, v.imag)),
            )
        }

    # Defining Class Helper Function
    def _pick_up_indices(
        self,
        coords: dict[complex, str],
        q2i: dict[complex, int],
        *labels: str,
    ) -> list:
        """
        Return q2i indices matching any of the provided labels in coords.
        """
        label_set = set(labels)
        return [q2i[q] for q, t in coords.items() if t in label_set]

    def get_data_indexes(self, type: str):
        # check Validity of Type
        if type not in ["surface", "color"]:
            raise ValueError("Type must be either 'surface' or 'color'")

        # Return Data Indexes for Surface
        if type == "surface":
            return self._pick_up_indices(self.surface_dict, self.q2i, "DATA")

        # Return Data Indexes for Color
        elif type == "color":
            return self._pick_up_indices(self.color_dict, self.q2i, "DATA")

    def get_x_stabilizer_indexes(self, type: str):
        # check Validity of Type
        if type not in ["surface", "color"]:
            raise ValueError("Type must be either 'surface' or 'color'")

        # Return X Stabilizer Indexes for Surface
        if type == "surface":
            return self._pick_up_indices(self.surface_dict, self.q2i, "X")

        # Return X Stabilizer Indexes for Color
        elif type == "color":
            return self._pick_up_indices(self.color_dict, self.q2i, "X")

    def get_z_stabilizer_indexes(self, type: str):
        # check Validity of Type
        if type not in ["surface", "color"]:
            raise ValueError("Type must be either 'surface' or 'color'")

        # Return Z Stabilizer Indexes for Surface
        if type == "surface":
            return self._pick_up_indices(self.surface_dict, self.q2i, "Z")

        # Return Z Stabilizer Indexes for Color
        elif type == "color":
            return self._pick_up_indices(self.color_dict, self.q2i, "Z")
