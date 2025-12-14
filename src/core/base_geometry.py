"""
Central Geometry Class
======================

Houses central functions used across different code geometries.
"""


class BaseGeometry:
    def _get_q2i(self):
        """
        Returns the index mapping of lattice points.
        """

        index = {}
        for idx, coord in enumerate(self.coords.keys()):
            index[coord] = int(idx)
        return index

    def _get_i2q(self):
        """
        Returns the coordinate mapping of indices.
        """

        index = {}
        for idx, coord in enumerate(self.coords.keys()):
            index[int(idx)] = coord

        return index

    def _get_specific_coords(self, type: str) -> dict[complex, str]:
        """
        Returns only the qubit coordinates with type corresponding to the label.

        Returns:
            dict[complex, str]
                Dictionary with data qubit coordinates as keys and labels as values
        """

        data_qubits = {coord: label for coord, label in self.coords.items() if type in label}
        return data_qubits
