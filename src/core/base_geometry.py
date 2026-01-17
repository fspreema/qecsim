from abc import ABC, abstractmethod

"""
Central Geometry Class
======================

Houses central functions used across different code geometries.
"""


class BaseGeometry(ABC):
    def __init__(self, **kwargs):
        self.coords = self.get_coords(**kwargs)
        self.q2i = self._get_q2i()
        self.i2q = self._get_i2q()

    @abstractmethod
    def get_coords(self, **kwargs) -> dict[complex, str]:
        """
        Abstract Method to get all coordinates with their labels.

        Returns:
            dict[complex, str]
                Dictionary with coordinates as keys and labels as values
        """
        pass

    @abstractmethod
    def _get_central_labels(self, **kwargs) -> dict[complex, str]:
        """
        Abstract Method to get only the qubit coordinates without boundaries.

        Returns:
            dict[complex, str]
                Dictionary with data qubit coordinates as keys and labels as values
        """
        pass

    @abstractmethod
    def _get_boundary_labels(self, **kwargs) -> dict[complex, str]:
        """
        Abstract Method to get only the boundary qubit coordinates.

        Returns:
            dict[complex, str]
                Dictionary with boundary qubit coordinates as keys and labels as values
        """
        pass

    def _get_q2i(self) -> dict[complex, int]:
        """
        Returns the index mapping of lattice points.
        """

        index = {}
        for idx, coord in enumerate(self.coords.keys()):
            index[coord] = int(idx)
        return index

    def _get_i2q(self) -> dict[int, complex]:
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

    def _get_specific_indices(
        self,
        type: str,
        valid_coords: dict[complex, str] = None,
    ) -> list[int]:
        """
        Returns only the qubit indices with label exactly matching the type.

        Args:
            type (str): Label type to filter for exact matches
            valid_coords (dict[complex, str], optional): If given, only consider these coordinates
                instead of all coordinates in self.coords.
        """

        if valid_coords is None:
            source = self.coords
        else:
            source = valid_coords

        return [self.q2i[coord] for coord, label in source.items() if label == type]
