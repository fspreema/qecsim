from src.core.base_get_stab_pairings import BasePairings

Coord = complex
Label = str
Pair = tuple[Coord, Coord]

__all__ = ["XZZXPairings"]


class XZZXPairings(BasePairings):
    def __init__(
        self,
        patch: dict[complex, str],
    ):
        self.patch = patch
        self.stab_to_data: dict[Pair, str] = {}
        self.stab_to_data = self.get_schedule(self.stab_to_data)

    def get_schedule(self, stab_to_data: dict[Pair, str]) -> dict[Pair, str]:
        """
        Populate schedule for XZZX-style stabilizers
        """

        # Add interior and boundary CXs
        self._attach_interior_cx(stab_to_data)
        self._attach_boundary_cx(stab_to_data)

        return stab_to_data

    def _attach_interior_cx(self, stab_to_data: dict[Pair, str]):
        """
        Adds the 4-body CX Schedule for the *interior* stabilizers
        """

        for coords, string in self.patch.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, -1, -1)
            q2 = self._neighbours(coords, +1, -1)
            q3 = self._neighbours(coords, -1, +1)
            q4 = self._neighbours(coords, +1, +1)

            if string in {"STAB-VER", "STAB-HOR"}:
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                self._assign_orders(
                    stab_to_data,
                    pairs,
                    ["1-CX", "2-CZ", "3-CZ", "4-CX"],
                )

        return stab_to_data

    def _attach_boundary_cx(self, stab_to_data: dict[Pair, str]):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        for coords, string in self.patch.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, -1, -1)
            q2 = self._neighbours(coords, +1, -1)
            q3 = self._neighbours(coords, -1, +1)
            q4 = self._neighbours(coords, +1, +1)

            if string == "STAB-BOUND-L-Hor":
                pairs = [(q2, coords), (q4, coords)]
                self._assign_orders(
                    stab_to_data,
                    pairs,
                    ["3-CZ", "4-CX"],
                )
            elif string == "STAB-BOUND-R-Hor":
                pairs = [(q1, coords), (q3, coords)]
                self._assign_orders(
                    stab_to_data,
                    pairs,
                    ["1-CX", "2-CZ"],
                )
            elif string == "STAB-BOUND-A-Ver":
                pairs = [(q3, coords), (q4, coords)]
                self._assign_orders(
                    stab_to_data,
                    pairs,
                    ["3-CZ", "4-CX"],
                )
            elif string == "STAB-BOUND-B-Ver":
                pairs = [(q1, coords), (q2, coords)]
                self._assign_orders(
                    stab_to_data,
                    pairs,
                    ["1-CX", "2-CZ"],
                )

        return stab_to_data
