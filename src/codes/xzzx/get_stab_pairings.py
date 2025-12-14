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

    def get_schedule(self) -> dict[Pair, str]:
        """
        Populate schedule for XZZX-style stabilizers
        """

        stab_to_data: dict[Pair, str] = {}

        def _attach_interior_cx():
            """
            Adds the 4-body CX Schedule for the *interior* stabilizers
            """

            for coords, string in self.patch.items():
                # Get neigbouring data coords
                q1 = self._neighbours(coords, -1, -1)
                q2 = self._neighbours(coords, +1, -1)
                q3 = self._neighbours(coords, -1, +1)
                q4 = self._neighbours(coords, +1, +1)

                if string in {"STAB-Ver", "STAB-Hor"}:
                    pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                    self._assign_orders(
                        stab_to_data,
                        pairs,
                        ["1-CX", "2-CZ", "3-CZ", "4-CX"],
                    )

        def _attach_boundary_cx():
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

        # Add interior and boundary CXs
        _attach_interior_cx()
        _attach_boundary_cx()

        return stab_to_data
