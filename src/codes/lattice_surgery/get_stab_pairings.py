from src.core.base_get_stab_pairings import BasePairings

Coord = complex
Label = str
Pair = tuple[Coord, Coord]

__all__ = ["LatticeSurgeryPairings"]


class StandardLatticePairings(BasePairings):
    def __init__(
        self,
        qubit_coords: dict[complex, str],
        merging: bool,
        merging_type: str | None = None,
    ):
        self.merging = merging
        self.merging_type = merging_type
        self.patch = qubit_coords
        self.stab_to_data: dict[Pair, str] = {}

    def get_schedule_standard(self) -> dict[Pair, str]:
        """
        Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

        * X-Stabs: control ist the **Data** qubit -> (data, stab)
        * Z-Stabs: control is the **stab** qubit -> (stab, data)
        """

        # Running the Functions to populate the Schedule
        self._attach_interior_standard()
        self._attach_boundary_standard()

        return self.stab_to_data

    def _attach_interior_standard(self):
        """
        Adds the 4-body CX Schedule for the *interior* stabilizers
        """

        for coords, string in self.patch.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if string == "X-STAB":
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX", "3-CX", "4-CX"],
                )
            elif string == "Z-STAB":
                pairs = [(coords, q1), (coords, q3), (coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX", "3-CX", "4-CX"],
                )

    def _attach_boundary_standard(self):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        for coords, string in self.patch.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if string == "Z-STAB-BOUND-L-A":
                pairs = [(coords, q1), (coords, q3)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX"],
                )

            elif string == "Z-STAB-BOUND-R-A":
                pairs = [(coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["3-CX", "4-CX"],
                )

            elif string == "X-STAB-BOUND-A-A":
                pairs = [(q4, coords), (q3, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "3-CX"],
                )

            elif string == "X-STAB-BOUND-B-A":
                pairs = [(q2, coords), (q1, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["2-CX", "1-CX"],
                )

            elif string == "Z-STAB-BOUND-L-T":
                pairs = [(coords, q1), (coords, q3)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["5-CX", "6-CX"],
                )

            elif string == "Z-STAB-BOUND-R-T":
                pairs = [(coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["5-CX", "6-CX"],
                )

            elif string == "X-STAB-BOUND-A-T":
                pairs = [(q4, coords), (q3, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["5-CX", "6-CX"],
                )

            elif string == "X-STAB-BOUND-B-T":
                pairs = [(q2, coords), (q1, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["5-CX", "6-CX"],
                )

            elif string == "Z-STAB-BOUND-L-C":
                pairs = [(coords, q1), (coords, q3)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["5-CX", "6-CX"],
                )

            elif string == "Z-STAB-BOUND-R-C":
                pairs = [(coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["5-CX", "6-CX"],
                )

            elif string == "X-STAB-BOUND-A-C":
                pairs = [(q4, coords), (q3, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["5-CX", "6-CX"],
                )

            elif string == "X-STAB-BOUND-B-C":
                pairs = [(q2, coords), (q1, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["5-CX", "6-CX"],
                )


class ACMergingLatticePairings(BasePairings):
    def __init__(
        self,
        qubit_coords: dict[complex, str],
        merging: bool,
        merging_type: str | None = None,
    ):
        self.merging = merging
        self.merging_type = merging_type
        self.patch = qubit_coords
        self.stab_to_data: dict[Pair, str] = {}

    def get_schedule_ac(self) -> dict[Pair, str]:
        """
        Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

        * X-Stabs: control ist the **Data** qubit -> (data, stab)
        * Z-Stabs: control is the **stab** qubit -> (stab, data)
        """

        # Running the Functions to populate the Schedule
        self._attach_interior_ac()
        self._attach_boundary_ac()

        return self.stab_to_data

    def _attach_interior_ac(self):
        """
        Adds the 4-body CX Schedule for the *interior* stabilizers
        """

        for coords, string in self.patch.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if string in {"X-STAB", "X-STAB-BOUND-A-C"}:
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX", "3-CX", "4-CX"],
                )

            elif string in {"Z-STAB", "Z-STAB-SURGERY-M"}:
                pairs = [(coords, q1), (coords, q3), (coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX", "3-CX", "4-CX"],
                )

    def _attach_boundary_ac(self):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        for coords, string in self.patch.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if string in {"Z-STAB-BOUND-L-A", "Z-STAB-BOUND-L-C", "Z-STAB-SURGERY-L"}:
                pairs = [(coords, q1), (coords, q3)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX"],
                )

            elif string in {"Z-STAB-BOUND-R-A", "Z-STAB-BOUND-R-C"}:
                pairs = [(coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["3-CX", "4-CX"],
                )

            elif string == "X-STAB-BOUND-A-A":
                pairs = [(q4, coords), (q3, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "3-CX"],
                )

            elif string == "X-STAB-BOUND-B-C":
                pairs = [(q2, coords), (q1, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["2-CX", "1-CX"],
                )


class ATMergingLatticePairings(BasePairings):
    def __init__(
        self,
        qubit_coords: dict[complex, str],
        merging: bool,
        merging_type: str | None = None,
    ):
        self.merging = merging
        self.merging_type = merging_type
        self.patch = qubit_coords
        self.stab_to_data: dict[Pair, str] = {}

    def get_schedule_at(self) -> dict[Pair, str]:
        """
        Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

        * X-Stabs: control ist the **Data** qubit -> (data, stab)
        * Z-Stabs: control is the **stab** qubit -> (stab, data)
        """

        # Running the Functions to populate the Schedule
        self._attach_interior_at()
        self._attach_boundary_at()

        return self.stab_to_data

    def _attach_interior_at(self):
        """
        Adds the 4-body CX Schedule for the *interior* stabilizers
        """

        for coords, string in self.patch.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if string in {"X-STAB", "X-STAB-SURGERY-M"}:
                pairs = [(q1, coords), (q2, coords), (q3, coords), (q4, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX", "3-CX", "4-CX"],
                )

            elif string in {"Z-STAB", "Z-STAB-BOUND-L-T"}:
                pairs = [(coords, q1), (coords, q3), (coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX", "3-CX", "4-CX"],
                )

    def _attach_boundary_at(self):
        """
        Adds the 2-body CX Schedule for the *boundary* stabilizers
        """

        for coords, string in self.patch.items():
            # Get neigbouring data coords
            q1 = self._neighbours(coords, +1, -1)
            q2 = self._neighbours(coords, -1, -1)
            q3 = self._neighbours(coords, +1, +1)
            q4 = self._neighbours(coords, -1, +1)

            if string == "Z-STAB-BOUND-L-A":
                pairs = [(coords, q1), (coords, q3)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["1-CX", "2-CX"],
                )

            elif string == "Z-STAB-BOUND-R-T":
                pairs = [(coords, q2), (coords, q4)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["3-CX", "4-CX"],
                )

            elif string in {"X-STAB-BOUND-A-A", "X-STAB-BOUND-A-T"}:
                pairs = [(q4, coords), (q3, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["4-CX", "3-CX"],
                )

            elif string in {"X-STAB-BOUND-B-A", "X-STAB-BOUND-B-T", "X-STAB-SURGERY-B"}:
                pairs = [(q2, coords), (q1, coords)]
                self._assign_orders(
                    self.stab_to_data,
                    pairs,
                    ["2-CX", "1-CX"],
                )


class LatticeSurgeryPairings(
    StandardLatticePairings,
    ACMergingLatticePairings,
    ATMergingLatticePairings,
):
    def __init__(
        self,
        qubit_coords: dict[complex, str],
        merging: bool,
        merging_type: str | None = None,
    ):
        if merging_type not in {"AC", "AT", None}:
            raise ValueError("merging_type must be 'AC', 'AT' or None")

        self.merging = merging
        self.merging_type = merging_type
        self.patch = qubit_coords

    def get_schedule(self) -> dict[Pair, str]:
        """
        Returns the CX-Schedule {(data_coord, stab_coord): order} for a given lattice

        * X-Stabs: control ist the **Data** qubit -> (data, stab)
        * Z-Stabs: control is the **stab** qubit -> (stab, data)
        """

        # Initlizing Dict
        self.stab_to_data: dict[Pair, str] = {}

        # Running the Functions to populate the Schedule
        if self.merging and self.merging_type == "AC":
            self.get_schedule_ac()
        elif self.merging and self.merging_type == "AT":
            self.get_schedule_at()
        else:
            self.get_schedule_standard()

        return self.stab_to_data
