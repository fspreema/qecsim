Coord = complex
Label = str
Pair = tuple[Coord, Coord]


# Create the abstract base class every pairing class will inherit from
class BasePairings:
    def __init__(self, patch: dict[complex, str], distance: int):
        self.patch = patch
        self.distance = distance

    def _neighbours(
        self,
        c: complex,
        dx: int,
        dy: int,
    ) -> complex:
        """
        Returns the neigbouring complex number with the real distance of dx
        and imag distance of dy
        """

        return (c.real + dx) + (c.imag + dy) * 1j

    def _assign_orders(
        self,
        table: dict[complex, str],
        pairs: list[complex],
        orders: list[str],
    ) -> None:
        """
        Adds the list of pairs into stab to data dict with the corresponding orders
        """

        for (a, b), order in zip(pairs, orders, strict=True):
            table[(a, b)] = order
