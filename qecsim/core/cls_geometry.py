import stim

"""
Geometry Class which build all Coordinates and converts them to Indices for the Circuits

-> Can later be added together by the __add__ function if no coordinates overlap
-> If overlap, offset can be adjusted
"""

class Geometry:
    def __init__(self, distance: int, type: str, offset: complex = 0 + 0j):
        self.distance = distance
        self.type = type
        self.offset = offset

    