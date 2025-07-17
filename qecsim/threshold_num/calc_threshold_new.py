from __future__ import division 
import numpy as np
import matplotlib.pyplot as plt

__all__ = ["calc_threshold_new"]

def _interpolated_intercept(x, y1, y2):
    """Find the intercept of two curves, given by the same x data"""

    def _intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the second line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """    

        def _line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def _intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = _line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = _line([point3[0],point3[1]], [point4[0],point4[1]])

        R = _intersection(L1, L2)

        return R

    idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xc, yc = _intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
    return xc,yc

def calc_threshold_new(data_stats : list):

    """
    Calculates the crossing and therefore the threshold of Sinter.Stats with two distances

    Variables:
    -> data_stats : list (Sinter.Stats)

    Returns:
    -> float : Point of cossing i.e. threhshold level in percentage
    
    Info:
    -> If more than 2 distances are given, chooses the highest distances for the most accurate approximation
    """

    distances = sorted(set(stat.json_metadata["distance"] for stat in data_stats))

    ################################################
    # Settings up Dict and choosing highest distance
    ################################################

    d1 = distances[-1]
    d2 = distances[-2]
    x_dots : dict = {}
    y_dots : dict = {}

    ###########################################
    # Filter out Datasets to d1 and d2 distance
    ###########################################

    for dist in [d1, d2]:
        physical_p = [
            stat.json_metadata["p"]
            for stat in data_stats
            if stat.json_metadata["distance"] == dist
        ]

        logical_p = [
            stat.errors / stat.shots
            for stat in data_stats
            if stat.json_metadata["distance"] == dist
        ]

        """
        It is enough to check one Dataset for minimum distance
        """

        if len(physical_p) < 2:
            raise ValueError("Dataset too small")

        ########################################################
        # Sort and convert and convert to log log for linear fit
        ########################################################

        filtered = [(p, l) for p, l in zip(physical_p, logical_p) if l > 0 and p > 0]
        
        if len(filtered) < 2:
            raise ValueError("Not enough valid points i.e. p > 0")

        physical_p, logical_p = zip(*sorted(filtered))

        x_dots[dist] = np.log10(physical_p)
        y_dots[dist] = np.log10(logical_p)

    # d1 or d2 is irrelevant as they share the same x values
    x  = x_dots[d1]
    y1 = y_dots[d1]
    y2 = y_dots[d2]

    # Nearest corssing point with interpolation
    xc, yc = _interpolated_intercept(x,y1,y2)

    min_yc = 1e8
    min_xc = 0

    # Search lowest "crossing point"
    for x_val, y_val in zip(xc, yc):
        if y_val < min_yc:
            min_yc = y_val
            min_xc = x_val

    # Return to threshold value i.e. convert
    threshold_val = 10 ** min_xc

    return threshold_val[0]

