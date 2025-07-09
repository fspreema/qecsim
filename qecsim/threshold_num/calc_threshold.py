from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
import numpy as np
from collections import defaultdict
import sinter

__all__ = ["calc_threshold"]

#-----------------
# Global Function:
#-----------------

def threshold_linear(data_stats: list[sinter.TaskStats], p_min : float = 0.0025, p_max : float = 0.015) -> float:

    """
    Returns: float
    -> Calculated threshold

    Function selects stats with highest distance and calculates corssing point by interpolation
    -> Threshold is defined as the limit of corssings between d & d+2 distances
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
            if stat.json_metadata["distance"] == dist and p_min < stat.json_metadata["p"] < p_max
        ]
        logical_p = [
            stat.errors / stat.shots
            for stat in data_stats
            if stat.json_metadata["distance"] == dist and p_min < stat.json_metadata["p"] < p_max
        ]

        """
        It is enough to check one Dataset for minim distance
        """
        if len(physical_p) < 2:
            raise ValueError("Dataset too small")

        ########################################################
        # Sort and convert and convert to log log for linear fit
        ########################################################
        physical_p, logical_p = zip(*sorted(zip(physical_p, logical_p)))

        x_dots[dist] = np.log10(physical_p)
        y_dots[dist] = np.log10(logical_p)
    
    ##########################################
    # Intepolate Data for root_scalar function
    ##########################################

    f1_interp = interp1d(x_dots[d1], y_dots[d1], kind='linear', bounds_error=False)
    f2_interp = interp1d(x_dots[d2], y_dots[d2], kind='linear', bounds_error=False)

    ##################################################################
    # Define function for root_scalar and boundaries for search region
    ##################################################################

    def diff(x):
        return f1_interp(x) - f2_interp(x)

    x_min = max(min(x_dots[d1]), min(x_dots[d2]))
    x_max = min(max(x_dots[d1]), max(x_dots[d2]))

    #######################################
    # Run root scalar and check convergence
    #######################################

    sol = root_scalar(diff, bracket=[x_min, x_max], method='brentq')

    if not sol.converged:
        raise RuntimeError("Root finding did not converge")

    #################################
    # Convert back from log10(x) to x
    #################################

    threshold = 10 ** sol.root
    return threshold