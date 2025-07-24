from scipy.optimize import minimize_scalar
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import numpy as np
import sinter

__all__ = ["calc_threshold"]

#-----------------
# Global Function:
#-----------------

def threshold_approx(data_stats: list[sinter.TaskStats], p_min : float = 0.0025, p_max : float = 0.015) -> float:

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
        filtered_stats = [stat for stat in data_stats
                        if stat.json_metadata["distance"] == dist
                        and p_min < stat.json_metadata["p"] < p_max
]
        physical_p = [stat.json_metadata["p"] for stat in filtered_stats]
        logical_p = [stat.errors / stat.shots for stat in filtered_stats]

        """
        It is enough to check one Dataset for minim distance
        """
        if len(physical_p) < 2:
            raise ValueError("Dataset too small")

        ########################################################
        # Sort and convert and convert to log log for linear fit
        ########################################################

        filtered = [(p, l) for p, l in zip(physical_p, logical_p) if l > 0 and p > 0]
        if len(filtered) < 2:
            raise ValueError("Not enough valid points (logical_p > 0) for interpolation")

        physical_p, logical_p = zip(*sorted(filtered))

        x_dots[dist] = np.log10(physical_p)
        y_dots[dist] = np.log10(logical_p)
    
    ##########################################
    # Intepolate Data for root_scalar function
    ##########################################

    f1_interp = UnivariateSpline(x_dots[d1], y_dots[d1], s=1e-4)
    f2_interp = UnivariateSpline(x_dots[d2], y_dots[d2], s=1e-4)

    ##################################################################
    # Define function for root_scalar and boundaries for search region
    ##################################################################

    def diff(x):
        return (f1_interp(x) - f2_interp(x))**2

    x_min = max(min(x_dots[d1]), min(x_dots[d2]))
    x_max = min(max(x_dots[d1]), max(x_dots[d2]))

    #######################################
    # Run root scalar and check convergence
    #######################################

    sol = minimize_scalar(diff, bounds=(x_min, x_max), method='bounded')

    if not sol.success:
        raise RuntimeError("Minimization did not converge")

    #################################
    # Convert back from log10(x) to x
    #################################

    threshold = 10 ** sol.x
    return threshold