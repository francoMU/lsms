from pprint import pprint

import bisect
import pandas as pd
import numpy as np
from scipy.integrate import trapezoid

if __name__ == "__main__":

    radius_sphere = 2.71929

    data = pd.read_csv("interpol_2.dat",
                       delim_whitespace=True,
                       header=None,
                       names=['R', 'rho'])

    pprint(data)

    R = data.loc[:, 'R'].to_numpy()
    rho = data.loc[:, 'rho'].to_numpy()

    element = bisect.bisect_left(R, radius_sphere)

    trapez_result = trapezoid(rho[:element], R[:element])

    print(trapez_result)