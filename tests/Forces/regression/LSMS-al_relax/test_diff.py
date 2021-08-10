from pprint import pprint

from numpy import genfromtxt

import numpy as np
import pandas as pd
import sys

comp_1 = sys.argv[1]
comp_2 = sys.argv[2]

df_1 = pd.read_csv(comp_1, sep='\s+', header=None)
array_1 = df_1.to_numpy()

df_2 = pd.read_csv(comp_2, sep='\s+', header=None)
array_2 = df_2.to_numpy()

for i, j in zip(array_1, array_2):
    np.testing.assert_allclose(i, j)
