from pprint import pprint
from numba import jit
import numpy as np


def vec(r_1, r_2, T):
    delta = r1 - r2 - T
    return delta / np.sqrt(np.sum(delta ** 2)) ** 3


Z1 = 25
Z2 = 13

r1 = np.array([0.0, 0.0, 0.0])
r2 = np.array([0.5, 0.5, 0.5])

force = np.zeros(3)

t = 40

for a in range(-t, t):
    for b in range(-t, t):
        for c in range(-t, t):
            # With the periodic image
            T = np.array([a, b, c], dtype=np.float64)

            res = vec(r1, r2, T)

            force += - Z1 * Z2 * res

pprint(force)

r1 = np.array([0.0, 0.0, 0.0])
r2 = np.array([0.5, 0.5, 0.6])

force = np.zeros(3)

for a in range(-t, t):
    for b in range(-t, t):
        for c in range(-t, t):
            # With the periodic image
            T = np.array([a, b, c], dtype=np.float64)

            res = vec(r1, r2, T)

            force += - Z1 * Z2 * res

pprint(force)
