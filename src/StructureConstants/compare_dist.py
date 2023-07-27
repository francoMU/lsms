
from pymatgen.core.lattice import Lattice
from pymatgen.util.coord import pbc_shortest_vectors

import numpy as np

lattice = Lattice(
    [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
)

f1 = [0, 0, 0]
f2 = [0.9, 0.9, 0.9]


v, d2 = pbc_shortest_vectors(lattice, f1, f2, return_d2=True)

print(v)
print(d2)
print(np.sqrt(d2))
print("-----")
