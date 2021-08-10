import matplotlib.pyplot as plt
import numpy as np

pos = [0.5, 0.505, 0.5025, 0.51, 0.515]

forces = [
    -1.2148627259717075,
    -1.1839872523561565,
    -1.1992885235556070,
    -1.1541824759833621,
    -1.1254017883815202
]

plt.plot(pos, forces - np.min(forces), 'r*')
plt.show()
