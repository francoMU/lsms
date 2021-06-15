from pyshtools.expand import spharm_lm

import matplotlib.pyplot as plt
import numpy as np


@np.vectorize
def spherical_harmonics(phi, theta, l, m):
    return spharm_lm(l, m, theta, phi,
                     normalization='ortho',
                     degrees=False,
                     csphase=-1)


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Create the mesh in polar coordinates and compute corresponding Z.
r = np.linspace(0, 1.25, 50)
phi = np.linspace(0, 2 * np.pi, 50)
R, P = np.meshgrid(r, phi)
theta = np.pi / 2

Z = np.zeros_like(R)

Z += spherical_harmonics(P, theta, 0, 0)
Z += spherical_harmonics(P, theta, 1, -1)
Z += spherical_harmonics(P, theta, 1, 0)
Z += spherical_harmonics(P, theta, 1, 1)

# Express the mesh in the cartesian system.
X, Y = R * np.cos(P), R * np.sin(P)

# Plot the surface.
ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)

# Tweak the limits and add latex math labels.
ax.set_zlim(0, 1)
ax.set_xlabel(r'$\phi_\mathrm{real}$')
ax.set_ylabel(r'$\phi_\mathrm{im}$')
ax.set_zlabel(r'$V(\phi)$')

plt.show()
