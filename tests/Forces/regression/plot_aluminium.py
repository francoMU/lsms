import matplotlib.pyplot as plt
import numpy as np

pos = [0.5, 0.505, 0.51, 0.52]

forces = [0.0,
    -0.0392367777299251,
    -0.0786509694578151,
    -0.1578459851115376]

forces_mst = [
    -0.33982783e-1,
    -0.67813213e-1,
    -0.13573077,
]


forces_numeric = [0]*3

numeric_delta = [-0.001, 0.001]
numeric_component = [-1931.724839370260270, -1931.724767254148219]
delta = (numeric_component[1] - numeric_component[0])/(numeric_delta[1] - numeric_delta[0])
print(delta)
forces_numeric[0] = delta

numeric_delta = [-0.001, 0.001]
numeric_component = [-1931.721590436598035, -1931.721493352251855]
delta = (numeric_component[1] - numeric_component[0])/(numeric_delta[1] - numeric_delta[0])
print(delta)
forces_numeric[1] = delta

numeric_delta = [-0.001, 0.001]
numeric_component = [-1931.712114198335257, -1931.711941598344083]
delta = (numeric_component[1] - numeric_component[0])/(numeric_delta[1] - numeric_delta[0])
print(delta)
forces_numeric[2] = delta


plt.plot(pos, forces, 'r*', label="LSMS")
plt.plot(pos[1:4], -np.array(forces_numeric), 'b+', label="Numeric LSMS")
plt.plot(pos[1:4], np.array(forces_mst), 'cx', label="MST FP")
plt.legend()
plt.show()
