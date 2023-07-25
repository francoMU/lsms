from numpy import dot, pi, transpose
from numpy.linalg import inv
import numpy as np
import sys

matrix = np.array(
    [[2.0, 0.0, 0.0],
     [0.1, 1.8, 0.0],
     [0.1, 0.2, 0.9]])

a = matrix.copy().T

delta = 0.75

b = np.zeros((3, 3))  # Vectors after the Gram-Schmidt process
u = np.zeros((3, 3))  # Gram-Schmidt coefficients
m = np.zeros(3)  # These are the norm squared of each vec.

b[:, 0] = a[:, 0]
m[0] = dot(b[:, 0], b[:, 0])
for i in range(1, 3):
    u[i, 0:i] = dot(a[:, i].T, b[:, 0:i]) / m[0:i]
    b[:, i] = a[:, i] - dot(b[:, 0:i], u[i, 0:i].T)
    m[i] = dot(b[:, i], b[:, i])

print(u)
print(b)
print(m)


k = 2

mapping = np.identity(3, dtype=np.double)
while k <= 3:

    print("Size reduction")

    # Size reduction.
    for i in range(k - 1, 0, -1):
        q = round(u[k - 1, i - 1])
        if q != 0:
            # Reduce the k-th basis vector.
            a[:, k - 1] = a[:, k - 1] - q * a[:, i - 1]
            mapping[:, k - 1] = mapping[:, k - 1] - q * mapping[:, i - 1]
            uu = list(u[i - 1, 0: (i - 1)])
            uu.append(1)
            # Update the GS coefficients.
            u[k - 1, 0:i] = u[k - 1, 0:i] - q * np.array(uu)


    print(a)
    print(mapping)
    print(u)

    # Check the Lovasz condition.
    if dot(b[:, k - 1], b[:, k - 1]) >= (delta - abs(u[k - 1, k - 2]) ** 2) * dot(b[:, (k - 2)], b[:, (k - 2)]):
        # Increment k if the Lovasz condition holds.
        print("Step 1")
        k += 1
    else:
        print("Step 2")
        # If the Lovasz condition fails,
        # swap the k-th and (k-1)-th basis vector
        v = a[:, k - 1].copy()
        a[:, k - 1] = a[:, k - 2].copy()
        a[:, k - 2] = v

        v_m = mapping[:, k - 1].copy()
        mapping[:, k - 1] = mapping[:, k - 2].copy()
        mapping[:, k - 2] = v_m



        # Update the Gram-Schmidt coefficients
        for s in range(k - 1, k + 1):
            u[s - 1, 0: (s - 1)] = dot(a[:, s - 1].T, b[:, 0: (s - 1)]) / m[0: (s - 1)]
            b[:, s - 1] = a[:, s - 1] - dot(b[:, 0: (s - 1)], u[s - 1, 0: (s - 1)].T)
            m[s - 1] = dot(b[:, s - 1], b[:, s - 1])

        if k > 2:
            k -= 1
            print("ppppppppppppp")
        else:
            # We have to do p/q, so do lstsq(q.T, p.T).T instead.
            p = dot(a[:, k:3].T, b[:, (k - 2): k])
            q = np.diag(m[(k - 2): k])  # type: ignore
            # pylint: disable=E1101
            result = np.linalg.lstsq(q.T, p.T, rcond=None)[0].T  # type: ignore
            u[k:3, (k - 2): k] = result


print(" ----------------- ")
print(a.T)
print("====================")
print(mapping.T)

print("\n\n\n")
print("====================")

k = 2
a = np.array(
    [[2.0, 0.0, 0.0],
     [0.1, 1.8, 0.0],
     [0.1, 0.2, 0.9]])

b = a.T
b[0, 0] = 3.0

m = [0.1, 0.2, 0.4]

p = dot(a[:, k:3].T, b[:, (k - 2): k])
q = np.diag(m[(k - 2): k])  # type: ignore

print(p)
print(q)

print("...................")
print("\n\n\n")
print("====================")