import numpy as np

from scipy.integrate import quad

R = 0.5
kappa = 1.0 - 0.1j
nu = 0.1
l = -1


def func(chi):
    return np.real(np.exp(-R ** 2 * chi ** 2 + kappa ** 2 / (4 * chi ** 2)))


res = quad(func, 0.5 * np.sqrt(nu), np.inf, epsabs=1.0e-12, epsrel=1.0e-12, limit=1000)
print(res)


def func(chi):
    return np.imag(np.exp(-R ** 2 * chi ** 2 + kappa ** 2 / (4 * chi ** 2)))


res = quad(func, 0.5 * np.sqrt(nu), np.inf, epsabs=1.0e-12, epsrel=1.0e-12, limit=1000)
print(res)


def func(chi):
    return chi ** (2 * l) * np.real(np.exp(-R ** 2 * chi ** 2 + kappa ** 2 / (4 * chi ** 2)))


res = quad(func, 0.5 * np.sqrt(nu), np.inf, epsabs=1.0e-12, epsrel=1.0e-12, limit=1000)
print(res)


def func(chi):
    return chi ** (2 * l) * np.imag(np.exp(-R ** 2 * chi ** 2 + kappa ** 2 / (4 * chi ** 2)))


res = quad(func, 0.5 * np.sqrt(nu), np.inf, epsabs=1.0e-12, epsrel=1.0e-12, limit=1000)
print(res)

l = 2
def func(chi):
    return chi ** (2 * l) * np.real(np.exp(-R ** 2 * chi ** 2 + kappa ** 2 / (4 * chi ** 2)))


res = quad(func, 0.5 * np.sqrt(nu), np.inf, epsabs=1.0e-12, epsrel=1.0e-12, limit=1000)
print(res)


def func(chi):
    return chi ** (2 * l) * np.imag(np.exp(-R ** 2 * chi ** 2 + kappa ** 2 / (4 * chi ** 2)))


res = quad(func, 0.5 * np.sqrt(nu), np.inf, epsabs=1.0e-12, epsrel=1.0e-12, limit=1000)
print(res)
