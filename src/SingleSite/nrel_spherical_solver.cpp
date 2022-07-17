#include "nrel_spherical_solver.hpp"

#include <complex>
#include <cmath>
#include <cstdlib>


/**
 * @brief Radial schroedinger equation
 *
 * @param r
 * @param rp
 * @param E
 * @param vpot
 * @param l
 * @param p
 * @param q
 * @param f
 */
static inline void cschroed_rad_func(const double r, const double rp,
                                     const std::complex<double> E, const double vpot,
                                     const int l, const std::complex<double> p,
                                     const std::complex<double> q,
                                     std::complex<double> f[2]) {
    f[0] = rp * q;
    f[1] = rp * 2 * ((vpot - E + l * (l + 1) / (2 * r * r)) * p);
}

/**
 * @brief Solves radial Schrödinger equation for the central potential
 *
 * @param Z
 * @param E
 * @param l
 * @param P
 * @param Q
 * @param R
 * @param Rp
 * @param V
 * @param end
 * @param sign_correction
 */
void nrel_spher_solver(double Z, std::complex<double> E, int l, std::complex<double> *P,
                       std::complex<double> *Q, const double *__restrict__ R,
                       const double *__restrict__ Rp, const double *__restrict__ V, int end,
                       int sign_correction) {
    int idx;
    std::complex<double> yp[2], f0[2], f1[2], f2[2], f3[2];
    std::complex<double> delta, c5, lam, r5, b5, i0, i1;

    /* Starting conditions
     *
     * - l == 0:
     *     P = r0 - Z*r0^2
     *     Q = 1 - 2*r0*Z
     *
     * - l != 0:
     *     P = r0^(l+1)
     *     Q = (l+1) * r0 ** l
     */

    if (l == 0) {
        if (Z <= 0.01) {
            P[0] = (2.0 + 2.0 * E * R[0] * R[0]) * pow(-1, sign_correction);
            Q[0] = E * R[0] * pow(-1, sign_correction);
        } else {
            P[0] = (R[0] - Z * R[0] * R[0]) * pow(-1, sign_correction);
            Q[0] = (1 - 2 * R[0] * Z) * pow(-1, sign_correction);
        }
    } else {
        P[0] = pow(R[0], (l + 1.0)) * pow(-1, sign_correction);
        Q[0] = (l + 1.0) * pow(R[0], l) * pow(-1, sign_correction);
    }

    // PC Euler
    cschroed_rad_func(R[0], Rp[0], E, V[0], l, P[0], Q[0], f0);
    yp[0] = P[0] + f0[0];
    yp[1] = Q[0] + f0[1];
    cschroed_rad_func(R[1], Rp[1], E, V[1], l, yp[0], yp[1], f1);
    P[1] = P[0] + (f0[0] + f1[0]) / 2.0;
    Q[1] = Q[0] + (f0[1] + f1[1]) / 2.0;

    // Adams Moulton 3
    cschroed_rad_func(R[1], Rp[1], E, V[1], l, P[1], Q[1], f1);

    i0 = P[1] + 1.0 / 12.0 * (-f0[0] + 8.0 * f1[0]);
    i1 = Q[1] + 1.0 / 12.0 * (-f0[1] + 8.0 * f1[1]);
    lam = 5.0 / 12.0;

    c5 = -2.0 * Rp[2] * (E - V[2] - l * (l + 1.0) / (2 * R[2] * R[2]));
    delta = 1.0 - lam * lam * Rp[2] * c5;

    P[2] = (i0 + lam * Rp[2] * i1) / delta;
    Q[2] = (lam * c5 * i0 + i1) / delta;

    // Adams Moulton 4
    cschroed_rad_func(R[2], Rp[2], E, V[2], l, P[2], Q[2], f2);

    i0 = P[2] + 1.0 / 24 * (f0[0] - 5.0 * f1[0] + 19.0 * f2[0]);
    i1 = Q[2] + 1.0 / 24 * (f0[1] - 5.0 * f1[1] + 19.0 * f2[1]);
    lam = 9.0 / 24.0;

    r5 = R[3];
    b5 = Rp[3];
    c5 = -2.0 * Rp[3] * (E - V[3] - l * (l + 1.0) / (2.0 * r5 * r5));
    delta = 1.0 - lam * lam * b5 * c5;

    P[3] = (i0 + lam * b5 * i1) / delta;
    Q[3] = (lam * c5 * i0 + i1) / delta;

    for (idx = 3; idx < end - 1; ++idx) {
        cschroed_rad_func(R[idx], Rp[idx], E, V[idx], l, P[idx], Q[idx], f3);

        i0 = P[idx] +
             1.0 / 720.0 *
             (-19.0 * f0[0] + 106.0 * f1[0] - 264.0 * f2[0] + 646.0 * f3[0]);

        i1 = Q[idx] +
             1.0 / 720.0 *
             (-19.0 * f0[1] + 106.0 * f1[1] - 264.0 * f2[1] + 646.0 * f3[1]);
        lam = 251.0 / 720.0;

        r5 = R[idx + 1];
        b5 = Rp[idx + 1];

        c5 = -2.0 * Rp[idx + 1] * (E - V[idx + 1] - l * (l + 1.0) / (2.0 * r5 * r5));

        delta = 1.0 - lam * lam * b5 * c5;

        P[idx + 1] = (i0 + lam * b5 * i1) / delta;
        Q[idx + 1] = (lam * c5 * i0 + i1) / delta;

        f0[0] = f1[0];
        f0[1] = f1[1];

        f1[0] = f2[0];
        f1[1] = f2[1];

        f2[0] = f3[0];
        f2[1] = f3[1];
    }
}
