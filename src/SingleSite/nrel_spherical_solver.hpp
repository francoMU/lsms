
#ifndef LSMS_NREL_SPHERICAL_SOLVER_HPP
#define LSMS_NREL_SPHERICAL_SOLVER_HPP

#include <complex>

void nrel_spher_solver(double Z, std::complex<double> E, int l, std::complex<double> *P,
                       std::complex<double> *Q, const double *R,
                       const double *Rp, const double *V, int end,
                       int sign_correction);

#endif //LSMS_NREL_SPHERICAL_SOLVER_HPP
