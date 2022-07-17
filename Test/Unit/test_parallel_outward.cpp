//
// Created by F.Moitzi on 18.09.2021.
//

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <complex>
#include <vector>
#include <chrono>

#include "nrel_spherical_solver.hpp"

extern "C"
{
void rwave_(bool *srelflg, bool *irregflg,
            std::complex<double> *rr,
            std::complex<double> *ri,
            std::complex<double> *drr,
            std::complex<double> *dri,
            std::complex<double> *matom,
            int *l,
            std::complex<double> *energy,
            double *pot,
            int *nr,
            double *rad,
            std::complex<double> *potshift,
            int *jmt, // largest grid point
            int *kmax,
            double *emach,
            double *eunit,
            std::complex<double> *ws,
            int *wsdim,
            double *r_sph);
}

int main(int argc, char *argv[]) {

    int N = 1001;
    double rmin = 0.0001;
    double a = 2.7e6;
    double rmax = 2.95;
    int nsites = 1;
    int energy_points = 1;
    int lmax = 1;

    // Initialize all values
    std::vector<std::complex<double>> energies(energy_points, 0.1);
    std::vector<int> l_qns(lmax + 1);

    for (int l = 0; l <= lmax; l++) {
        l_qns[l] = l;
    }

    // Generate grid
    std::vector<double> R(N);
    std::vector<double> Rp(N);

    double beta = log(a) / (N - 2);
    double alpha = (rmax - rmin) / (exp(beta * (N - 1)) - 1.0);

    for (int i = 0; i < N; i++) {
        R[i] = alpha * (exp(beta * i) - 1.0) + rmin;
        Rp[i] = alpha * exp(beta * i) * beta;
    }

    // Generate potential
    std::vector<double> V(N);
    for (int i = 0; i < N; i++) {
        V[i] = -5.0 / R[i];
    }


    int total_points = nsites * energy_points * (1 + lmax);

    std::vector<std::vector<std::complex<double>>> P(total_points);
    std::vector<std::vector<std::complex<double>>> Q(total_points);
    std::vector<std::vector<std::complex<double>>> P_irr(total_points);
    std::vector<std::vector<std::complex<double>>> Q_irr(total_points);

    for (int ipoints = 0; ipoints < total_points; ipoints++) {
        P[ipoints].resize(N);
        Q[ipoints].resize(N);
        P_irr[ipoints].resize(N);
        Q_irr[ipoints].resize(N);
    }

    int ipoints = 0;

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();

//#pragma omp parallel for shared(P, Q, R, Rp, V, energies, N) collapse(3)
    for (int isite = 0; isite < nsites; isite++) {

        for (int ie = 0; ie < energy_points; ie++) {

            for (int l = 0; l <= lmax; l++) {

                nrel_spher_solver(5.0, energies[ie],
                                  l_qns[l], P[ipoints].data(), Q[ipoints].data(),
                                  R.data(), Rp.data(), V.data(), N, 1);

                ipoints++;
                break;

            }

        }

    }

    auto t2 = high_resolution_clock::now();

    /* Getting number of milliseconds as a double. */
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << ms_double.count() << "ms\n";

    auto td1 = high_resolution_clock::now();

    bool srelflg = true;
    bool irregflg = false;
    std::complex<double> matom{};
    int nr = 1;
    std::complex<double> potshift = 0.0;
    double emach = 1.0e-14;
    double eunit = 1.0;
    int ws_dim = 2 * (lmax + 1);
    int jmt = N - 5;
    int kmax = N - 1;
    double r_sph = R[jmt];
    std::complex<double> bes(2 * (lmax + 1));

    ipoints = 0;

    for (int isite = 0; isite < nsites; isite++) {

        for (int ie = 0; ie < energy_points; ie++) {

            for (int l = 0; l <= lmax; l++) {

                rwave_(&srelflg, &irregflg,
                       P[ipoints].data(),
                       P_irr[ipoints].data(),
                       Q[ipoints].data(),
                       P_irr[ipoints].data(),
                       &matom,
                       &l,
                       &energies[ie],
                       V.data(),
                       &nr,
                       R.data(),
                       &potshift,
                       &jmt, // largest grid point
                       &kmax,
                       &emach,
                       &eunit, &bes,
                       &ws_dim,
                       &r_sph);

                ipoints++;
                break;
            }

        }

    }

    auto td2 = high_resolution_clock::now();

    /* Getting number of milliseconds as a double. */
    duration<double, std::milli> ms_doubled = td2 - td1;

    std::cout << ms_doubled.count() << "ms\n";


    return EXIT_SUCCESS;
}