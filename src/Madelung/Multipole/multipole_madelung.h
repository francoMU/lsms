
#ifndef MUST_MULTIPOLE_MADELUNG_H
#define MUST_MULTIPOLE_MADELUNG_H

#include <complex>

// Reduced madelung matrix
extern "C" std::complex<double> *DL_matrix_ptr;
// Prefactor of madelung matrix
extern "C" double *DL_factor_ptr;
// Scaling factor
extern "C" double *alat_ptr;

extern "C" double variable_2;

#if defined (__cplusplus)
extern "C" {
#endif

void test_calling_from_cpp();

void testMatrix(double **matrix);

void initMadelung(int num_atoms,
                  int num_loca_atoms,
                  int *gindex,
                  int lmax_rho,
                  int lmax_pot,
                  double *bravais,
                  double *pos,
                  int iprint);

void getMadelungMatrix(int local_atom_index,
                       double **madelung_matrix,
                       int *madelung_matrix_size);

void getDLMatrix(int local_atom_index,
                 std::complex<double> **matrix,
                 int *atom_size,
                 int *k_size,
                 double *factor);


void getDLFactor(int j_index,
                 double **matrix,
                 int *matrix_size);

void endMadelung();

#if defined (__cplusplus)
}
#endif


#endif //MUST_MULTIPOLE_MADELUNG_H
