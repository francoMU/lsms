
#ifndef MUST_MULTIPOLE_MADELUNG_H
#define MUST_MULTIPOLE_MADELUNG_H


#if defined (__cplusplus)
extern "C" {
#endif

void initMadelung(int num_atoms,
                  int num_loca_atoms,
                  int *gindex,
                  int lmax_rho,
                  int lmax_pot,
                  double *bravais,
                  double *pos,
                  int iprint);

#if defined (__cplusplus)
}
#endif


#endif //MUST_MULTIPOLE_MADELUNG_H
