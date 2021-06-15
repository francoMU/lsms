/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_CORESTATES_H
#define LSMS_CORESTATES_H

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

void getCoreStates(LSMSSystemParameters &lsms, AtomData &atom);


extern "C"
{
void deepst_(int *nqn, int *lqn, int *kqn, Real *en, Real *rv, Real *r, Real *rf, Real *h, Real *z, Real *c,
             int *nitmax, Real *tol, int *nws, int *nlast, int *iter, int *iprpts, int *ipdeq);
void semcst_(int *nqn, int *lqn, int *kqn, Real *en, Real *rv, Real *r, Real *rf, Real *h, Real *z, Real *c,
             int *nitmax, Real *tol, int *nmt, int *nws, int *nlast, int *iter, int *iprpts, int *ipdeq);
void newint_(int *nr, Real *r, Real *f, Real *g, int *ip0);


/*
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     subroutine getcor(IN n_spin_pola,IN mtasa,
!    >                  IN jmt,IN jws,IN r_mesh,IN h,IN xstart,IN vr,
!    >                  IN numc,IN nc,IN lc,IN kc,INOUT ec,
!    >                  IN ztotss,IN zsemss,IN zcorss,
!    >                  OUT ecorv,OUT esemv,OUT corden,OUT semcor,
!    >                  IN nrelc,
!    >                  OUT qcpsc_mt,OUT qcpsc_ws,OUT mcpsc_mt,OUT mcpsc_ws,
!    >                  IN iprpts, IN ipcore,
!    >                  IN iprint,IN istop)
c     ================================================================
*/

void getcor_(int *n_spin_pola, int *mtasa,
             int *jmt, int *jws, double *r_mesh, double *h, double *xstart, double *vr,
             int *numc, int *nc, int *lc, int *kc, double *ec,
             double *ztotss, double *zsemss, double *zcorss,
             double *ecorv, double *esemv, double *corden, double *semcor,
             int *nrelc,
             double *qcpsc_mt, double *qcpsc_ws, double *mcpsc_mt, double *mcpsc_ws,
             int *iprpts, int *ipcore,
             int *iprint, char *istop, int itop_length);
}

#endif
