/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_LDA_C_VWN_4   31   /* Vosko, Wilk, & Nusair (4)   */

#include "maple2c/lda_c_vwn_4.c"

#define func maple2c_func
#include "work_lda.c"

const xc_func_info_type xc_func_info_lda_c_vwn_4 = {
  XC_LDA_C_VWN_4,
  XC_CORRELATION,
  "Vosko, Wilk & Nusair (VWN4)",
  XC_FAMILY_LDA,
  {&xc_ref_Vosko1980_1200, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC | XC_FLAGS_HAVE_FXC | XC_FLAGS_HAVE_KXC,
  1e-24,
  0, NULL, NULL,
  NULL, NULL,
  work_lda, NULL, NULL
};