/*
  This file was generated automatically with ./scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2020 (X86 64 LINUX)
  Maple source      : ./maple/lda_exc/lda_x_sloc.mpl
  Type of functional: lda_exc
*/

#define maple2c_order 4
#define MAPLE2C_FLAGS (XC_FLAGS_I_HAVE_EXC | XC_FLAGS_I_HAVE_VXC | XC_FLAGS_I_HAVE_FXC | XC_FLAGS_I_HAVE_KXC | XC_FLAGS_I_HAVE_LXC)


static inline void
func_unpol(const xc_func_type *p, int order, const double *rho , double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{

#ifndef XC_DONT_COMPILE_EXC
  double t1, t4, t5, t7, t8, t10;

#ifndef XC_DONT_COMPILE_FXC
  double t16, t17, t21;

#ifndef XC_DONT_COMPILE_KXC
  double t26, t27, t31;

#ifndef XC_DONT_COMPILE_LXC
  double t37, t46;
#endif

#endif

#endif

#endif


  lda_x_sloc_params *params;

  assert(p->params != NULL);
  params = (lda_x_sloc_params * )(p->params);

  t1 = params->b + 0.1e1;
  t4 = params->a / t1 / 0.2e1;
  t5 = pow(rho[0], params->b);
  t7 = pow(p->zeta_threshold, t1);
  t8 = my_piecewise3(0.1e1 <= p->zeta_threshold, t7, 1);
  t10 = t4 * t5 * t8;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = -0.2e1 * t10;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = -0.2e1 * t4 * t5 * params->b * t8 - 0.2e1 * t10;

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t16 = t4 * t5;
  t17 = 0.1e1 / rho[0];
  t21 = params->b * params->b;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = -0.2e1 * t16 * t21 * t17 * t8 - 0.2e1 * t16 * params->b * t17 * t8;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t26 = rho[0] * rho[0];
  t27 = 0.1e1 / t26;
  t31 = t21 * params->b;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = -0.2e1 * t16 * t31 * t27 * t8 + 0.2e1 * t16 * params->b * t27 * t8;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t37 = 0.1e1 / t26 / rho[0];
  t46 = t21 * t21;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = 0.2e1 * t16 * t21 * t37 * t8 + 0.4e1 * t16 * t31 * t37 * t8 - 0.2e1 * t16 * t46 * t37 * t8 - 0.4e1 * t16 * params->b * t37 * t8;

#ifndef XC_DONT_COMPILE_MXC

  if(order < 5) return;


#endif

#endif

#endif

#endif

#endif


}


static inline void
func_pol(const xc_func_type *p, int order, const double *rho , double *zk LDA_OUT_PARAMS_NO_EXC(XC_COMMA double *, ))
{

#ifndef XC_DONT_COMPILE_EXC
  double t1, t3, t4, t5, t6, t7, t8, t9;
  double t10, t11, t12, t13, t14, t15, t16, t17;
  double t18, t19;

#ifndef XC_DONT_COMPILE_VXC
  double t22, t24, t25, t26, t27, t28, t29, t30;
  double t31, t32, t35, t36, t37, t38, t41, t42;
  double t45, t48, t49, t52, t53;

#ifndef XC_DONT_COMPILE_FXC
  double t56, t57, t59, t61, t63, t64, t66, t68;
  double t70, t71, t72, t73, t74, t75, t78, t79;
  double t81, t86, t87, t88, t89, t90, t91, t93;
  double t98, t99, t103, t105, t106, t107, t114, t115;
  double t116, t123, t124, t129, t130, t133, t138, t139;
  double t140, t142, t147, t148;

#ifndef XC_DONT_COMPILE_KXC
  double t151, t153, t155, t158, t160, t161, t163, t165;
  double t168, t170, t171, t172, t174, t175, t177, t182;
  double t183, t184, t186, t194, t195, t196, t198, t199;
  double t201, t206, t214, t215, t221, t224, t226, t229;
  double t230, t231, t233, t234, t238, t239, t240, t249;
  double t250, t251, t258, t259, t260, t262, t263, t267;
  double t268, t269, t278, t279, t280, t287, t288, t294;
  double t296, t297, t298, t300, t305, t309, t313, t319;
  double t320, t321, t323, t328, t330, t334, t340, t341;
  double t348, t349, t351, t352, t358, t366, t367, t368;
  double t370, t371, t376, t384, t385;

#ifndef XC_DONT_COMPILE_LXC
  double t390, t394, t396, t399, t402, t404, t407, t411;
  double t413, t416, t419, t421, t422, t423, t424, t425;
  double t426, t428, t433, t434, t439, t446, t448, t459;
  double t460, t461, t462, t463, t464, t465, t467, t472;
  double t473, t478, t494, t495, t500, t503, t512, t514;
  double t515, t522, t523, t526, t527, t541, t554, t561;
  double t562, t564, t573, t575, t576, t583, t584, t587;
  double t588, t602, t615, t622, t623, t628, t631, t637;
  double t638, t640, t641, t643, t646, t649, t651, t654;
  double t655, t660, t662, t673, t675, t677, t680, t697;
  double t706, t708, t710, t712, t717, t719, t730, t732;
  double t734, t737, t754, t763, t765, t771, t779, t782;
  double t790, t792, t794, t812, t820, t824, t830, t840;
  double t841, t843, t861, t869, t877, t887, t888, t898;
  double t899, t902, t907, t908, t913, t919, t930, t931;
  double t932, t933, t936, t941, t942, t947, t963, t964;
#endif

#endif

#endif

#endif

#endif


  lda_x_sloc_params *params;

  assert(p->params != NULL);
  params = (lda_x_sloc_params * )(p->params);

  t1 = params->b + 0.1e1;
  t3 = 0.1e1 / t1 / 0.2e1;
  t4 = params->a * t3;
  t5 = rho[0] + rho[1];
  t6 = pow(t5, params->b);
  t7 = rho[0] - rho[1];
  t8 = 0.1e1 / t5;
  t9 = t7 * t8;
  t10 = 0.1e1 + t9;
  t11 = t10 <= p->zeta_threshold;
  t12 = pow(p->zeta_threshold, t1);
  t13 = pow(t10, t1);
  t14 = my_piecewise3(t11, t12, t13);
  t15 = 0.1e1 - t9;
  t16 = t15 <= p->zeta_threshold;
  t17 = pow(t15, t1);
  t18 = my_piecewise3(t16, t12, t17);
  t19 = t14 + t18;
  if(zk != NULL && (p->info->flags & XC_FLAGS_HAVE_EXC))
    zk[0] = -t4 * t6 * t19;

#ifndef XC_DONT_COMPILE_VXC

  if(order < 1) return;


  t22 = t6 * params->b;
  t24 = t4 * t22 * t19;
  t25 = t5 * params->a;
  t26 = t3 * t6;
  t27 = t13 * t1;
  t28 = t5 * t5;
  t29 = 0.1e1 / t28;
  t30 = t7 * t29;
  t31 = t8 - t30;
  t32 = 0.1e1 / t10;
  t35 = my_piecewise3(t11, 0, t27 * t31 * t32);
  t36 = t17 * t1;
  t37 = -t31;
  t38 = 0.1e1 / t15;
  t41 = my_piecewise3(t16, 0, t36 * t37 * t38);
  t42 = t35 + t41;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[0] = -t25 * t26 * t42 - t24 + (-t4 * t6 * t19);

  t45 = -t8 - t30;
  t48 = my_piecewise3(t11, 0, t27 * t45 * t32);
  t49 = -t45;
  t52 = my_piecewise3(t16, 0, t36 * t49 * t38);
  t53 = t48 + t52;
  if(vrho != NULL && (p->info->flags & XC_FLAGS_HAVE_VXC))
    vrho[1] = -t25 * t26 * t53 - t24 + (-t4 * t6 * t19);

#ifndef XC_DONT_COMPILE_FXC

  if(order < 2) return;


  t56 = t4 * t6;
  t57 = params->b * t8;
  t59 = t56 * t57 * t19;
  t61 = t4 * t6 * t42;
  t63 = params->b * params->b;
  t64 = t63 * t8;
  t66 = t56 * t64 * t19;
  t68 = t4 * t22 * t42;
  t70 = t1 * t1;
  t71 = t13 * t70;
  t72 = t31 * t31;
  t73 = t10 * t10;
  t74 = 0.1e1 / t73;
  t75 = t72 * t74;
  t78 = 0.1e1 / t28 / t5;
  t79 = t7 * t78;
  t81 = -0.2e1 * t29 + 0.2e1 * t79;
  t86 = my_piecewise3(t11, 0, t27 * t81 * t32 - t27 * t75 + t71 * t75);
  t87 = t17 * t70;
  t88 = t37 * t37;
  t89 = t15 * t15;
  t90 = 0.1e1 / t89;
  t91 = t88 * t90;
  t93 = -t81;
  t98 = my_piecewise3(t16, 0, t36 * t93 * t38 - t36 * t91 + t87 * t91);
  t99 = t86 + t98;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[0] = -t25 * t26 * t99 - t59 - 0.2e1 * t61 - t66 - 0.2e1 * t68;

  t103 = t4 * t6 * t53;
  t105 = t4 * t22 * t53;
  t106 = t31 * t74;
  t107 = t106 * t45;
  t114 = my_piecewise3(t11, 0, 0.2e1 * t27 * t79 * t32 - t27 * t107 + t71 * t107);
  t115 = t37 * t90;
  t116 = t115 * t49;
  t123 = my_piecewise3(t16, 0, -0.2e1 * t36 * t79 * t38 - t36 * t116 + t87 * t116);
  t124 = t114 + t123;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[1] = -t25 * t26 * t124 - t103 - t105 - t59 - t61 - t66 - t68;

  t129 = t45 * t45;
  t130 = t129 * t74;
  t133 = 0.2e1 * t29 + 0.2e1 * t79;
  t138 = my_piecewise3(t11, 0, t27 * t133 * t32 - t27 * t130 + t71 * t130);
  t139 = t49 * t49;
  t140 = t139 * t90;
  t142 = -t133;
  t147 = my_piecewise3(t16, 0, t36 * t142 * t38 - t36 * t140 + t87 * t140);
  t148 = t138 + t147;
  if(v2rho2 != NULL && (p->info->flags & XC_FLAGS_HAVE_FXC))
    v2rho2[2] = -t25 * t26 * t148 - 0.2e1 * t103 - 0.2e1 * t105 - t59 - t66;

#ifndef XC_DONT_COMPILE_KXC

  if(order < 3) return;


  t151 = params->b * t29;
  t153 = t56 * t151 * t19;
  t155 = t56 * t57 * t42;
  t158 = t4 * t6 * t99;
  t160 = t63 * params->b;
  t161 = t160 * t29;
  t163 = t56 * t161 * t19;
  t165 = t56 * t64 * t42;
  t168 = t4 * t22 * t99;
  t170 = t70 * t1;
  t171 = t13 * t170;
  t172 = t72 * t31;
  t174 = 0.1e1 / t73 / t10;
  t175 = t172 * t174;
  t177 = t106 * t81;
  t182 = t28 * t28;
  t183 = 0.1e1 / t182;
  t184 = t7 * t183;
  t186 = 0.6e1 * t78 - 0.6e1 * t184;
  t194 = my_piecewise3(t11, 0, t27 * t186 * t32 + t171 * t175 + 0.2e1 * t27 * t175 - 0.3e1 * t71 * t175 - 0.3e1 * t27 * t177 + 0.3e1 * t71 * t177);
  t195 = t17 * t170;
  t196 = t88 * t37;
  t198 = 0.1e1 / t89 / t15;
  t199 = t196 * t198;
  t201 = t115 * t93;
  t206 = -t186;
  t214 = my_piecewise3(t16, 0, t36 * t206 * t38 + t195 * t199 + 0.2e1 * t36 * t199 - 0.3e1 * t87 * t199 - 0.3e1 * t36 * t201 + 0.3e1 * t87 * t201);
  t215 = t194 + t214;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[0] = -t25 * t26 * t215 + t153 - 0.3e1 * t155 - 0.3e1 * t158 - t163 - 0.3e1 * t165 - 0.3e1 * t168;

  t221 = t56 * t57 * t53;
  t224 = 0.2e1 * t4 * t6 * t124;
  t226 = t56 * t64 * t53;
  t229 = 0.2e1 * t4 * t22 * t124;
  t230 = t72 * t174;
  t231 = t230 * t45;
  t233 = t81 * t74;
  t234 = t233 * t45;
  t238 = t71 * t31;
  t239 = t74 * t7;
  t240 = t239 * t78;
  t249 = t27 * t7;
  t250 = t78 * t74;
  t251 = t250 * t31;
  t258 = my_piecewise3(t11, 0, -0.6e1 * t27 * t184 * t32 + 0.2e1 * t27 * t78 * t32 + t171 * t231 + 0.2e1 * t27 * t231 - 0.3e1 * t71 * t231 - t27 * t234 + t71 * t234 + 0.4e1 * t238 * t240 - 0.4e1 * t249 * t251);
  t259 = t88 * t198;
  t260 = t259 * t49;
  t262 = t93 * t90;
  t263 = t262 * t49;
  t267 = t87 * t37;
  t268 = t90 * t7;
  t269 = t268 * t78;
  t278 = t36 * t7;
  t279 = t78 * t90;
  t280 = t279 * t37;
  t287 = my_piecewise3(t16, 0, 0.6e1 * t36 * t184 * t38 - 0.2e1 * t36 * t78 * t38 + t195 * t260 + 0.2e1 * t36 * t260 - 0.3e1 * t87 * t260 - t36 * t263 + t87 * t263 - 0.4e1 * t267 * t269 + 0.4e1 * t278 * t280);
  t288 = t258 + t287;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[1] = -t25 * t26 * t288 + t153 - 0.2e1 * t155 - t158 - t163 - 0.2e1 * t165 - t168 - t221 - t224 - t226 - t229;

  t294 = t4 * t6 * t148;
  t296 = t4 * t22 * t148;
  t297 = t31 * t174;
  t298 = t297 * t129;
  t300 = t71 * t45;
  t305 = t106 * t133;
  t309 = -0.2e1 * t78 - 0.6e1 * t184;
  t313 = t27 * t45;
  t319 = my_piecewise3(t11, 0, t27 * t309 * t32 + t171 * t298 + 0.4e1 * t300 * t240 - 0.4e1 * t313 * t240 + 0.2e1 * t27 * t298 - t27 * t305 - 0.3e1 * t71 * t298 + t71 * t305);
  t320 = t37 * t198;
  t321 = t320 * t139;
  t323 = t87 * t49;
  t328 = t115 * t142;
  t330 = -t309;
  t334 = t36 * t49;
  t340 = my_piecewise3(t16, 0, t36 * t330 * t38 + t195 * t321 - 0.4e1 * t323 * t269 + 0.4e1 * t334 * t269 + 0.2e1 * t36 * t321 - 0.3e1 * t87 * t321 - t36 * t328 + t87 * t328);
  t341 = t319 + t340;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[2] = -t25 * t26 * t341 + t153 - t155 - t163 - t165 - 0.2e1 * t221 - t224 - 0.2e1 * t226 - t229 - t294 - t296;

  t348 = t129 * t45;
  t349 = t348 * t174;
  t351 = t45 * t74;
  t352 = t351 * t133;
  t358 = -0.6e1 * t78 - 0.6e1 * t184;
  t366 = my_piecewise3(t11, 0, t27 * t358 * t32 + t171 * t349 + 0.2e1 * t27 * t349 - 0.3e1 * t27 * t352 - 0.3e1 * t71 * t349 + 0.3e1 * t71 * t352);
  t367 = t139 * t49;
  t368 = t367 * t198;
  t370 = t49 * t90;
  t371 = t370 * t142;
  t376 = -t358;
  t384 = my_piecewise3(t16, 0, t36 * t376 * t38 + t195 * t368 + 0.2e1 * t36 * t368 - 0.3e1 * t36 * t371 - 0.3e1 * t87 * t368 + 0.3e1 * t87 * t371);
  t385 = t366 + t384;
  if(v3rho3 != NULL && (p->info->flags & XC_FLAGS_HAVE_KXC))
    v3rho3[3] = -t25 * t26 * t385 + t153 - t163 - 0.3e1 * t221 - 0.3e1 * t226 - 0.3e1 * t294 - 0.3e1 * t296;

#ifndef XC_DONT_COMPILE_LXC

  if(order < 4) return;


  t390 = t56 * t63 * t78 * t19;
  t394 = 0.2e1 * t56 * params->b * t78 * t19;
  t396 = t56 * t151 * t42;
  t399 = t56 * t57 * t99;
  t402 = t4 * t6 * t215;
  t404 = t63 * t63;
  t407 = t56 * t404 * t78 * t19;
  t411 = 0.2e1 * t56 * t160 * t78 * t19;
  t413 = t56 * t161 * t42;
  t416 = t56 * t64 * t99;
  t419 = t4 * t22 * t215;
  t421 = t70 * t70;
  t422 = t13 * t421;
  t423 = t72 * t72;
  t424 = t73 * t73;
  t425 = 0.1e1 / t424;
  t426 = t423 * t425;
  t428 = t230 * t81;
  t433 = t81 * t81;
  t434 = t433 * t74;
  t439 = t106 * t186;
  t446 = t7 / t182 / t5;
  t448 = -0.24e2 * t183 + 0.24e2 * t446;
  t459 = t27 * t448 * t32 - 0.6e1 * t171 * t426 + 0.6e1 * t171 * t428 - 0.6e1 * t27 * t426 + 0.12e2 * t27 * t428 - 0.3e1 * t27 * t434 - 0.4e1 * t27 * t439 + t422 * t426 + 0.11e2 * t71 * t426 - 0.18e2 * t71 * t428 + 0.3e1 * t71 * t434 + 0.4e1 * t71 * t439;
  t460 = my_piecewise3(t11, 0, t459);
  t461 = t17 * t421;
  t462 = t88 * t88;
  t463 = t89 * t89;
  t464 = 0.1e1 / t463;
  t465 = t462 * t464;
  t467 = t259 * t93;
  t472 = t93 * t93;
  t473 = t472 * t90;
  t478 = t115 * t206;
  t494 = -t36 * t448 * t38 - 0.6e1 * t195 * t465 + 0.6e1 * t195 * t467 - 0.6e1 * t36 * t465 + 0.12e2 * t36 * t467 - 0.3e1 * t36 * t473 - 0.4e1 * t36 * t478 + t461 * t465 + 0.11e2 * t87 * t465 - 0.18e2 * t87 * t467 + 0.3e1 * t87 * t473 + 0.4e1 * t87 * t478;
  t495 = my_piecewise3(t16, 0, t494);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[0] = t390 - t394 + 0.4e1 * t396 - 0.6e1 * t399 - 0.4e1 * t402 - t407 + t411 - 0.4e1 * t413 - 0.6e1 * t416 - 0.4e1 * t419 - t25 * t26 * (t460 + t495);

  t500 = t4 * t22 * t288;
  t503 = t186 * t74 * t45;
  t512 = t172 * t425 * t45;
  t514 = t171 * t31;
  t515 = t174 * t45;
  t522 = t174 * t7;
  t523 = t522 * t78;
  t526 = t71 * t81;
  t527 = t515 * t31;
  t541 = 0.24e2 * t27 * t446 * t32;
  t554 = t239 * t183;
  t561 = 0.12e2 * t249 * t78 * t174 * t72 + 0.18e2 * t249 * t183 * t74 * t31 + 0.6e1 * t171 * t72 * t523 - 0.12e2 * t27 * t183 * t32 - 0.6e1 * t249 * t250 * t81 + 0.6e1 * t313 * t297 * t81 + 0.3e1 * t514 * t515 * t81 - 0.18e2 * t71 * t72 * t523 - 0.6e1 * t171 * t512 - 0.18e2 * t238 * t554 + 0.6e1 * t526 * t240 - 0.6e1 * t27 * t251 + 0.6e1 * t71 * t251 - t27 * t503 - 0.6e1 * t27 * t512 + t422 * t512 + t71 * t503 + 0.11e2 * t71 * t512 - 0.9e1 * t526 * t527 + t541;
  t562 = my_piecewise3(t11, 0, t561);
  t564 = t206 * t90 * t49;
  t573 = t196 * t464 * t49;
  t575 = t195 * t37;
  t576 = t198 * t49;
  t583 = t198 * t7;
  t584 = t583 * t78;
  t587 = t87 * t93;
  t588 = t576 * t37;
  t602 = 0.24e2 * t36 * t446 * t38;
  t615 = t268 * t183;
  t622 = -0.18e2 * t278 * t183 * t90 * t37 - 0.12e2 * t278 * t78 * t198 * t88 + 0.12e2 * t36 * t183 * t38 - 0.6e1 * t195 * t88 * t584 + 0.6e1 * t278 * t279 * t93 + 0.6e1 * t334 * t320 * t93 + 0.3e1 * t575 * t576 * t93 + 0.18e2 * t87 * t88 * t584 - 0.6e1 * t195 * t573 + 0.18e2 * t267 * t615 - 0.6e1 * t587 * t269 + 0.6e1 * t36 * t280 - 0.6e1 * t87 * t280 - t36 * t564 - 0.6e1 * t36 * t573 + t461 * t573 + t87 * t564 + 0.11e2 * t87 * t573 - 0.9e1 * t587 * t588 - t602;
  t623 = my_piecewise3(t16, 0, t622);
  t628 = t4 * t6 * t288;
  t631 = t56 * t161 * t53;
  t637 = t56 * t57 * t124;
  t638 = 0.3e1 * t637;
  t640 = t56 * t64 * t124;
  t641 = 0.3e1 * t640;
  t643 = t56 * t151 * t53;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[1] = -t419 - 0.3e1 * t500 - t25 * t26 * (t562 + t623) - t402 - 0.3e1 * t628 + t390 - t407 - t631 + 0.3e1 * t396 - 0.3e1 * t399 - 0.3e1 * t413 - 0.3e1 * t416 - t638 - t641 - t394 + t411 + t643;

  t646 = t4 * t22 * t341;
  t649 = t81 * t174 * t129;
  t651 = t7 * t7;
  t654 = t651 / t182 / t28;
  t655 = t654 * t74;
  t660 = t233 * t133;
  t662 = t106 * t309;
  t673 = t72 * t425 * t129;
  t675 = t230 * t133;
  t677 = t171 * t649 + t171 * t675 + 0.2e1 * t27 * t649 - 0.8e1 * t27 * t655 - t27 * t660 - 0.2e1 * t27 * t662 + t422 * t673 - 0.3e1 * t71 * t649 + 0.8e1 * t71 * t655 + t71 * t660 + 0.2e1 * t71 * t662 + t541;
  t680 = t351 * t78;
  t697 = t522 * t78 * t31;
  t706 = 0.8e1 * t514 * t515 * t79 - 0.6e1 * t171 * t673 - 0.6e1 * t27 * t673 + 0.2e1 * t27 * t675 - 0.4e1 * t27 * t680 - 0.12e2 * t300 * t554 - 0.24e2 * t300 * t697 + 0.12e2 * t313 * t554 + 0.16e2 * t313 * t697 + 0.11e2 * t71 * t673 - 0.3e1 * t71 * t675 + 0.4e1 * t71 * t680;
  t708 = my_piecewise3(t11, 0, t677 + t706);
  t710 = t93 * t198 * t139;
  t712 = t654 * t90;
  t717 = t262 * t142;
  t719 = t115 * t330;
  t730 = t88 * t464 * t139;
  t732 = t259 * t142;
  t734 = t195 * t710 + t195 * t732 + 0.2e1 * t36 * t710 - 0.8e1 * t36 * t712 - t36 * t717 - 0.2e1 * t36 * t719 + t461 * t730 - 0.3e1 * t87 * t710 + 0.8e1 * t87 * t712 + t87 * t717 + 0.2e1 * t87 * t719 - t602;
  t737 = t370 * t78;
  t754 = t583 * t78 * t37;
  t763 = -0.8e1 * t575 * t576 * t79 - 0.6e1 * t195 * t730 + 0.12e2 * t323 * t615 + 0.24e2 * t323 * t754 - 0.12e2 * t334 * t615 - 0.16e2 * t334 * t754 - 0.6e1 * t36 * t730 + 0.2e1 * t36 * t732 + 0.4e1 * t36 * t737 + 0.11e2 * t87 * t730 - 0.3e1 * t87 * t732 - 0.4e1 * t87 * t737;
  t765 = my_piecewise3(t16, 0, t734 + t763);
  t771 = t4 * t6 * t341;
  t779 = t56 * t57 * t148;
  t782 = t56 * t64 * t148;
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[2] = -0.2e1 * t500 - 0.2e1 * t646 - t25 * t26 * (t708 + t765) - 0.2e1 * t628 - 0.2e1 * t771 + t390 - t407 - 0.2e1 * t631 + 0.2e1 * t396 - t399 - 0.2e1 * t413 - t416 - 0.4e1 * t637 - 0.4e1 * t640 - t779 - t394 + t411 + 0.2e1 * t643 - t782;

  t790 = t4 * t6 * t385;
  t792 = t4 * t22 * t385;
  t794 = t31 * t425 * t348;
  t812 = t351 * t309;
  t820 = t106 * t358;
  t824 = 0.12e2 * t183 + 0.24e2 * t446;
  t830 = t27 * t133;
  t840 = -0.9e1 * t300 * t174 * t133 * t31 + 0.6e1 * t71 * t7 * t250 * t133 + 0.6e1 * t171 * t129 * t523 + 0.12e2 * t27 * t129 * t523 - 0.18e2 * t71 * t129 * t523 + 0.3e1 * t514 * t515 * t133 + t27 * t824 * t32 - 0.6e1 * t171 * t794 - 0.6e1 * t830 * t240 - 0.6e1 * t27 * t794 - 0.3e1 * t27 * t812 - t27 * t820 + t422 * t794 + 0.6e1 * t830 * t527 + 0.11e2 * t71 * t794 + 0.3e1 * t71 * t812 + t71 * t820;
  t841 = my_piecewise3(t11, 0, t840);
  t843 = t37 * t464 * t367;
  t861 = t370 * t330;
  t869 = t115 * t376;
  t877 = t36 * t142;
  t887 = -0.9e1 * t323 * t198 * t142 * t37 - 0.6e1 * t87 * t7 * t279 * t142 - 0.6e1 * t195 * t139 * t584 - 0.12e2 * t36 * t139 * t584 + 0.18e2 * t87 * t139 * t584 + 0.3e1 * t575 * t576 * t142 - t36 * t824 * t38 - 0.6e1 * t195 * t843 + 0.6e1 * t877 * t269 - 0.6e1 * t36 * t843 - 0.3e1 * t36 * t861 - t36 * t869 + t461 * t843 + 0.6e1 * t877 * t588 + 0.11e2 * t87 * t843 + 0.3e1 * t87 * t861 + t87 * t869;
  t888 = my_piecewise3(t16, 0, t887);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[3] = t390 - t394 + t396 + 0.3e1 * t643 - t638 - 0.3e1 * t779 - 0.3e1 * t771 - t407 + t411 - t413 - 0.3e1 * t631 - t641 - 0.3e1 * t782 - 0.3e1 * t646 - t790 - t792 - t25 * t26 * (t841 + t888);

  t898 = t129 * t129;
  t899 = t898 * t425;
  t902 = t129 * t174 * t133;
  t907 = t133 * t133;
  t908 = t907 * t74;
  t913 = t351 * t358;
  t919 = 0.24e2 * t183 + 0.24e2 * t446;
  t930 = t27 * t919 * t32 - 0.6e1 * t171 * t899 + 0.6e1 * t171 * t902 - 0.6e1 * t27 * t899 + 0.12e2 * t27 * t902 - 0.3e1 * t27 * t908 - 0.4e1 * t27 * t913 + t422 * t899 + 0.11e2 * t71 * t899 - 0.18e2 * t71 * t902 + 0.3e1 * t71 * t908 + 0.4e1 * t71 * t913;
  t931 = my_piecewise3(t11, 0, t930);
  t932 = t139 * t139;
  t933 = t932 * t464;
  t936 = t139 * t198 * t142;
  t941 = t142 * t142;
  t942 = t941 * t90;
  t947 = t370 * t376;
  t963 = -t36 * t919 * t38 - 0.6e1 * t195 * t933 + 0.6e1 * t195 * t936 - 0.6e1 * t36 * t933 + 0.12e2 * t36 * t936 - 0.3e1 * t36 * t942 - 0.4e1 * t36 * t947 + t461 * t933 + 0.11e2 * t87 * t933 - 0.18e2 * t87 * t936 + 0.3e1 * t87 * t942 + 0.4e1 * t87 * t947;
  t964 = my_piecewise3(t16, 0, t963);
  if(v4rho4 != NULL && (p->info->flags & XC_FLAGS_HAVE_LXC))
    v4rho4[4] = t390 - t394 + 0.4e1 * t643 - 0.6e1 * t779 - 0.4e1 * t790 - t407 + t411 - 0.4e1 * t631 - 0.6e1 * t782 - 0.4e1 * t792 - t25 * t26 * (t931 + t964);

#ifndef XC_DONT_COMPILE_MXC

  if(order < 5) return;


#endif

#endif

#endif

#endif

#endif


}
