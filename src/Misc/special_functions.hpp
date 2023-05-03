
#include <complex>

namespace lsms {

void sph_bessel_jl(std::complex<double> z, int lmax,
                   std::complex<double> *jl, std::complex<double> *jlp) {

  int l = 0;

  int n_max, k;
  std::complex<double> tk, cf1, dcf1, den, c, d, z_1;
  double dummy;

  constexpr double eps = std::numeric_limits<double>::epsilon();
  constexpr double min = std::numeric_limits<double>::min();
  constexpr int MAX_BESSEL_ITER = 2000;

  if (abs(z) < eps) {

    for (l = 0; l <= lmax; l++) {
      jl[l] = 0.0;
      jlp[l] = 0.0;
    }

    jl[0] = 1.0;
    jlp[1] = 1.0;
    return;
  }

  z_1 = 1.0 / z;

//  // Modified Steed's algorithm
  tk = z_1 * 2.0 * (double) n_max * z_1 + z_1 * 3.0;
  cf1 = z_1 * (double) n_max;

  if (std::abs(cf1) < small) {
    cf1 = small;
  }

  den = 1.0;

  c = cf1;
  d = 0.0;


  for(k = 0; k < MAX_BESSEL_ITER; k++) {

    c = tk - 1.0 / c;
    d = tk - d;

    if (abs(c) < small) then
          c = small
      endif
      if (abs(d) < small) then
            d = small
        endif

            d = 1 / d;
        dcf1 = d * c;
    cf1 = cf1 * dcf1;

    if (real(d) < 0.0_dp) then
          den = -den;
      endif

      if (abs(dcf1 - 1.0_dp) < eps) then
            break
        endif

            tk = tk + 2 * z_1;
  }


}

}
//
//
//   !-----------------------------------------------------------------------------------
//   !> Complex version of the same function.
//   !-----------------------------------------------------------------------------------
//   pure subroutine sph_bessel_jl_c(z, l_max, jl, jl_p)
//      use constants, only: nan64
//
//      !> Argument
//      complex(dp), intent(in) :: z
//
//      !> Maximum order `l`
//      integer, intent(in) :: l_max
//
//      !> Resulting values of j_l(z)
//      complex(dp), intent(out) :: jl(0:)
//
//      !> Resulting values of j'_l(z)
//      complex(dp), intent(out), optional :: jl_p(0:)
//
//      ! Local variables
//      real(dp), parameter :: small = sqrt(tiny(1.0_dp))
//      real(dp), parameter :: eps = epsilon(0.01_dp)
//      integer :: n_max, k
//      complex(kind=dp) :: tk, cf1, dcf1, den, c, d, z_1
//      real(kind=dp) :: dummy
//
//      !
//      !--- Extract the actual length of the output arrays.
//      !    The calculation is done up to `l = min(n_max-1, l_max)`.
//      !    This is done to avoid extra checks on input.
//      !
//      if (present(jl_p)) then
//         n_max = min(l_max, ubound(jl, 1), ubound(jl_p, 1))
//      else
//         n_max = min(l_max, ubound(jl, 1))
//      endif
//
//      !
//      !--- Quick return for zero argument
//      !
//      if (abs(z) < eps) then
//         jl = 0.0_dp
//         jl(0) = 1.0_dp
//         if (present(jl_p)) then
//            jl_p = 0.0_dp
//            jl_p(1) = 1.0_dp
//         endif
//         return
//      endif
//
//      z_1 = 1 / z
//
//      !
//      !--- Modified Steed's algorithm
//      !
//      tk = 2 * n_max * z_1 + 3 * z_1
//      cf1 = n_max * z_1
//
//      if (abs(cf1) < small) then
//         cf1 = small
//      endif
//
//      den = 1.0_dp
//
//      c = cf1
//      d = 0.0_dp
//
//      do k = 1, MAX_BESS_ITER
//         c = tk - 1 / c
//         d = tk - d
//
//         if (abs(c) < small) then
//            c = small
//         endif
//         if (abs(d) < small) then
//            d = small
//         endif
//
//         d = 1 / d
//         dcf1 = d * c
//         cf1 = cf1 * dcf1
//
//         if (real(d) < 0.0_dp) then
//            den = -den
//         endif
//
//         if (abs(dcf1 - 1.0_dp) < eps) then
//            exit
//         endif
//
//         tk = tk + 2 * z_1
//      enddo
//
//      !
//      !--- If the iterations have not converged, return NaN
//      !
//      if (k > MAX_BESS_ITER) then
//         jl = (nan64, nan64)
//         if (present(jl_p)) then
//            jl_p = (nan64, nan64)
//         endif
//         dummy = nan64 !! To silence the inappropriate warning in gfortran
//         return
//      endif
//
//      !
//      !--- Get the values at `l = n_max`
//      !
//      jl(n_max) = den
//
//      if (present(jl_p)) then
//         jl_p(n_max) = cf1 * den
//
//         !
//         !--- Downward recursion
//         !
//         c = n_max * z_1
//         do k = n_max, 1, -1
//            jl(k - 1) = (c + z_1) * jl(k) + jl_p(k)
//            c = c - z_1
//            jl_p(k - 1) = c * jl(k - 1) - jl(k)
//         enddo
//      else
//         tk = cf1 * den
//
//         !
//         !--- Downward recursion
//         !
//         c = n_max * z_1
//         do k = n_max, 1, -1
//            jl(k - 1) = (c + z_1) * jl(k) + tk
//            c = c - z_1
//            tk = c * jl(k - 1) - jl(k)
//         enddo
//      endif
//
//      !
//      !--- Rescale the entire series to correct for the
//      !    accumulated error.
//      !
//      den = jl(0)
//
//      d = sin(z) * z_1 / den
//      jl = jl * d
//
//      if (present(jl_p)) then
//         jl_p = jl_p * d
//      endif
//   end subroutine sph_bessel_jl_c