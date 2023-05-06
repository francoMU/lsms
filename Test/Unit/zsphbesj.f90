module sp_mod
   use ISO_C_BINDING

contains

   subroutine zsphbesj(lmax, x, as, ac, x1, bj, bscrat) bind(c)
      !_______________________________________________________________________
      ! calculate spherical bessel functions j_l
      ! input: lmax integer scalar, max l for j_l
      !        x    compex*16 scalar, argument of j_l
      ! returns: as complex*16 scalar, sin(x)
      !          ac complex*16 scalar, cos(x)
      !          x1 complex*16 scalar, 1/x
      !          bj complex*16 array of (0:lmax), j_l
      !          bscrat complex*16 array of (0:lmax), scratch space

      !  xgz ornl 1994

      implicit complex(kind=8)(a-h, o-z)
      complex(kind=8) bj(0:lmax), bscrat(0:lmax)

      as = sin(x)
      ac = cos(x)
      x1 = 1.0d0 / x
      bj(0) = as * x1
      !   Forward recursion for small L
      1 k = ishft(int(0.75d0 * (abs(dreal(x)) + abs(dimag(x)))) + 2, -2)
      bscrat(0) = as
      if(k>=1) then
         if(k>lmax) k = lmax
         bscrat(1) = bj(0) - ac
         bj(1) = bscrat(1) * x1
         do l = 2, k
            bscrat(l) = (2 * l - 1) * bj(l - 1) - bscrat(l - 2)
            bj(l) = bscrat(l) * x1
         end do
         if(k==lmax) return
      endif
      !   Backward recursion from very large L down to L=k
      a = bscrat(k)
      nm = ishft(lmax, 4)
      aj = 2 * nm + 3
      x2 = x * x
      do l = nm, lmax + 2, -1
         aj = (2 * l + 1) - x2 / aj
      end do
      bscrat(lmax) = (2 * lmax + 3) * aj * x1 - x
      bj(lmax) = bscrat(lmax) * x1
      bscrat(lmax - 1) = (2 * lmax + 1) * bj(lmax) - aj
      do l = lmax - 1, k + 1, -1
         bj(l) = bscrat(l) * x1
         bscrat(l - 1) = (2 * l + 1) * bj(l) - bscrat(l + 1)
      end do
      ! scale to give correct bj(k)
      aj = a / bscrat(k)
      do l = k + 1, lmax
         bj(l) = aj * bj(l)
      end do
      return
   end

end module sp_mod