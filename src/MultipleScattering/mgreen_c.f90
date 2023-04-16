!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mgreen(zj_flag, kkrsz, lmax, jws, iprpts, &
   zlr_left, zlr_right, jlr, green, w1, &
   iprint, istop)
   !     ================================================================
   !
   implicit none
   !
   character  istop*32
   character  sname*20
   !
   integer    zj_flag
   integer    kkrsz
   integer    lmax
   integer    jws
   integer    iprpts
   integer    iprint
   integer    ir
   integer    i
   integer    l, m

   !
   complex*16 zlr_left(iprpts, 0:lmax)
   complex*16 zlr_right(iprpts, 0:lmax)
   complex*16 jlr(iprpts, 0:lmax)
   complex*16 w1(kkrsz, kkrsz)
   complex*16 green(jws)
   complex*16 wt, ctmp
   !
   parameter (sname = 'mgreen')
   !
   !     ****************************************************************
   !     calculates the Green's function.................................
   !     ****************************************************************
   !
   !     ================================================================
   !     calculate the green function for this (energy,species,sub-lat,
   !     spin) :: Spherical part only in this code.......................
   !     ----------------------------------------------------------------
   call zeroout(green, 2 * jws)
   !     ----------------------------------------------------------------
   ! backward sum trying to cancel singularities in high l's first.
   i = kkrsz
   do l = lmax, 0, -1
      wt = 0.d0
      do m = -l, l
         wt = wt + w1(i, i)
         i = i - 1
      enddo
      do ir = 1, jws
         if(zj_flag.eq.1) then
            ctmp = zlr_right(ir, l) * wt - jlr(ir, l) * (2 * l + 1)
         else
            ctmp = zlr_right(ir, l) * wt
         endif
         !green(ir) = green(ir) + zlr_left(ir, l) * ctmp
      enddo
   enddo

   i = 0
   do l = 0, lmax, 1
      do m = -l, l

         i = i + 1

         do ir = 1, jws
            green(ir) = green(ir) + zlr_left(ir, l) * (zlr_right(ir, l) * w1(i, i) - jlr(ir, l))
         end do


         end do
   end do



   !
   !     ================================================================
   if(istop.eq.sname) then
      call fstop(sname)
   else
      return
   endif
   !
end
