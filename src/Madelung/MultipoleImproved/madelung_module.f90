! ********************************************************************
! *  madelung_mod                                                  *
! *  ==============                                                  *
! *     Purpose: compute Madelung matrix                             *
! *                                                                  *
! *     Public Functions:                                            *
! *            initMadelung -- initialize the module                 *
! *            endMadelung  -- clean up the arrays allocated in the  *
! *                            module                                *
! *       getMadelungMatrix -- returns M_{0,0}(1:Na,i)               *
! *             getDLMatrix -- returns alat^l * DL(1:Na,1:jmax,i)    *
! *             getDLFactor -- returns the prefactor of DL matrix.   *
! *                            The Madelung matrix is a product of   *
! *                            this factor and DL matrix.            *
! *                                                                  *
! ********************************************************************

module madelung_mod
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO, SQRTm1, TEN2m6, TEN2m8, &
         HALF, ONE, TWO, FOUR, PI, PI2, PI4, Y0inv
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, &
         StopHandler
   use IntegerFactorsModule, only : lofk, mofk, lofj, mofj, kofj
   use LatticeModule, only : numlat
   use IntegerFactorsModule, only : initIntegerFactors, endIntegerFactors
   use MathParamModule, only : THIRD, THREE
   use SphericalHarmonicsModule, only : initSphericalHarmonics
   use GauntFactorsModule, only : isGauntInitialized=>isInitialized
   use GauntFactorsModule, only : initGauntFactors
   use GauntFactorsModule, only : getGauntFactor
   use, intrinsic :: iso_c_binding, only : c_ptr, c_loc, c_double
   implicit none
   private

contains

   !
   ! @TODO: make a subroutine that is callable from C++
   !
   subroutine calculate_madelung_matrix(num_local_atoms, num_atoms, gindex, &
         lmax_rho_in, lmax_pot_in, bravais, posi, iprint)

      integer (kind=IntKind), intent(in) :: num_atoms
      integer (kind=IntKind), intent(in) :: num_local_atoms
      integer (kind=IntKind), intent(in) :: gindex(num_local_atoms)
      integer (kind=IntKind), intent(in) :: iprint
      !
      integer (kind=IntKind), intent(in) :: lmax_rho_in
      integer (kind=IntKind), intent(in) :: lmax_pot_in
      !
      integer (kind=IntKind) :: i, n1, l, m, kl, jl, jl_pot
      integer (kind=IntKind) :: l_pot, kl_pot, l_rho, kl_rho
      integer (kind=IntKind) :: l_sum, m_dif
      !
      real (kind=RealKind), intent(in) :: bravais(3, 3)
      real (kind=RealKind), target, intent(in) :: posi(3, num_atoms)
      !
      real (kind=RealKind), pointer :: factmat(:)
      real (kind=RealKind) :: gaunt
      real (kind=RealKind) :: sfac
      real (kind=RealKind) :: factor
      real (kind=RealKind) :: vatom
      real (kind=RealKind) :: vbrar(3, 3)
      real (kind=RealKind) :: vbrak(3, 3)

      integer (kind=IntKind) :: lmax_rho=0
      integer (kind=IntKind) :: jmax_rho=1
      integer (kind=IntKind) :: kmax_rho=1
      integer (kind=IntKind) :: lmax_pot=0
      integer (kind=IntKind) :: jmax_pot=1
      integer (kind=IntKind) :: kmax_pot=1
      integer (kind=IntKind) :: lmax_mad=0
      integer (kind=IntKind) :: jmax_mad=1
      integer (kind=IntKind) :: kmax_mad=1

      real (kind=RealKind), allocatable, target :: madmat(:, :)
      real (kind=RealKind), allocatable, target :: rspace(:)

      !
      !
      lmax_pot=lmax_pot_in
      jmax_pot=(lmax_pot + 1) * (lmax_pot + 2) / 2
      kmax_pot=(lmax_pot + 1) * (lmax_pot + 1)
      lmax_rho=lmax_rho_in
      jmax_rho=(lmax_rho + 1) * (lmax_rho + 2) / 2
      kmax_rho=(lmax_rho + 1) * (lmax_rho + 1)
      lmax_mad=1 * max(lmax_pot_in, lmax_rho_in)
      jmax_mad=(lmax_mad + 1) * (lmax_mad + 2) / 2
      kmax_mad=(lmax_mad + 1) * (lmax_mad + 1)

      !
      !
      allocate(madmat(num_atoms, num_local_atoms))
      allocate(rspace(1:lmax_mad + 1))

      ! calculate the scaling factor of the bravais lattice




   end subroutine calculate_madelung_matrix



   subroutine calScalingFactor(Bravais, sfac, eta, nrslat, nknlat)
      !
      real (kind=RealKind), intent(in) :: Bravais(3, 3)
      real (kind=RealKind), intent(out) :: sfac
      real (kind=RealKind), intent(out) :: eta
      integer (kind=IntKind), intent(out) :: nrslat
      integer (kind=IntKind), intent(out) :: nknlat
      !
      integer (kind=IntKind) :: nm1, nm2, nm3, iter
      integer (kind=IntKind), parameter :: max_iter=10000
      !
      real (kind=RealKind) :: a1, a2, a3
      real (kind=RealKind) :: vbrar(3, 3), vbrak(3, 3)
      real (kind=RealKind) :: volr, vfac
      real (kind=RealKind) :: rscut, kncut
      !
      real (kind=RealKind), parameter :: tfac=1.2d0
      real (kind=RealKind), parameter :: fstep=0.02d0
      !
      logical :: done=.false.
      !
      a1=sqrt(Bravais(1, 1) * Bravais(1, 1) + Bravais(2, 1) * Bravais(2, 1) + &
         Bravais(3, 1) * Bravais(3, 1))
      a2=sqrt(Bravais(1, 2) * Bravais(1, 2) + Bravais(2, 2) * Bravais(2, 2) + &
         Bravais(3, 2) * Bravais(3, 2))
      a3=sqrt(Bravais(1, 3) * Bravais(1, 3) + Bravais(2, 3) * Bravais(2, 3) + &
         Bravais(3, 3) * Bravais(3, 3))
      !
      sfac=min(a1, a2, a3)
      eta=HALF + 0.1d0 * max(a1, a2, a3) / sfac
      !
      sfac=sfac / PI2
      !
      done=.false.
      iter=0
      do while (.not. done)
         iter=iter + 1
         !     ================================================================
         !     scale the Bravais lattice and the reciprical lattice.
         !     ================================================================
         vbrar(1:3, 1:3)=Bravais(1:3, 1:3) / sfac
         volr=(vbrar(2, 1) * vbrar(3, 2) - vbrar(3, 1) * vbrar(2, 2)) * vbrar(1, 3) + &
            (vbrar(3, 1) * vbrar(1, 2) - vbrar(1, 1) * vbrar(3, 2)) * vbrar(2, 3) + &
            (vbrar(1, 1) * vbrar(2, 2) - vbrar(2, 1) * vbrar(1, 2)) * vbrar(3, 3)
         vfac=PI2 / abs(volr)
         vbrak(1, 1)=vfac * (vbrar(2, 2) * vbrar(3, 3) - vbrar(3, 2) * vbrar(2, 3))
         vbrak(2, 1)=vfac * (vbrar(3, 2) * vbrar(1, 3) - vbrar(1, 2) * vbrar(3, 3))
         vbrak(3, 1)=vfac * (vbrar(1, 2) * vbrar(2, 3) - vbrar(2, 2) * vbrar(1, 3))
         vbrak(1, 2)=vfac * (vbrar(2, 3) * vbrar(3, 1) - vbrar(3, 3) * vbrar(2, 1))
         vbrak(2, 2)=vfac * (vbrar(3, 3) * vbrar(1, 1) - vbrar(1, 3) * vbrar(3, 1))
         vbrak(3, 2)=vfac * (vbrar(1, 3) * vbrar(2, 1) - vbrar(2, 3) * vbrar(1, 1))
         vbrak(1, 3)=vfac * (vbrar(2, 1) * vbrar(3, 2) - vbrar(3, 1) * vbrar(2, 2))
         vbrak(2, 3)=vfac * (vbrar(3, 1) * vbrar(1, 2) - vbrar(1, 1) * vbrar(3, 2))
         vbrak(3, 3)=vfac * (vbrar(1, 1) * vbrar(2, 2) - vbrar(2, 1) * vbrar(1, 2))
         !
         !     ================================================================
         !     calculate rscut, the radius of real space truncation sphere.....
         !     ----------------------------------------------------------------

         !call getrscut(vbrar(1:3, 1), vbrar(1:3, 2), vbrar(1:3, 3), rscut, nm1, nm2, nm3)
         call numlat(vbrar, rscut, nm1, nm2, nm3, nrslat)

         !     ----------------------------------------------------------------
         !
         !     ================================================================
         !     calculate rscut, the radius of real space truncation sphere.
         !     ----------------------------------------------------------------

         ! call getkncut(vbrak(1:3, 1), vbrak(1:3, 2), vbrak(1:3, 3), kncut, nm1, nm2, nm3)
         call numlat(vbrak, kncut, nm1, nm2, nm3, nknlat)

         !     ----------------------------------------------------------------
         !     write(6,'(a,3i8)')'iter, nrslat, nknlat = ',iter,nrslat,nknlat
         if (iter > max_iter .or. sfac <= 0.1d0) then
            write(6, '(a,3i8)')'iter, nrslat, nknlat = ', iter, nrslat, nknlat
            !        =============================================================
            !        If this message shows up, reduce fstep value.
            !        -------------------------------------------------------------
            call WarningHandler('calScalingFactor', &
               'The scaling factor may not be optimal', sfac)
            done=.true.
         else if (nknlat < nrslat / 2) then
            !        sfac = sfac/tfac
            sfac=sfac - fstep
         else if (nrslat < nknlat / 2) then
            sfac=sfac + fstep
         else
            done=.true.
         endif
      enddo
      !
   end subroutine calScalingFactor

   subroutine getrscut(a1, a2, a3, rscut, nm1, nm2, nm3)
      !
      integer (kind=IntKind), intent(out) :: nm1
      integer (kind=IntKind), intent(out) :: nm2
      integer (kind=IntKind), intent(out) :: nm3
      integer (kind=IntKind) :: i
      integer (kind=IntKind) :: j
      integer (kind=IntKind) :: k
      !
      real (kind=RealKind), intent(in) :: a1(3)
      real (kind=RealKind), intent(in) :: a2(3)
      real (kind=RealKind), intent(in) :: a3(3)
      real (kind=RealKind), intent(out) :: rscut
      real (kind=RealKind) :: r(3)
      real (kind=RealKind) :: rm
      real (kind=RealKind) :: term
      real (kind=RealKind), pointer :: gamma_l(:)
      real (kind=RealKind), parameter :: epsi=1.0d-14
      !
      gamma_l=>rspace(1:lmax_mad + 1)
      !
      !  ===================================================================
      !  calculate nm1,nm2,nm3...........................................
      !  ===================================================================
      r(1)=sqrt(a1(1) * a1(1) + a1(2) * a1(2) + a1(3) * a1(3))
      term=ONE
      nm1=0
      do while(term.gt.HALF * epsi)
         nm1=nm1 + 1
         rm=nm1 * r(1)
         !     term=erfc(rm/eta)/rm**lp1
         call calGammaFunc(rm / eta, lmax_mad, gamma_l(1:lmax_mad + 1))
         term=gamma_l(lmax_mad + 1) / (rm / TWO)**(lmax_mad + 1)
      enddo
      !
      r(2)=sqrt(a2(1) * a2(1) + a2(2) * a2(2) + a2(3) * a2(3))
      term=ONE
      nm2=0
      do while(term.gt.HALF * epsi)
         nm2=nm2 + 1
         rm=nm2 * r(2)
         !     term=erfc(rm/eta)/rm**lp1
         call calGammaFunc(rm / eta, lmax_mad, gamma_l(1:lmax_mad + 1))
         term=gamma_l(lmax_mad + 1) / (rm / TWO)**(lmax_mad + 1)
      enddo
      !
      r(3)=sqrt(a3(1) * a3(1) + a3(2) * a3(2) + a3(3) * a3(3))
      term=ONE
      nm3=0
      do while(term.gt.HALF * epsi)
         nm3=nm3 + 1
         rm=nm3 * r(3)
         !     term=erfc(rm/eta)/rm**lp1
         call calGammaFunc(rm / eta, lmax_mad, gamma_l(1:lmax_mad + 1))
         term=gamma_l(lmax_mad + 1) / (rm / TWO)**(lmax_mad + 1)
      enddo
      !
      !  ===================================================================
      !  calculate rscut.................................................
      !  ===================================================================
      rscut=r(1) * nm1
      do i=-1, 1
         r(1)=i * a1(1) * nm1
         r(2)=i * a1(2) * nm1
         r(3)=i * a1(3) * nm1
         do j=-1, 1
            r(1)=r(1) + j * a2(1) * nm2
            r(2)=r(2) + j * a2(2) * nm2
            r(3)=r(3) + j * a2(3) * nm2
            do k=-1, 1
               r(1)=r(1) + k * a3(1) * nm3
               r(2)=r(2) + k * a3(2) * nm3
               r(3)=r(3) + k * a3(3) * nm3
               rm=sqrt(r(1) * r(1) + r(2) * r(2) + r(3) * r(3))
               rscut=max(rscut, rm)
            enddo
         enddo
      enddo
      !
   end subroutine getrscut

end module madelung_mod
