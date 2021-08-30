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
   use KindParamModule, only: IntKind, RealKind, CmplxKind
   use MathParamModule, only: ZERO, CZERO, SQRTm1, TEN2m6, TEN2m8, &
      HALF, ONE, TWO, FOUR, PI, PI2, PI4, Y0inv
   use ErrorHandlerModule, only: ErrorHandler, WarningHandler, &
      StopHandler
   use IntegerFactorsModule, only: lofk, mofk, lofj, mofj, kofj
   use LatticeModule, only: numlat
   use IntegerFactorsModule, only: initIntegerFactors, endIntegerFactors
   use MathParamModule, only: THIRD, THREE
   use SphericalHarmonicsModule, only: initSphericalHarmonics
   use GauntFactorsModule, only: isGauntInitialized => isInitialized
   use GauntFactorsModule, only: initGauntFactors
   use GauntFactorsModule, only: getGauntFactor
   use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_double
   implicit none
   private

   integer (kind = IntKind), parameter :: PrintLevel = 1

contains

   !
   ! @TODO: make a subroutine that is callable from C++
   !
   subroutine calculate_madelung_matrix(num_local_atoms, num_atoms, gindex, &
      lmax_rho_in, lmax_pot_in, bravais, posi, iprint)

      integer (kind = IntKind), intent(in) :: num_atoms
      integer (kind = IntKind), intent(in) :: num_local_atoms
      integer (kind = IntKind), intent(in) :: gindex(num_local_atoms)
      integer (kind = IntKind), intent(in) :: iprint
      !
      integer (kind = IntKind), intent(in) :: lmax_rho_in
      integer (kind = IntKind), intent(in) :: lmax_pot_in
      !
      integer (kind = IntKind) :: i, n1, l, m, kl, jl, jl_pot
      integer (kind = IntKind) :: l_pot, kl_pot, l_rho, kl_rho
      integer (kind = IntKind) :: l_sum, m_dif
      !
      real (kind = RealKind), intent(in) :: bravais(3, 3)
      real (kind = RealKind), target, intent(in) :: posi(3, num_atoms)
      !
      real (kind = RealKind), pointer :: factmat(:)
      real (kind = RealKind) :: gaunt
      real (kind = RealKind) :: sfac, eta, a0, alat
      real (kind = RealKind) :: factor, omegbra
      real (kind = RealKind) :: vatom
      real (kind = RealKind) :: vbrar(3, 3)
      real (kind = RealKind) :: vbrak(3, 3)

      integer (kind = IntKind) :: lmax_rho
      integer (kind = IntKind) :: jmax_rho
      integer (kind = IntKind) :: kmax_rho
      integer (kind = IntKind) :: lmax_pot
      integer (kind = IntKind) :: jmax_pot
      integer (kind = IntKind) :: kmax_pot
      integer (kind = IntKind) :: lmax_mad
      integer (kind = IntKind) :: jmax_mad
      integer (kind = IntKind) :: kmax_mad, nrslat, nknlat

      complex (kind=CmplxKind), allocatable, target :: cspace(:)
      complex (kind=CmplxKind), allocatable :: Ylm(:), ctmp(:)

      real (kind = RealKind), allocatable :: atom_position(:, :)

      real (kind = RealKind), allocatable, target :: madmat(:, :)
      real (kind = RealKind), allocatable, target :: rspace(:)

      real (kind = RealKind), allocatable :: rslat (:, :)
      real (kind = RealKind), allocatable :: rslatsq (:)

      real (kind = RealKind), allocatable :: knlat (:,:)
      real (kind = RealKind), allocatable :: knlatsq (:)

      complex (kind = CmplxKind), allocatable :: DL_matrix(:, :, :)
      complex (kind = CmplxKind), allocatable :: DL_factor(:, :)

      !
      !
      lmax_pot = lmax_pot_in
      jmax_pot = (lmax_pot + 1) * (lmax_pot + 2) / 2
      kmax_pot = (lmax_pot + 1) * (lmax_pot + 1)
      lmax_rho = lmax_rho_in
      jmax_rho = (lmax_rho + 1) * (lmax_rho + 2) / 2
      kmax_rho = (lmax_rho + 1) * (lmax_rho + 1)
      lmax_mad = 1 * max(lmax_pot_in, lmax_rho_in)
      jmax_mad = (lmax_mad + 1) * (lmax_mad + 2) / 2
      kmax_mad = (lmax_mad + 1) * (lmax_mad + 1)

      !
      !
      allocate(madmat(num_atoms, num_local_atoms))
      allocate(rspace(1:lmax_mad + 1))
      allocate(atom_position, source = posi)

      ! calculate the scaling factor of the bravais lattice
      call calScalingFactor(Bravais, rspace, lmax_mad, sfac, eta, nrslat, nknlat)

      a0 = sfac
      !
      vbrar(1:3, 1:3) = bravais(1:3, 1:3)

      !  change units so that both vbrar and atom_posi_* are in
      !  in the units of a0
      call dscal(9, ONE / a0, vbrar, 1)
      call dscal(num_atoms, ONE / a0, atom_position(1, 1), 3)
      call dscal(num_atoms, ONE / a0, atom_position(2, 1), 3)
      call dscal(num_atoms, ONE / a0, atom_position(3, 1), 3)
      !  -------------------------------------------------------------------
      !
      omegbra = (vbrar(2, 1) * vbrar(3, 2) - vbrar(3, 1) * vbrar(2, 2)) * vbrar(1, 3) + &
         (vbrar(3, 1) * vbrar(1, 2) - vbrar(1, 1) * vbrar(3, 2)) * vbrar(2, 3) + &
         (vbrar(1, 1) * vbrar(2, 2) - vbrar(2, 1) * vbrar(1, 2)) * vbrar(3, 3)

      factor = PI2 / omegbra
      omegbra = abs(omegbra)
      vatom = omegbra / real(num_atoms, RealKind)
      alat = a0 * (THREE * vatom / PI4)**THIRD

      !
      ! Generate basis vector for reciprocal lattice
      !
      vbrak(1, 1) = factor * (vbrar(2, 2) * vbrar(3, 3) - vbrar(3, 2) * vbrar(2, 3))
      vbrak(2, 1) = factor * (vbrar(3, 2) * vbrar(1, 3) - vbrar(1, 2) * vbrar(3, 3))
      vbrak(3, 1) = factor * (vbrar(1, 2) * vbrar(2, 3) - vbrar(2, 2) * vbrar(1, 3))
      vbrak(1, 2) = factor * (vbrar(2, 3) * vbrar(3, 1) - vbrar(3, 3) * vbrar(2, 1))
      vbrak(2, 2) = factor * (vbrar(3, 3) * vbrar(1, 1) - vbrar(1, 3) * vbrar(3, 1))
      vbrak(3, 2) = factor * (vbrar(1, 3) * vbrar(2, 1) - vbrar(2, 3) * vbrar(1, 1))
      vbrak(1, 3) = factor * (vbrar(2, 1) * vbrar(3, 2) - vbrar(3, 1) * vbrar(2, 2))
      vbrak(2, 3) = factor * (vbrar(3, 1) * vbrar(1, 2) - vbrar(1, 1) * vbrar(3, 2))
      vbrak(3, 3) = factor * (vbrar(1, 1) * vbrar(2, 2) - vbrar(2, 1) * vbrar(1, 2))


      call genLattice(vbrar, vbrak, rslat, rslatsq, knlat, knlatsq, nknlat, nrslat, &
            rspace, eta, lmax_mad)

      if(PrintLevel > 1) then
         write(6, '(/)')
         write(6, *) "Madelung:: scaling factor: ", sfac
         write(6, '(12x,a)')                                              &
            '    n                    rslat                  rslatsq'
         write(6, '(12x,56(''=''))')
         write(6, '(12x,1i5,2x,4f12.5)')                                  &
            (n1, rslat(n1, 1), rslat(n1, 2), rslat(n1, 3), rslatsq(n1), n1=1, nrslat)
         write(6, '(/)')
         write(6, '(12x,a)')                                              &
            '    n                    knlat                  knlatsq'
         write(6, '(12x,56(''=''))')
         write(6, '(12x,1i5,2x,4f12.5)')                                  &
            (n1, knlat(n1, 1), knlat(n1, 2), knlat(n1, 3), knlatsq(n1), n1=1, nknlat)
      endif

      if (jmax_mad > 1) then

         if (.not.isGauntInitialized()) then
            !        -------------------------------------------------------------
            call initSphericalHarmonics(lmax_mad * 2)
            call initGauntFactors(lmax_mad, 'xxxx', iprint)
            !        -------------------------------------------------------------
         endif
         !
         call initIntegerFactors(lmax_mad)

         allocate(DL_matrix(num_atoms, 1:kmax_mad, num_local_atoms))
         allocate(DL_factor(kmax_mad, jmax_mad))

         allocate(Ylm(1:kmax_mad), cspace(1:kmax_mad), ctmp(1:kmax_mad))
         DL_matrix=CZERO
         DL_factor=CZERO
         !
         factmat=>rspace(1:lmax_mad + 1)
         !
         factmat(1)=ONE
         do l=1, lmax_mad
            factmat(l + 1)=factmat(l) / (2 * l + ONE)
         enddo

         do jl_pot=1, jmax_mad
            l_pot=lofj(jl_pot)
            kl_pot=kofj(jl_pot)
            do kl_rho=1, kmax_mad
               l_rho=lofk(kl_rho)
               l_sum=l_pot + l_rho
               m_dif=mofk(kl_rho) - mofj(jl_pot)
               kl=(l_sum + 1) * (l_sum + 1) - l_sum + m_dif
               gaunt=getGauntFactor(kl_pot, kl_rho, kl)
               DL_factor(kl_rho, jl_pot)=gaunt * factmat(l_pot + 1) * factmat(l_rho + 1)
            enddo
         enddo
         nullify(factmat)
      end if

      !
      ! Generalized madlung matrix
      !


      !
      !
      !


      if (jmax_mad > 1) then
         deallocate(Ylm, cspace, ctmp)
      endif

      if (jmax_mad > 1) then
         call endIntegerFactors()
      endif

   end subroutine calculate_madelung_matrix

   subroutine calScalingFactor(Bravais, rspace, lmax_mad, sfac, eta, nrslat, nknlat)
      !
      real (kind = RealKind), intent(in) :: Bravais(3, 3)
      real (kind = RealKind), intent(in) :: rspace(:)
      integer (kind = IntKind), intent(in) :: lmax_mad
      real (kind = RealKind), intent(out) :: sfac
      real (kind = RealKind), intent(out) :: eta
      integer (kind = IntKind), intent(out) :: nrslat
      integer (kind = IntKind), intent(out) :: nknlat
      !
      integer (kind = IntKind) :: nm1, nm2, nm3, iter
      integer (kind = IntKind), parameter :: max_iter = 10000
      !
      real (kind = RealKind) :: a1, a2, a3
      real (kind = RealKind) :: vbrar(3, 3), vbrak(3, 3)
      real (kind = RealKind) :: volr, vfac
      real (kind = RealKind) :: rscut, kncut
      !
      real (kind = RealKind), parameter :: tfac = 1.2d0
      real (kind = RealKind), parameter :: fstep = 0.02d0
      !
      logical :: done = .false.
      !
      a1 = sqrt(Bravais(1, 1) * Bravais(1, 1) + Bravais(2, 1) * Bravais(2, 1) + &
         Bravais(3, 1) * Bravais(3, 1))
      a2 = sqrt(Bravais(1, 2) * Bravais(1, 2) + Bravais(2, 2) * Bravais(2, 2) + &
         Bravais(3, 2) * Bravais(3, 2))
      a3 = sqrt(Bravais(1, 3) * Bravais(1, 3) + Bravais(2, 3) * Bravais(2, 3) + &
         Bravais(3, 3) * Bravais(3, 3))
      !
      sfac = min(a1, a2, a3)
      eta = HALF + 0.1d0 * max(a1, a2, a3) / sfac
      !
      sfac = sfac / PI2
      !
      done = .false.
      iter = 0
      do while (.not. done)
         iter = iter + 1
         !     ================================================================
         !     scale the Bravais lattice and the reciprical lattice.
         !     ================================================================
         vbrar(1:3, 1:3) = Bravais(1:3, 1:3) / sfac
         volr = (vbrar(2, 1) * vbrar(3, 2) - vbrar(3, 1) * vbrar(2, 2)) * vbrar(1, 3) + &
            (vbrar(3, 1) * vbrar(1, 2) - vbrar(1, 1) * vbrar(3, 2)) * vbrar(2, 3) + &
            (vbrar(1, 1) * vbrar(2, 2) - vbrar(2, 1) * vbrar(1, 2)) * vbrar(3, 3)
         vfac = PI2 / abs(volr)
         vbrak(1, 1) = vfac * (vbrar(2, 2) * vbrar(3, 3) - vbrar(3, 2) * vbrar(2, 3))
         vbrak(2, 1) = vfac * (vbrar(3, 2) * vbrar(1, 3) - vbrar(1, 2) * vbrar(3, 3))
         vbrak(3, 1) = vfac * (vbrar(1, 2) * vbrar(2, 3) - vbrar(2, 2) * vbrar(1, 3))
         vbrak(1, 2) = vfac * (vbrar(2, 3) * vbrar(3, 1) - vbrar(3, 3) * vbrar(2, 1))
         vbrak(2, 2) = vfac * (vbrar(3, 3) * vbrar(1, 1) - vbrar(1, 3) * vbrar(3, 1))
         vbrak(3, 2) = vfac * (vbrar(1, 3) * vbrar(2, 1) - vbrar(2, 3) * vbrar(1, 1))
         vbrak(1, 3) = vfac * (vbrar(2, 1) * vbrar(3, 2) - vbrar(3, 1) * vbrar(2, 2))
         vbrak(2, 3) = vfac * (vbrar(3, 1) * vbrar(1, 2) - vbrar(1, 1) * vbrar(3, 2))
         vbrak(3, 3) = vfac * (vbrar(1, 1) * vbrar(2, 2) - vbrar(2, 1) * vbrar(1, 2))
         !
         !     ================================================================
         !     calculate rscut, the radius of real space truncation sphere.....
         !     ----------------------------------------------------------------

         call getrscut(vbrar(1:3, 1), vbrar(1:3, 2), vbrar(1:3, 3), &
            rspace(1:lmax_mad + 1), eta, lmax_mad, &
            nm1, nm2, nm3, rscut)
         call numlat(vbrar, rscut, nm1, nm2, nm3, nrslat)

         !     ----------------------------------------------------------------
         !
         !     ================================================================
         !     calculate rscut, the radius of real space truncation sphere.
         !     ----------------------------------------------------------------

         call getkncut(vbrak(1:3, 1), vbrak(1:3, 2), vbrak(1:3, 3), &
            eta, lmax_mad, &
            kncut, nm1, nm2, nm3)
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
            done = .true.
         else if (nknlat < nrslat / 2) then
            !        sfac = sfac/tfac
            sfac = sfac - fstep
         else if (nrslat < nknlat / 2) then
            sfac = sfac + fstep
         else
            done = .true.
         endif
      enddo
      !
   end subroutine calScalingFactor

   subroutine getrscut(a1, a2, a3, &
      gamma_l, eta, lmax_mad, &
      nm1, nm2, nm3, rscut)

      real (kind = RealKind), intent(in) :: a1(3)
      real (kind = RealKind), intent(in) :: a2(3)
      real (kind = RealKind), intent(in) :: a3(3)

      real (kind = RealKind), intent(in) :: eta
      real (kind = RealKind), intent(in) :: gamma_l(:)
      integer (kind = IntKind), intent(in) :: lmax_mad

      integer (kind = IntKind), intent(out) :: nm1
      integer (kind = IntKind), intent(out) :: nm2
      integer (kind = IntKind), intent(out) :: nm3
      real (kind = RealKind), intent(out) :: rscut

      integer (kind = IntKind) :: i, k, j
      real (kind = RealKind) :: r(3)
      real (kind = RealKind) :: rm
      real (kind = RealKind) :: term
      real (kind = RealKind), parameter :: epsi = 1.0d-14
      !
      !
      !  ===================================================================
      !  calculate nm1,nm2,nm3...........................................
      !  ===================================================================
      r(1) = sqrt(a1(1) * a1(1) + a1(2) * a1(2) + a1(3) * a1(3))
      term = ONE
      nm1 = 0
      do while(term.gt.HALF * epsi)
         nm1 = nm1 + 1
         rm = nm1 * r(1)
         !     term=erfc(rm/eta)/rm**lp1
         call calGammaFunc(rm / eta, lmax_mad, gamma_l(1:lmax_mad + 1))
         term = gamma_l(lmax_mad + 1) / (rm / TWO)**(lmax_mad + 1)
      enddo
      !
      r(2) = sqrt(a2(1) * a2(1) + a2(2) * a2(2) + a2(3) * a2(3))
      term = ONE
      nm2 = 0
      do while(term.gt.HALF * epsi)
         nm2 = nm2 + 1
         rm = nm2 * r(2)
         !     term=erfc(rm/eta)/rm**lp1
         call calGammaFunc(rm / eta, lmax_mad, gamma_l(1:lmax_mad + 1))
         term = gamma_l(lmax_mad + 1) / (rm / TWO)**(lmax_mad + 1)
      enddo
      !
      r(3) = sqrt(a3(1) * a3(1) + a3(2) * a3(2) + a3(3) * a3(3))
      term = ONE
      nm3 = 0
      do while(term.gt.HALF * epsi)
         nm3 = nm3 + 1
         rm = nm3 * r(3)
         !     term=erfc(rm/eta)/rm**lp1
         call calGammaFunc(rm / eta, lmax_mad, gamma_l(1:lmax_mad + 1))
         term = gamma_l(lmax_mad + 1) / (rm / TWO)**(lmax_mad + 1)
      enddo
      !
      !  ===================================================================
      !  calculate rscut.................................................
      !  ===================================================================
      rscut = r(1) * nm1
      do i = -1, 1
         r(1) = i * a1(1) * nm1
         r(2) = i * a1(2) * nm1
         r(3) = i * a1(3) * nm1
         do j = -1, 1
            r(1) = r(1) + j * a2(1) * nm2
            r(2) = r(2) + j * a2(2) * nm2
            r(3) = r(3) + j * a2(3) * nm2
            do k = -1, 1
               r(1) = r(1) + k * a3(1) * nm3
               r(2) = r(2) + k * a3(2) * nm3
               r(3) = r(3) + k * a3(3) * nm3
               rm = sqrt(r(1) * r(1) + r(2) * r(2) + r(3) * r(3))
               rscut = max(rscut, rm)
            enddo
         enddo
      enddo
      !
   end subroutine getrscut

   subroutine getkncut(a1, a2, a3, eta, lmax_mad, kncut, nm1, nm2, nm3)
      !
      integer (kind = IntKind), intent(out) :: nm1
      integer (kind = IntKind), intent(out) :: nm2
      integer (kind = IntKind), intent(out) :: nm3
      integer (kind = IntKind), intent(in) :: lmax_mad
      real (kind = RealKind), intent(in) :: eta
      integer (kind = IntKind) :: i, j, k
      !
      real (kind = RealKind), intent(in) :: a1(3)
      real (kind = RealKind), intent(in) :: a2(3)
      real (kind = RealKind), intent(in) :: a3(3)
      real (kind = RealKind), intent(out) :: kncut
      real (kind = RealKind) :: r(3)
      real (kind = RealKind) :: rm2
      real (kind = RealKind) :: term
      real (kind = RealKind) :: fac
      real (kind = RealKind), parameter :: epsi = 1.0d-14
      !
      fac = eta * eta / FOUR
      !
      !  ===================================================================
      !  calculate nm1,nm2,nm3...........................................
      !  ===================================================================
      r(1) = a1(1) * a1(1) + a1(2) * a1(2) + a1(3) * a1(3)
      term = ONE
      nm1 = 0
      do while(term.gt.HALF * epsi)
         nm1 = nm1 + 1
         rm2 = nm1 * nm1 * r(1)
         !     term=exp(-fac*rm2)/rm2
         term = exp(-fac * rm2) * (sqrt(rm2))**(lmax_mad - 2)
      enddo
      !
      r(2) = a2(1) * a2(1) + a2(2) * a2(2) + a2(3) * a2(3)
      term = ONE
      nm2 = 0
      do while(term.gt.HALF * epsi)
         nm2 = nm2 + 1
         rm2 = nm2 * nm2 * r(2)
         !     term=exp(-fac*rm2)/rm2
         term = exp(-fac * rm2) * (sqrt(rm2))**(lmax_mad - 2)
      enddo
      !
      r(3) = a3(1) * a3(1) + a3(2) * a3(2) + a3(3) * a3(3)
      term = ONE
      nm3 = 0
      do while(term.gt.HALF * epsi)
         nm3 = nm3 + 1
         rm2 = nm3 * nm3 * r(3)
         !     term=exp(-fac*rm2)/rm2
         term = exp(-fac * rm2) * (sqrt(rm2))**(lmax_mad - 2)
      enddo
      !
      !  ===================================================================
      !  calculate kncut.................................................
      !  ===================================================================
      kncut = sqrt(r(1)) * nm1
      do i = -1, 1
         r(1) = i * a1(1) * nm1
         r(2) = i * a1(2) * nm1
         r(3) = i * a1(3) * nm1
         do j = -1, 1
            r(1) = r(1) + j * a2(1) * nm2
            r(2) = r(2) + j * a2(2) * nm2
            r(3) = r(3) + j * a2(3) * nm2
            do k = -1, 1
               r(1) = r(1) + k * a3(1) * nm3
               r(2) = r(2) + k * a3(2) * nm3
               r(3) = r(3) + k * a3(3) * nm3
               rm2 = sqrt(r(1) * r(1) + r(2) * r(2) + r(3) * r(3))
               kncut = max(kncut, rm2)
            enddo
         enddo
      enddo
      !
   end subroutine getkncut

   subroutine genLattice(vbrar, vbrak, rslat, rslatsq, knlat, knlatsq, nknlat, nrslat, &
      rspace, eta, lmax_mad)
      use LatticeModule, only: createLattice
      !
      integer (kind = IntKind) :: nm1, nm2, nm3
      integer (kind = IntKind) :: nr, nk, ipmax
      !
      real (kind = RealKind), intent(in) :: vbrar(3, 3)
      real (kind = RealKind), intent(in) :: vbrak(3, 3)

      real (kind = RealKind), intent(inout), allocatable, target :: rslat(:, :)
      real (kind = RealKind), intent(inout), allocatable, target :: rslatsq(:)

      real (kind = RealKind), intent(inout), allocatable, target :: knlat(:, :)
      real (kind = RealKind), intent(inout), allocatable, target :: knlatsq(:)

      integer (kind = IntKind), intent(out) :: nknlat
      integer (kind = IntKind), intent(out) :: nrslat

      real (kind = RealKind), intent(in) :: rspace(:)
      real (kind = RealKind), intent(in) :: eta
      integer (kind = IntKind), intent(out) :: lmax_mad

      !
      real (kind = RealKind), pointer :: vec_x(:)
      real (kind = RealKind), pointer :: vec_y(:)
      real (kind = RealKind), pointer :: vec_z(:)
      real (kind = RealKind), pointer :: vecsq(:)
      !
      real (kind = RealKind) :: rscut
      real (kind = RealKind) :: kncut

      !
      !  *******************************************************************
      !  Sets up real space Bravais lattice vectors
      !  *******************************************************************
      !
      !  ===================================================================
      !  calculate rscut, the radius of real space truncation sphere.....
      !  -------------------------------------------------------------------
      call getrscut(vbrar(1:3, 1), vbrar(1:3, 2), vbrar(1:3, 3), &
         rspace(1:lmax_mad + 1), eta, lmax_mad, &
         nm1, nm2, nm3, rscut)

      call numlat(vbrar, rscut, nm1, nm2, nm3, nr)
      !  -------------------------------------------------------------------
      nrslat = nr
      allocate(rslat(1:nrslat, 1:3), rslatsq(1:nrslat))
      rslat = ZERO
      rslatsq = ZERO

      !
      !  ===================================================================
      !  generate the real space lattice vectors.........................
      !  ===================================================================
      vec_x => rslat(1:nrslat, 1)
      vec_y => rslat(1:nrslat, 2)
      vec_z => rslat(1:nrslat, 3)
      vecsq => rslatsq(1:nrslat)
      !
      ipmax = nrslat
      !  -------------------------------------------------------------------
      call createLattice(vbrar, rscut, nm1, nm2, nm3, vec_x, vec_y, vec_z, vecsq, nr, ipmax)
      !  -------------------------------------------------------------------
      !
      if(PrintLevel.ge.1) then
         write(6, '(/,'' Real Space Lattice:: nm1,nm2,nm3   = '',3i5)')nm1, nm2, nm3
         write(6, '(  ''                      Rs cut radius = '',1f10.5)') rscut
         write(6, '(  ''                      Number of Rs  = '',i5)') nrslat
      endif
      !
      !  *******************************************************************
      !  Sets up receprocal space Bravais lattice vectors
      !  *******************************************************************
      !
      !  ===================================================================
      !  calculate kncut, the radius of k-space truncation sphere........
      !  -------------------------------------------------------------------
      call getkncut(vbrak(1:3, 1), vbrak(1:3, 2), vbrak(1:3, 3), &
         eta, lmax_mad, kncut, nm1, nm2, nm3)
      call numlat(vbrak, kncut, nm1, nm2, nm3, nr)
      !  -------------------------------------------------------------------
      nknlat = nr
      allocate(knlat(1:nknlat, 1:3), knlatsq(1:nknlat))
      knlat = ZERO
      knlatsq = ZERO
      !
      !  ===================================================================
      !  generate the reciprocal space lattice vectors...................
      !  ===================================================================
      vec_x => knlat(1:nknlat, 1)
      vec_y => knlat(1:nknlat, 2)
      vec_z => knlat(1:nknlat, 3)
      vecsq => knlatsq(1:nknlat)
      !
      ipmax = nknlat
      !  -------------------------------------------------------------------
      call createLattice(vbrak, kncut, nm1, nm2, nm3, vec_x, vec_y, vec_z, vecsq, nk, ipmax)
      !  -------------------------------------------------------------------
      !
      if(PrintLevel.ge.1) then
         write(6, '(/,'' Reciprocal Lattice:: nm1,nm2,nm3   = '',3i5)')nm1, nm2, nm3
         write(6, '(  ''                      Kn cut radius = '',1f10.5)') kncut
         write(6, '(  ''                      Number of Kn  = '',i5)') nknlat
      endif

   end subroutine genLattice

   subroutine madewd(id, myatom, cspace, kmax_end, kmax_mad, a0, alat, eta, &
      GlobalNumAtoms, jmax_end, jmax_mad, omegbra, atom_position, madmat, DL_matrix)
      !  ===================================================================
      !
      !  *******************************************************************
      !     program for calculating the madelung constants
      !
      !     input :
      !     =====
      !             id = local atom index
      !             myatom = global atom index
      !
      !  *******************************************************************
#ifdef DIRECT_SUM
   use SphericalHarmonicsModule, only : calYlm
#endif
      !
      implicit none
      !
      integer (kind=IntKind), intent(in) :: id
      integer (kind=IntKind), intent(in) :: myatom
      complex (kind=CmplxKind), intent(in), target :: cspace(:)
      integer (kind=IntKind), intent(in) :: kmax_end
      integer (kind=IntKind), intent(in) :: kmax_mad
      real (kind=RealKind), intent(in) :: a0
      real (kind=RealKind), intent(in) :: alat
      real (kind=RealKind), intent(in) :: eta
      integer (kind=IntKind), intent(in) :: GlobalNumAtoms
      integer (kind=IntKind), intent(in) :: jmax_end
      integer (kind=IntKind), intent(in) :: jmax_mad
      real (kind=RealKind), intent(in) :: omegbra
      real (kind=RealKind), intent(in) :: atom_position(:,:)
      complex (kind=CmplxKind), intent(inout) :: madmat(:,:)
      complex (kind=CmplxKind), intent(inout) :: DL_matrix(:,:,:)

      !
      integer (kind=IntKind) :: ibegin
      integer (kind=IntKind) :: n
      integer (kind=IntKind) :: l, kl
      !
      real (kind=RealKind) :: aij(3)
      real (kind=RealKind) :: r0tm
      real (kind=RealKind) :: term0
      real (kind=RealKind) :: term12
#ifdef DIRECT_SUM
   real (kind=RealKind) :: vec(3), vm
   real (kind=RealKind), pointer :: factmat(:)
   complex (kind=CmplxKind), pointer :: dirsum(:)
#endif
      !
      complex (kind=CmplxKind), pointer :: dlm(:)
      !
      !  *******************************************************************
      !  calculate Madelung constant matrix:
      !
      !     for i <> j,
      !                                   2   2       -> ->
      !            4*pi          1    -eta *Kq /4 - i*Kq*aij
      !     M   =  ---- * sum  ----- e
      !      ij    tau    q<>0    2
      !                         Kq
      !
      !                          ->   ->                  2
      !                 1 - erf(|Rn + aij|/eta)     pi*eta
      !          + sum ------------------------- - ---------
      !             n         ->   ->                 tau
      !                      |Rn + aij|
      !
      !     for i = j,
      !                                   2   2
      !            4*pi          1    -eta *Kq /4
      !     M   =  ---- * sum  ----- e
      !      ii    tau    q<>0    2
      !                         Kq
      !
      !                                             2
      !                  1 - erf(Rn/eta)      pi*eta           2
      !          + sum  ----------------- - ---------- - --------------
      !            n<>0         Rn             tau        sqrt(pi)*eta
      !
      !     eta is the Ewald parameter;
      !     tau, atom_posi_*, rslat_* and knlat_* are in the units of a0;
      !     madmat is in the units of a0=alat;
      !  *******************************************************************
      !
      term0=-PI * eta * eta / omegbra
      !
      !  ===================================================================
      !  start madmat calculation
      !  ===================================================================
      do n=1, GlobalNumAtoms
         !     ================================================================
         !     aij is in the units of a0 ...................................
         !     ================================================================
         aij(1)=atom_position(1, myatom) - atom_position(1, n)
         aij(2)=atom_position(2, myatom) - atom_position(2, n)
         aij(3)=atom_position(3, myatom) - atom_position(3, n)
         !
         if(n .eq. myatom) then
            ibegin=2
            r0tm=-TWO / sqrt(PI) / eta
         else
            ibegin=1
            r0tm=ZERO
         endif
         !     ================================================================
         !     perform the reciprocal-space sum and real-space sum
         !     ----------------------------------------------------------------
         call madsum(ibegin, aij, term12)
         !     ----------------------------------------------------------------
         madmat(n, id)=term12 + r0tm + term0
         !     ================================================================
         !     The unit of madmat is finally resumed to a0 = alat
         !     ================================================================
         madmat(n, id)=madmat(n, id) / a0
         !
         if (jmax_mad > 1) then
            DL_matrix(n, 1, id)=madmat(n, id) * Y0inv
#ifdef DIRECT_SUM
         dirsum => cspace(1:kmax_mad)
         dirsum(:) = CZERO
         do i=nrslat,ibegin,-1
            vec(1) = aij(1) - rslat_x(i)
            vec(2) = aij(2) - rslat_y(i)
            vec(3) = aij(3) - rslat_z(i)
            vm = a0*sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))
!           ----------------------------------------------------------
            call calYlm(vec,lmax_mad,Ylm)
!           ----------------------------------------------------------
            do kl = kmax_mad, 2, -1
               dirsum(kl) = dirsum(kl) + Ylm(kl)/vm**(lofk(kl)+1)
            enddo
         enddo
!
         factmat => rspace(1:lmax_mad+1)
         factmat(1) = ONE
         do l = 1, lmax_mad
            factmat(l+1) = factmat(l)*(2*l+ONE)
         enddo
!
         do kl = 2, kmax_mad
            DL_matrix(n,kl,id) =                                         &
                      dirsum(kl)*PI4*factmat(lofk(kl))*alat**(lofk(kl))
         enddo
!
         nullify(dirsum, factmat)
#else
            dlm=>cspace(1:kmax_mad)
            dlm=CZERO
            !        -------------------------------------------------------------
            call dlsum(ibegin, aij, dlm)
            !        -------------------------------------------------------------
            do kl=2, kmax_mad
               l=lofk(kl)
               DL_matrix(n, kl, id)=dlm(kl) * (alat / a0)**l / a0
            enddo
            nullify(dlm)
#endif
         endif
      enddo                                  ! end do loop over n
      !
   end subroutine madewd

end module madelung_mod
