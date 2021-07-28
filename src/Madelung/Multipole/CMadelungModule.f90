module CMadelungModule
   use iso_c_binding, only : c_int, c_double
   implicit none
   private
   public :: c_initMadelung

contains

   subroutine c_initMadelung(num_local_atoms, num_atoms, gindex, &
         lmax_rho, lmax_pot, bravais, posi, iprint) bind(C, name="initMadelung")
      use MadelungModule, only : initMadelung

      integer (kind=c_int), value, intent(in) :: num_atoms !< Local number of atoms
      integer (kind=c_int), value, intent(in) :: num_local_atoms !< Global number of atoms
      integer (kind=c_int), intent(in) :: gindex(num_local_atoms) !< Global indexing
      integer (kind=c_int), value, intent(in) :: lmax_rho
      integer (kind=c_int), value, intent(in) :: lmax_pot
      real (kind=c_double), intent(in) :: bravais(3, 3)
      real (kind=c_double), intent(in) :: posi(3, num_atoms)
      integer (kind=c_int), value, intent(in) :: iprint

      call initMadelung(num_local_atoms, num_atoms, gindex, lmax_rho, lmax_pot, bravais, posi, iprint)

   end subroutine

end module CMadelungModule