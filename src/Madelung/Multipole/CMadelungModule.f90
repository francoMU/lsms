module CMadelungModule
   use iso_c_binding, only : c_int, c_double, c_ptr, c_loc
   implicit none
   private
   public :: c_initMadelung
   public :: test_calling_from_cpp

   real(c_double) :: variable
   real(c_double), bind(C) :: variable_2

contains

   subroutine test_calling_from_cpp() bind(C, name="test_calling_from_cpp")

      print *, "Method was called"

      variable=1
      variable_2=2

   end subroutine test_calling_from_cpp

   subroutine testMatrix(matrix) bind(C, name="testMatrix")

      type(c_ptr), intent(out) :: matrix

      real(kind=c_double), pointer :: tmp_matrix(:, :, :)

      allocate(tmp_matrix(2, 3, 4))

      tmp_matrix=1.0_8
      tmp_matrix(2, 1, 1)=2.0_8
      tmp_matrix(1, 2, 1)=3.0_8
      tmp_matrix(2, 3, 4)=4.0_8

      matrix=c_loc(tmp_matrix)

   end subroutine testMatrix

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

   subroutine c_endMadelung() bind(C, name="endMadelung")
      use MadelungModule, only : endMadelung

      call endMadelung()

   end subroutine c_endMadelung

end module CMadelungModule