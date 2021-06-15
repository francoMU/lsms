module forces_mod
   use iso_fortran_env, only : dp=>real64

   implicit none

   public :: calculate_nuclei_interactions
   public :: calculate_electron_nuclei_interactions

contains

   subroutine calculate_nuclei_interactions(lattice, frac_coords, charges, forces)
      real(dp), intent(in) :: lattice(3, 3)
      real(dp), intent(in) :: frac_coords(:, :)
      real(dp), intent(in) :: charges(size(frac_coords, 2))
      real(dp), intent(out) :: forces(size(frac_coords, 1), size(frac_coords, 2))

      real(dp) :: coords(size(frac_coords, 1), size(frac_coords, 2))

      real(dp) :: coord(3), translation_vec(3)

      integer :: idx, jdx, a_idx, b_idx, c_idx

      integer :: limit=150

      forces=0.0_dp

      ! Calculate cartesian coordinates
      coords=matmul(lattice, frac_coords)

      ! Intracell contributions
      do idx=1, size(frac_coords, 2)
         do jdx=1, size(frac_coords, 2)

            if (idx /= jdx) then

               coord=coords(:, idx) - coords(:, jdx)

               forces(:, idx)=forces(:, idx) + &
                     charges(idx) * charges(jdx) * coord / norm2(coord)**3
            end if

         end do
      end do

#ifdef DEBUG
      print *, "INTRACELL FORCES"
      print *, forces
#endif

      ! Contributions of periodic cells
      do idx=1, size(frac_coords, 2)

         do a_idx=- limit, limit + 1
            do b_idx=- limit, limit + 1
               do c_idx=- limit, limit + 1

                  ! Periodic images of the same site
                  if (a_idx /= 0 .and. b_idx /= 0 .and. c_idx /= 0) then

                     translation_vec=matmul(lattice, (/a_idx, b_idx, c_idx /)) &
                           - coords(:, idx)

                     forces(:, idx)=forces(:, idx) + &
                           charges(idx) * charges(idx) &
                                 * translation_vec / norm2(translation_vec)**3

                  end if

                  ! With periodic images of neighboring cells
                  do jdx=1, size(frac_coords, 2)

                     if (idx /= jdx) then

                        translation_vec=matmul(lattice, (/a_idx, b_idx, c_idx /)) &
                              - coords(:, jdx)

                        forces(:, idx)=forces(:, jdx) + &
                              charges(idx) * charges(jdx) &
                                    * translation_vec / norm2(translation_vec)**3
                     end if

                  end do

               end do
            end do
         end do

      end do

#ifdef DEBUG
      print *, "FORCES from periodic images"
      print *, forces
#endif

   end subroutine calculate_nuclei_interactions

   subroutine calculate_electron_nuclei_interactions(lattice, frac_coords, charge_moments, forces)
      real(dp), intent(in) :: lattice(3, 3)
      real(dp), intent(in) :: frac_coords(:, :)
      real(dp), intent(in) :: charge_moments(:, :)
      real(dp), intent(out) :: forces(size(frac_coords, 1), size(frac_coords, 2))

      real(dp) :: coords(size(frac_coords, 1), size(frac_coords, 2))

      real(dp) :: coord(3), translation_vec(3)

      integer :: idx, jdx, a_idx, b_idx, c_idx

      integer :: limit=200

      forces=0.0_dp

      ! Calculate cartesian coordinates
      coords=matmul(lattice, frac_coords)

      ! Intracell contributions
      do idx=1, size(frac_coords, 2)
         do jdx=1, size(frac_coords, 2)

            if (idx /= jdx) then

               coord=coords(:, idx) - coords(:, jdx)

               forces(:, idx)=forces(:, idx) + &
                     charges(idx) * charges(jdx) * coord / norm2(coord)**3
            end if

         end do
      end do

#ifdef DEBUG
      print *, "INTRACELL FORCES"
      print *, forces
#endif

      ! Contributions of periodic cells
      do idx=1, size(frac_coords, 2)

         do a_idx=- limit, limit + 1
            do b_idx=- limit, limit + 1
               do c_idx=- limit, limit + 1

                  ! Periodic images of the same site
                  if (a_idx /= 0 .and. b_idx /= 0 .and. c_idx /= 0) then

                     translation_vec=matmul(lattice, (/a_idx, b_idx, c_idx /)) &
                           - coords(:, idx)

                     forces(:, idx)=forces(:, idx) + &
                           charges(idx) * charges(idx) &
                                 * translation_vec / norm2(translation_vec)**3

                  end if

                  ! With periodic images of neighboring cells
                  do jdx=1, size(frac_coords, 2)

                     if (idx /= jdx) then

                        translation_vec=matmul(lattice, (/a_idx, b_idx, c_idx /)) &
                              - coords(:, jdx)

                        forces(:, idx)=forces(:, jdx) + &
                              charges(idx) * charges(jdx) &
                                    * translation_vec / norm2(translation_vec)**3
                     end if

                  end do

               end do
            end do
         end do

      end do

#ifdef DEBUG
      print *, "FORCES from periodic images"
      print *, forces
#endif

   end subroutine calculate_electron_nuclei_interactions

end module forces_mod

program forces_program
   use iso_fortran_env, only : dp=>real64
   use forces_mod
   implicit none

   real(dp) :: lattice(3, 3)
   real(dp) :: frac_coords(3, 2)
   real(dp) :: charges(2)
   real(dp) :: charge_moms(2)

   real(dp) :: forces_nn(3, 2)
   real(dp) :: forces_en(3, 2)

   lattice=reshape((/3.0, 0.0, 0.0, &
         0.0, 3.0, 0.0, &
         0.0, 0.0, 3.0/), (/3, 3/))

   charges=(/23.0, 23.0/)

   charge_mom(//)

   frac_coords=reshape((/0.0, 0.0, 0.0, &
         0.5, 0.5, 0.5/), (/3, 2/))

   call calculate_nuclei_interactions(lattice, frac_coords, charges, forces_nn)

   call calculate_electron_nuclei_interactions(lattice, frac_coords, charges_mom, forces_en)



end program