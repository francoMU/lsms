
#include "Structure.hpp"

SC::Structure::Structure(NDArray<double, 2> &lattice, NDArray<double, 2> &coords, std::vector<int> &species)
: lattice{lattice}, coords{coords}, species{species}
{

  /*
   * Perform some size checks
   *
   */





}

void SC::Structure::calculateLLLreduction() {

//  subroutine calculate_lll(lattice, lll, lll_mapping, lll_inverse, delta)
//  use linalg_mod, only: diag, lstsq, inverse
//  use utils_mod, only: option
//
//  !> Lattice (column-wise)
//  real(kind=dp), intent(in) :: lattice(3, 3)
//
//  !> Reduced lattice matrix
//  real(kind=dp), intent(out) :: lll(3, 3)
//
//  !> Mapping to get to that lattice
//  real(kind=dp), intent(out) :: lll_mapping(3, 3)
//
//  !> Inverse of the mapping
//  real(kind=dp), intent(out) :: lll_inverse(3, 3)
//
//  !> Reduction parameter
//  real(kind=dp), optional, intent(in) :: delta
//
//  !-- Auxiliar variables
//      integer :: k, i, s, q
//
//  real(kind=dp) :: a(0:2, 0:2)
//  real(kind=dp) :: b(0:2, 0:2)
//  real(kind=dp) :: u(0:2, 0:2)
//  real(kind=dp) :: m(0:2)
//  real(kind=dp) :: mu
//  real(kind=dp) :: mapping(0:2, 0:2)
//  real(kind=dp) :: uu(0:9), v(3), v_m(3)
//  real(kind=dp) :: qq(2, 2), pp(2)
//
//  real(kind=dp) :: delta_

double delta = 0.75;
//      delta_ = option(0.75_dp, delta)
//
//  uu = 0.0_dp

NDArray<double, 2> uu(3,3);

uu = 0.0;

//
//  Basis vectors needs to be in column representation
//a = lattice



//  b = 0.0_dp !< Vectors after the Gram-Schmidt process
//  u = 0.0_dp !< Gram-Schmidt coeffiencts
//  m = 0.0_dp !< These are the norm squared of each vector
//
//  b(:, 0) = a(:, 0)
//  m(0) = dot_product(b(:, 0), b(:, 0))
//
//  do i = 1, 2
//  u(i, 0:i - 1) = matmul(a(:, i), b(:, 0:i - 1)) / m(0:i - 1)
//  b(:, i) = a(:, i) - matmul(b(:, 0:i - 1), u(i, 0:i - 1))
//  m(i) = dot_product(b(:, i), b(:, i))
//  end do
//
//    !-- Checked
//  k = 2
//  mapping = 0.0_dp
//  mapping(0, 0) = 1.0_dp
//  mapping(1, 1) = 1.0_dp
//  mapping(2, 2) = 1.0_dp
//
//  do while (k <= 3)
//
//      !--- length reduction step
//  do i = k - 1, 1, -1
//
//  mu = u(k - 1, i - 1)
//
//  ! Reduce the k-th basis vector
//  if ((abs(mu)) > 0.5) then
//
//        q = nint(u(k - 1, i - 1))
//    a(:, k - 1) = a(:, k - 1) - q * a(:, i - 1)
//  mapping(:, k - 1) = mapping(:, k - 1) - q * mapping(:, i - 1)
//
//  uu = 0.0_dp
//  if (i >= 2) then
//    uu(0:(i - 2)) = u(i - 1, 0:(i - 2))
//  end if
//
//    uu(i - 1) = 1.0_dp
//
//  u(k - 1, 0:i - 1) = u(k - 1, 0:i - 1) - q * uu(0:i - 1)
//
//  end if
//
//    end do
//
//    !
//        !--- Swap step: Check the Lovasz condition.
//  !
//  if (dot_product(b(:, k - 1), b(:, k - 1)) &
//      >= (delta_ - abs(u(k - 1, k - 2)) ** 2) * dot_product(b(:, (k - 2)), b(:, (k - 2)))) then
//
//  !
//      ! Increment k if the Lovasz condition holds.
//  !
//      k = k + 1
//
//  else
//
//  !
//      ! If the Lovasz condition fails,
//  ! swap the k-th and (k-1)-th basis vector
//  !
//
//      v = a(:, k - 1)
//  a(:, k - 1) = a(:, k - 2)
//  a(:, k - 2) = v
//
//  v_m = mapping(:, k - 1)
//  mapping(:, k - 1) = mapping(:, k - 2)
//  mapping(:, k - 2) = v_m
//
//  do s = k - 1, k
//  u(s - 1, 0:(s - 2)) = matmul(a(:, s - 2), b(:, 0:(s - 2))) / m(0:(s - 2))
//  b(:, s - 1) = a(:, s - 1) - matmul(b(:, 0:(s - 2)), u(s - 1, 0:(s - 2)))
//  m(s - 1) = dot_product(b(:, s - 1), b(:, s - 1))
//  end do
//
//    if (k > 2) then
//          k = k - 1
//        else if (k == 2) then
//    pp = matmul(a(:, k), b(:, (k - 2):(k - 1)))
//  qq = transpose(diag(m((k - 2):(k - 1))))
//
//  call lstsq(qq, pp, u(k, (k - 2):(k - 1)))
//
//  end if
//
//    end if
//
//    end do
//
//    ! Just transfer it too better matrix representation
//  lll = transpose(a)
//
//  lll_mapping = transpose(mapping)
//  lll_inverse = inverse(lll_mapping)
//
//  end subroutine calculate_lll

}