! Created by fmoitzi on 8/29/21.

program test_integer_factors
   use IntegerFactorsModule
   use integer_factors_mod, only: IntegerFactors

   implicit none

   type(IntegerFactors) :: factors

   call factors % init(2)

   call initIntegerFactors(2)

   print *, all(lofk == factors % lofk)

   call endIntegerFactors()

end program test_integer_factors