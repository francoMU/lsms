/**
 *
 * higherOrderMadelung.hpp
 *
 * Created by Franco Moitzi on 8/3/21.
 *
 */
#ifndef MUST_HIGHERORDERMADELUNG_HPP
#define MUST_HIGHERORDERMADELUNG_HPP

#include <vector>
#include <algorithm>

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

namespace lsms {

    class HigherOrderMadelung {

    public:

        /*
         * Reduced madelung matrix
         */

        // G^L_0 = G^{l,m}_0 = {1,1}
        Matrix<Complex> G_11{};

        // G^L_0 = G^{l,m}_0 = {1,-1}
        Matrix<Complex> G_1m1{};

        // G^L_0 = G^{l,m}_0  = {1,0}
        Matrix<Complex> G_10{};

        /*
         * Prefactors for Madelung matrix
         *
         * M^{L,L'}_{n,n'} = f(L,L') G^{l+l',m-m'}_{n,n'}
         *
         */

        // f(L,L_prime) = f({1,0},{0,0})
        Real pre_factor_10_00{0.0};

        // f(L,L_prime) = f({1,1},{0,0})
        Real pre_factor_11_00{0.0};

        /*
         * Lattice scaling factor
         */

        Real lattice_scale{0.0};


        explicit HigherOrderMadelung(int num_local_atoms, int num_atoms);


    };

}


#endif //MUST_HIGHERORDERMADELUNG_HPP
