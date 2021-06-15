
#ifndef MUST_ENERGYSTRUCT_HPP
#define MUST_ENERGYSTRUCT_HPP


#include <Real.hpp>

namespace lsms {

    class EnergyStruct {

    public:

        /**
         * @brief struct representing total energy contribution
         *
         * @param eigenvalueSum
         * @param one_electron
         * @param kohn_sham
         * @param madelung
         * @param core_interaction
         * @param hartree
         * @param exchange_correlation
         */
        EnergyStruct(Real one_electron,
                     Real kohn_sham,
                     Real madelung,
                     Real core_interaction,
                     Real hartree,
                     Real exchange_correlation,
                     Real zero_point_energy) :
                one_electron{one_electron},
                kohn_sham{kohn_sham},
                madelung{madelung},
                core_interaction{core_interaction},
                hartree{hartree},
                exchange_correlation{exchange_correlation},
                zero_point_energy{zero_point_energy} {};

        EnergyStruct() = default;

        Real one_electron{0.0};
        Real kohn_sham{0.0};
        Real madelung{0.0};
        Real core_interaction{0.0};
        Real hartree{0.0};
        Real exchange_correlation{0.0};
        Real zero_point_energy{0.0};

        /**
         * @brief Get Kohn-Sham kinetic energy
         *
         * @return kohn sham kinetic energy
         */
        [[nodiscard]] Real kinetic_energy() const;

        /**
         * @brief All electrostatic interactions
         *
         * @return
         */
        [[nodiscard]] Real coloumb() const;

        /**
         * @brief Total energy
         *
         * @return
         */
        [[nodiscard]] Real total_energy() const;


    };

}


#endif //MUST_ENERGYSTRUCT_HPP
