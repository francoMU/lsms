
#include "energyStruct.hpp"

Real lsms::EnergyStruct::coloumb() const {
    return hartree + core_interaction + madelung;
}

Real lsms::EnergyStruct::kinetic_energy() const {
    return one_electron - kohn_sham;
}

Real lsms::EnergyStruct::total_energy() const {
    return kinetic_energy() + coloumb() + exchange_correlation;
}
