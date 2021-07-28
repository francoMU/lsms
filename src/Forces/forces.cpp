
#include "forces.hpp"

#include <cstdio>
#include <string>
#include <cmath>

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"

#include "utils.hpp"
#include "forces_params.hpp"


#ifdef DEBUG

#include <icecream.hpp>

#endif

void lsms::calculateForces(const LSMSCommunication &comm,
                           const LSMSSystemParameters &lsms,
                           LocalTypeInfo &local,
                           CrystalParameters &crystal,
                           ForceParameters parameter) {

    std::vector<Real> distance(3, 0.0);

    double radius_factor = 0.0;

    // Clear forces
    for (auto i = 0; i < local.num_local; i++) {
        std::fill(local.atom[i].force.begin(),
                  local.atom[i].force.end(), 0.0);
        radius_factor += lsms::pow(local.atom[i].rmt, 3);
    }


    constexpr double radius = 100;
    constexpr double radius_sq = radius * radius;

    auto m = parameter.i_x;
    auto n = parameter.i_y;
    auto u = parameter.i_z;

#ifdef DEBUG
    std::printf("Repeats: %d %d %d\n", m, n, u);

    double bravais_x[] = {
            crystal.bravais(0, 0),
            crystal.bravais(0, 1),
            crystal.bravais(0, 2),
    };

    double bravais_y[] = {
            crystal.bravais(1, 0),
            crystal.bravais(1, 1),
            crystal.bravais(1, 2),
    };

    double bravais_z[] = {
            crystal.bravais(2, 0),
            crystal.bravais(2, 1),
            crystal.bravais(2, 2),
    };


    std::printf("distance x: %f %f\n", norm(bravais_x), norm(bravais_x)*m);
    std::printf("distance y: %f %f\n", norm(bravais_y), norm(bravais_y)*n);
    std::printf("distance z: %f %f\n", norm(bravais_z), norm(bravais_z)*u);

#endif

    // Contribution from periodic cells
    for (auto i = 0; i < local.num_local; i++) {

        auto &atom_i = local.atom[i];

        for (auto j = 0; j < local.num_local; j++) {

            auto &atom_j = local.atom[j];

            for (int i_m = -m; i_m < m; i_m++) {
                for (int i_n = -n; i_n < n; i_n++) {
                    for (int i_u = -u; i_u < u; i_u++) {

                        /*
                         * Real-space sum around certain radius
                         */
                        distance[0] = crystal.position(0, i)
                                      - crystal.position(0, j)
                                      - i_m * crystal.bravais(0, 0)
                                      - i_n * crystal.bravais(0, 1)
                                      - i_u * crystal.bravais(0, 2);

                        distance[1] = crystal.position(1, i)
                                      - crystal.position(1, j)
                                      - i_m * crystal.bravais(1, 0)
                                      - i_n * crystal.bravais(1, 1)
                                      - i_u * crystal.bravais(1, 2);

                        distance[2] = crystal.position(2, i)
                                      - crystal.position(2, j)
                                      - i_m * crystal.bravais(2, 0)
                                      - i_n * crystal.bravais(2, 1)
                                      - i_u * crystal.bravais(2, 2);


                        auto distance_l2 = norm(distance);
                        auto distance_sq = distance_l2 * distance_l2;

                        if (distance_sq > radius_sq) {
                            continue;
                        }

                        auto distance_l2_3 = distance_l2 * distance_l2 * distance_l2;

                        auto charge_j = atom_j.ztotss   // nucleus charge
                                        - atom_j.qtotws // electron charge
                                        + atom_j.rhoInt * atom_j.omegaWS;

                        auto charge_i = atom_i.ztotss   // nucleus charge
                                        - atom_i.qtotws // electron charge
                                        + atom_i.rhoInt * atom_i.omegaWS;


#ifdef DEBUG
                        std::printf("rhoInt %f\n", atom_j.rhoInt);
                                      // + atom_j.qInt; // interstitial contribution
                                        //* lsms::pow(atom_j.rmt, 3) / radius_factor;

                        auto charge_1 = atom_j.rhoInt * atom_j.omegaMT;
                        auto charge_2 = atom_j.rhoInt * lsms.volumeTotal;
                        auto charge_3 = atom_j.qInt * atom_j.omegaMT / lsms.volumeInterstitial;

                        std::printf("rhoMT %f %f %f\n", charge_1, charge_2, charge_3);
                        std::printf("V %f %f %f\n",
                                    atom_i.omegaMT,
                                    lsms.volumeInterstitial,
                                    lsms.volumeTotal);

                        std::printf("Z %f\n", atom_j.ztotss);
                        std::printf("Val %f\n", atom_j.qvalmt);
                        std::printf("Core %f\n",atom_j.qcpsc_mt);
                        std::printf("Int %f\n", atom_j.qInt);
                        std::printf("Charge %f\n", charge);
#endif


                        if ((i_m != 0 && i_n != 0 && i_u != 0) || i != j) {
                            /*
                             * Forces in Rydberg per angstroem
                             */
                            atom_i.force[0] += 2.0 * charge_i * charge_j *
                                               distance[0] / distance_l2_3;
                            atom_i.force[1] += 2.0 * charge_i * charge_j *
                                               distance[1] / distance_l2_3;
                            atom_i.force[2] += 2.0 * charge_i * charge_j *
                                               distance[2] / distance_l2_3;

                        }


                    }
                }
            }
        }
    }


}


void lsms::displayForces(const LocalTypeInfo &local, const CrystalParameters &crystal) {


    std::string filler_string(80, '*');

    std::printf("%s\n", filler_string.c_str());


    for (int i = 0; i < local.num_local; i++) {
        auto &atom = local.atom[i];

        std::printf("** Atom: %4d\n", i);
        std::printf("** Force x: %24.16f\n", atom.force[0]);
        std::printf("** Force y: %24.16f\n", atom.force[1]);
        std::printf("** Force z: %24.16f\n", atom.force[2]);
    }

    std::printf("%s\n", filler_string.c_str());


}

void lsms::displayCharges(const LSMSCommunication &comm,
                          const LSMSSystemParameters &lsms,
                          const LocalTypeInfo &local,
                          const CrystalParameters &crystal) {


    std::string filler_string(80, '*');

    std::printf("%s\n", filler_string.c_str());

    double total_muffin_tin_charge = 0.0;
    double total_nuclei_charge = 0.0;
    double total_semi_core_charge = 0.0;
    double total_core_charge = 0.0;
    double total_valence_charge = 0.0;

    double total_muffin_tin_core_charge = 0.0;

    for (int i = 0; i < local.num_local; i++) {
        auto &atom = local.atom[i];


        std::printf("** Atom: %4d\n", i);
        std::printf("** Charge MT:    %24.16f\n", atom.qvalmt);
        std::printf("** Charge WS:    %24.16f\n", atom.qvalws);
        std::printf("** Charge Z:     %24.16f\n", atom.ztotss);
        std::printf("** C+S MT:       %24.16f\n", atom.mcpsc_mt);
        std::printf("** C+S WS:       %24.16f\n", atom.mcpsc_ws);
        std::printf("** C+S MT:       %24.16f\n", atom.qcpsc_mt);
        std::printf("** C+S WS:       %24.16f\n", atom.qcpsc_ws);
        std::printf("** MT Radius:    %24.16f\n", atom.rmt);

        total_muffin_tin_charge += atom.qvalmt;

        total_muffin_tin_core_charge += atom.mcpsc_mt;

        total_nuclei_charge += atom.ztotss;
        total_core_charge += atom.zcorss;
        total_semi_core_charge += atom.zsemss;
        total_valence_charge += atom.zvalss;

    }

    auto total_charge = total_muffin_tin_charge + local.atom[0].qInt;

    std::printf("** Total\n");
    std::printf("** Total MT Charge:       %24.16f\n", total_muffin_tin_charge);
    std::printf("** Total Charge:          %24.16f\n", total_charge);
    std::printf("** Total MT Core Charge:  %24.16f\n", total_muffin_tin_core_charge);
    std::printf("** Total:                 %24.16f\n", total_muffin_tin_core_charge
                                                       + total_charge);

    std::printf("**\n");
    std::printf("** Total Nuclei Charge:   %24.16f\n", total_nuclei_charge);
    std::printf("** Total Core Charge:     %24.16f\n", total_core_charge);
    std::printf("** Total Semicore Charge: %24.16f\n", total_semi_core_charge);
    std::printf("** Total Valence Charge:  %24.16f\n", total_valence_charge);
    std::printf("**\n");
    std::printf("%s\n", filler_string.c_str());

}
