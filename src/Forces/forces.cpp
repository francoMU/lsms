
#include "forces.hpp"

#include <cstdio>
#include <string>
#include <cmath>
#include <numeric>

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"

#include "utils.hpp"
#include "forces_params.hpp"


#ifdef DEBUG

#include <icecream.hpp>

#endif

void lsms::gatherStructure(const LSMSCommunication &comm,
                           const LSMSSystemParameters &lsms,
                           const CrystalParameters &crystal,
                           std::vector<int> &rcounts,
                           std::vector<int> &displs) {

    std::fill(rcounts.begin(),
              rcounts.end(), 0.0);

    std::fill(displs.begin(),
              displs.end(), 0.0);

    // number of elements in the message to receive from each process
    rcounts.resize(comm.size);

    for (int i = 0; i < lsms.num_atoms; i++) {
        auto index = crystal.types[i].node;
        rcounts[index] = crystal.types[i].local_id + 1;
    }

    // displacement vector
    displs.resize(comm.size);

    for (int i = 1; i < comm.size; i++) {
        displs[i] = std::accumulate(rcounts.begin(),
                                    rcounts.begin() + i, 0);
    }

}

void lsms::normalizeForces(const LSMSCommunication &comm,
                           const LSMSSystemParameters &lsms,
                           LocalTypeInfo &local,
                           const CrystalParameters &crystal) {

    //  Number of local atoms
    auto num_local_atoms = local.num_local;
    // Gives the total number of atoms
    auto num_atoms = lsms.num_atoms;


    std::array<Real, 3> local_force_sum{{0.0,
                                                0.0,
                                                0.0}};

    std::array<Real, 3> global_forces_sum{0.0};

    // Local force sum
    for (int i_local = 0; i_local < num_local_atoms; i_local++) {
        for (int i_vec = 0; i_vec < 3; i_vec++) {
            local_force_sum[i_vec] += local.atom[i_local].force[i_vec];
        }
    }

    MPI_Allreduce(local_force_sum.data(),
                  global_forces_sum.data(),
                  3, MPI_DOUBLE, MPI_SUM, comm.comm);

    for (int i_local = 0; i_local < num_local_atoms; i_local++) {
        for (int i_vec = 0; i_vec < 3; i_vec++) {
            local.atom[i_local].force[i_vec] -= global_forces_sum[i_vec] / num_atoms;
        }
    }


#ifdef DEBUG
    std::cout << comm.rank << " " << local_force_sum[0] << std::endl;
    std::cout << comm.rank << " " << local_force_sum[1] << std::endl;
    std::cout << comm.rank << " " << local_force_sum[2] << std::endl;
    std::cout << comm.rank << " " << global_forces_sum[0] << std::endl;
    std::cout << comm.rank << " " << global_forces_sum[1] << std::endl;
    std::cout << comm.rank << " " << global_forces_sum[2] << std::endl;
#endif

}


void lsms::gatherForces(const LSMSCommunication &comm,
                        const LSMSSystemParameters &lsms,
                        const LocalTypeInfo &local,
                        const CrystalParameters &crystal,
                        std::vector<std::array<Real, 3>> &global_forces,
                        int root) {

    //  Number of local atoms
    auto num_local_atoms = local.num_local;
    // Gives the total number of atoms
    auto num_atoms = lsms.num_atoms;

    MPI_Datatype mpi_force_type;
    auto status = MPI_Type_contiguous(3, MPI_DOUBLE, &mpi_force_type);
    MPI_Type_commit(&mpi_force_type);

    /*
     * Create array of local forces that is contiguous
     */
    std::vector<std::array<double, 3>> local_forces(num_local_atoms);
    global_forces.resize(num_atoms);

    for (int i_local = 0; i_local < num_local_atoms; i_local++) {
        //std::copy(local.atom[i_local].force.begin(),
        //          local.atom[i_local].force.end(),
        //          local_forces[i_local]);
        local_forces[i_local][0] = local.atom[i_local].force[0];
        local_forces[i_local][1] = local.atom[i_local].force[1];
        local_forces[i_local][2] = local.atom[i_local].force[2];
    }

    // number of elements in the message to receive from each process
    std::vector<int> rcounts(comm.size);
    std::vector<int> displs(comm.size, 0);

    gatherStructure(comm,
                    lsms,
                    crystal,
                    rcounts,
                    displs);

    status = MPI_Gatherv(local_forces.data(),
                         num_local_atoms,
                         mpi_force_type,
                         global_forces.data(),
                         rcounts.data(),
                         displs.data(),
                         mpi_force_type,
                         root,
                         comm.comm);

#ifdef DEBUG
    if (comm.rank == root) {
        for (auto &force : global_forces) {
            std::cout << " ***** " << std::endl;
            std::cout << force[0] << std::endl;
            std::cout << force[1] << std::endl;
            std::cout << force[2] << std::endl;
        }
    }
#endif

    status = MPI_Type_free(&mpi_force_type);
}


void lsms::calculateForces(const LSMSCommunication &comm,
                           const LSMSSystemParameters &lsms,
                           LocalTypeInfo &local,
                           CrystalParameters &crystal,
                           lsms::HigherOrderMadelung &madelung,
                           lsms::ForceParameters parameters) {

    //  Number of local atoms
    auto num_local_atoms = local.num_local;

    // Gives the total number of atoms
    auto num_atoms = lsms.num_atoms;
    auto latt_scale = madelung.lattice_scale;

    auto prefactor_11 = madelung.pre_factor_11_00;
    auto prefactor_10 = madelung.pre_factor_10_00;

    // Clear forces
    for (auto i = 0; i < local.num_local; i++) {
        std::fill(local.atom[i].force.begin(),
                  local.atom[i].force.end(), 0.0);
    }

    std::vector<Real> local_charges(num_local_atoms, 0.0);
    std::vector<Real> global_charges(num_atoms, 0.0);

    // Number of local atoms
    for (int i_local = 0; i_local < num_local_atoms; i_local++) {

        auto &atom_i = local.atom[i_local];

        // Charge multipole memonent i
        auto charge_i = atom_i.ztotss   // nucleus charge
                        - atom_i.qtotmt // electron charge
                        + atom_i.rhoInt * atom_i.omegaMT;

        local_charges[i_local] = charge_i;

        auto i_global = local.global_id[i_local];
        global_charges[i_global] = charge_i;

    }

    // number of elements in the message to receive from each process
    std::vector<int> rcounts(comm.size);

    for (int i = 0; i < num_atoms; i++) {
        auto index = crystal.types[i].node;
        rcounts[index] = crystal.types[i].local_id + 1;
    }

    // displacement vector
    std::vector<int> displs(comm.size, 0);

    for (int i = 1; i < comm.size; i++) {
        displs[i] = std::accumulate(rcounts.begin(),
                                    rcounts.begin() + i, 0);
    }

    MPI_Allgatherv(local_charges.data(),
                   num_local_atoms,
                   MPI_DOUBLE,
                   global_charges.data(),
                   rcounts.data(),
                   displs.data(),
                   MPI_DOUBLE,
                   comm.comm);


#ifdef DEBUG

    std::cout << " Lattice scale " << std::endl;
    std::cout << latt_scale << std::endl;

    std::cout << " Prefactors " << std::endl;
    std::cout << lsms.rank << " " << prefactor_11 << std::endl;
    std::cout << lsms.rank << " " << prefactor_10 << std::endl;

    std::cout << " Local charges " << std::endl;

    for (int i_local = 0; i_local < num_local_atoms; i_local++) {
        auto charge_i = local_charges[i_local];
        std::cout << lsms.rank << " " << i_local << " " << charge_i << std::endl;
    }

    std::cout << " Global charges " << std::endl;

    for (int j_global = 0; j_global < num_atoms; j_global++) {
        auto charge_j = global_charges[j_global];
        std::cout << lsms.rank << " " << j_global << " " << charge_j << std::endl;
    }

#endif


    /*
     *
     */
    for (int i_local = 0; i_local < num_local_atoms; i_local++) {

        auto &atom_i = local.atom[i_local];
        auto charge_i = atom_i.ztotss;

        for (int j_global = 0; j_global < num_atoms; j_global++) {
            auto &atom_j = local.atom[j_global];

            Complex G_11 = madelung.G_11(j_global, i_local);
            Complex G_10 = madelung.G_10(j_global, i_local);

#ifdef DEBUG
            std::cout << "G_11: " << G_11 << std::endl;
            std::cout << "G_10: " << G_10 << std::endl;
#endif

            auto charge_j = global_charges[j_global];

            charge_j *= std::sqrt(4.0 * M_PI);

            /*
             * F^(i)_x = - q^(i)_{eff} * sqrt(3 / (2 * M_PI)) *
             *              Re{ M^(ij)_{11,00} } * q^(j)_{eff} * sqrt(3 / (2 * M_PI)) / a_0
             */
            atom_i.force[0] += -std::sqrt(3.0 / (2.0 * M_PI))
                               * std::real(G_11 * prefactor_11) * charge_i / latt_scale
                               * charge_j;

            /*
             * F^(i)_y = + q^(i)_{eff} * sqrt(3 / (2 * M_PI)) *
             *              Im{ M^(ij)_{11,00} } * q^(j)_{eff} * sqrt(3 / (2 * M_PI)) / a_0
             */
            atom_i.force[1] += std::sqrt(3.0 / (2.0 * M_PI))
                               * std::imag(G_11 * prefactor_11) * charge_i / latt_scale
                               * charge_j;

            /*
             * F^(i)_z =  q^(i)_{eff} * sqrt(3 / (2 * M_PI)) *
             *               Re{ M^(ij)_{10,00} } * q^(j)_{eff} * sqrt(3 / (4 * M_PI)) / a_0
             */
            atom_i.force[2] += std::sqrt(3.0 / (4.0 * M_PI))
                               * std::real(G_10 * prefactor_10) * charge_i / latt_scale
                               * charge_j;

        }

    }

}


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


    std::printf("distance x: %f %f\n", norm(bravais_x), norm(bravais_x) * m);
    std::printf("distance y: %f %f\n", norm(bravais_y), norm(bravais_y) * n);
    std::printf("distance z: %f %f\n", norm(bravais_z), norm(bravais_z) * u);

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


void lsms::displayForces(const LSMSCommunication &comm,
                         const LSMSSystemParameters &lsms,
                         const LocalTypeInfo &local,
                         const CrystalParameters &crystal,
                         int rank) {

    std::vector<std::array<Real, 3>> global_forces(crystal.num_atoms);

    gatherForces(comm, lsms, local, crystal, global_forces);

    if (comm.rank == 0) {

        //std::string filler_string(80, '*');
        //std::printf("%s\n", filler_string.c_str());

        std::printf("\nForces:\n");

        std::array<Real, 3> global_forces_sum{0.0};

        for (int i_global = 0; i_global < crystal.num_atoms; i_global++) {

            auto &force = global_forces[i_global];

            for (int i_vec = 0; i_vec < 3; i_vec++) {
                global_forces_sum[i_vec] += force[i_vec];
            }

            std::printf("%6d: %24.16f %24.16f %24.16f\n",
                        i_global,
                        force[0],
                        force[1],
                        force[2]);
        }

        std::printf("%7s %24.16f %24.16f %24.16f\n",
                    "Global:",
                    global_forces_sum[0],
                    global_forces_sum[1],
                    global_forces_sum[2]);

        //std::printf("%s\n", filler_string.c_str());

    }


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
