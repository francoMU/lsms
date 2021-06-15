

cd LSMS-al_1_z3_1
mpirun -n 4 /home/fmoitzi/CLionProjects/lsms/cmake-build-release-gcc-9_mpi/bin/lsms < i_lsms | tee  outfile
cd ..

cd LSMS-al_1_z3_0
mpirun -n 4 /home/fmoitzi/CLionProjects/lsms/cmake-build-release-gcc-9_mpi/bin/lsms < i_lsms | tee  outfile
cd ..