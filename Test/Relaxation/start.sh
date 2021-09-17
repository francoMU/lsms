#!/bin/bash

LSMS=/home/fmoitzi/Development/CLionProjects/lsms/cmake-build-build_release_lsms/bin/lsms_relaxation
LSMS_LEGACY=/home/fmoitzi/Development/CLionProjects/lsms/cmake-build-build_release_lsms/bin/lsms


CWD=$(pwd)

cd $CWD || exit

echo $(date)

cd /home/fmoitzi/Development/CLionProjects/lsms/Test/Relaxation/000-check_conformity || exit
mpirun -n 4 $LSMS i_lsms.lua 2>&1 | tee output.out
mpirun -n 4 $LSMS_LEGACY i_lsms.lua 2>&1 | tee output_legacy.out

echo $(date)

cd $CWD || exit