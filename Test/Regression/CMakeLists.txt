
if (NOT USE_ACCELERATOR_CUDA_C AND NOT USE_ACCELERATOR_HIP)

    message(STATUS "Regression tests are used")

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/al_fcc
            DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/fe_bcc_libxc
            DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    find_program (BASH_PROGRAM bash)

    # 1. Regression test: Al fcc
    add_test(NAME regtest_al_fcc
            COMMAND
            ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=1 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:lsms_main> i_lsms.lua &&
            ${BASH_PROGRAM} ${CMAKE_CURRENT_BINARY_DIR}/al_fcc/test.sh
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/al_fcc
            )

    # 2. Regression test: Fe bcc (LIBXC)
    add_test(NAME regtest_fe_bcc_libxc
            COMMAND
            ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=1 ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:lsms_main> i_lsms.lua &&
            ${BASH_PROGRAM} ${CMAKE_CURRENT_BINARY_DIR}/fe_bcc_libxc/test.sh
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fe_bcc_libxc
            )

else()

    message(STATUS "Regression tests are not used")

endif ()