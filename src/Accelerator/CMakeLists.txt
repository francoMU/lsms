#
# Created by Franco P. Moitzi
#

target_sources(lsmscore
        PUBLIC

        inverse.cu
        inverse.h

        operations.h

        DeviceStorage.cu
        zblock_lu_cuda_c.cu


        )

target_sources(
        lsmscore
        PUBLIC
        accelerator_initialize.F
        accelerator_finalize.F
        Accelerator.cpp
        Accelerator.hpp
        fortran.c
        fortran.h
)


target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)