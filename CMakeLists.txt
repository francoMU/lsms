cmake_minimum_required(VERSION 3.18 FATAL_ERROR)


project(MuST
        LANGUAGES
        Fortran
        C
        CXX)


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)


#
# Add main source subdirectory
#
add_subdirectory(src)

add_subdirectory(include)

#
#
#


# Find MPI library
find_package(MPI REQUIRED)
target_link_libraries(lsmscore PUBLIC MPI::MPI_CXX)

# Find OpenMP library
find_package(OpenMP REQUIRED)
target_link_libraries(lsmscore PUBLIC OpenMP::OpenMP_CXX)

# Find HDF5
set(HDF5_PREFER_PARALLEL false)
find_package(HDF5 REQUIRED)
target_link_libraries(lsmscore PUBLIC HDF5::HDF5)

# Lapack
if (NOT DEFINED BLA_VENDOR)
    set(BLA_VENDOR "Generic")
endif ()
find_package(LAPACK REQUIRED)
target_link_libraries(lsmscore PUBLIC LAPACK::LAPACK)

# BLAS
find_package(BLAS REQUIRED)
target_link_libraries(lsmscore PUBLIC BLAS::BLAS)

#
if (NOT DEFINED LUA_INCLUDE_DIR)
    message(FATAL_ERROR "LUA_INCLUDE_DIR is not defined")
endif ()
if (NOT DEFINED LUA_LIB_DIR)
    message(FATAL_ERROR "LUA_LIB_DIR is not defined")
endif ()
target_include_directories(lsmscore PUBLIC ${LUA_INCLUDE_DIR})
target_link_directories(lsmscore PUBLIC ${LUA_LIB_DIR})
target_link_libraries(lsmscore PUBLIC lua)

#
target_include_directories(lsmscore PUBLIC ${CMAKE_SOURCE_DIR}/include)
add_subdirectory(mjson)
target_link_libraries(lsmscore PUBLIC mjson)

# LibXC
message(STATUS "LibXC: ${Libxc_DIR}")
find_library(LIBXC xc
        PATHS ${Libxc_DIR}/lib
        NO_DEFAULT_PATH
        REQUIRED
        )

target_include_directories(lsmscore PRIVATE ${Libxc_DIR}/include)
target_link_directories(lsmscore PRIVATE ${Libxc_DIR}/lib)
target_compile_definitions(lsmscore PRIVATE USE_LIBXC)
target_link_libraries(lsmscore PUBLIC ${LIBXC_LIBRARIES})


# DL
target_link_libraries(lsmscore PUBLIC ${CMAKE_DL_LIBS})

target_compile_definitions(lsmscore PUBLIC
        HAS_BACKTRACE)

#
# Enable Testing with CTest
#
option(BUILD_TESTING "Build the testing tree." ON)
include(CTest)
enable_testing()
message(STATUS "Testing is " ${BUILD_TESTING})
add_subdirectory(tests)

include(FetchContent)
set(WITH_C_INTERFACE yes CACHE INTERNAL "")  # Forces the value
set(CMAKE_BUILD_TYPE_OLD ${CMAKE_BUILD_TYPE})
set(CMAKE_BUILD_TYPE "Debug" CACHE INTERNAL "Set for release in dftatom")


add_executable(force dev_scripts/forces/forces.f90)
target_compile_options(force PUBLIC -O3 -mtune=native -cpp)
target_compile_definitions(force PUBLIC DEBUG)