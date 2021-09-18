#
# Created by Franco P. Moitzi
#

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
message(STATUS "CMake ${CMAKE_VERSION} with max. ${Ncpu} threads")

project(lsms LANGUAGES Fortran C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Default build type is RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE CACHE BOOL
            "Choose the type of build, options are: Debug Release"
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Build external libraries
message(STATUS "Build external libraries")

# Build LUA in this project
include(lua)


#
# Enable Testing with CTest
#
option(BUILD_TESTING "Build the testing tree." ON)
include(CTest)
enable_testing()
message(STATUS "Testing is " ${BUILD_TESTING})
add_subdirectory(unittests)


if (BUILD_LSMS)

    #
    # Add main source subdirectory
    #
    add_subdirectory(src)

    add_subdirectory(include)

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

    # Lua
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
    if (NOT DEFINED Libxc_DIR)
        message(FATAL_ERROR "Libxc_DIR is not defined")
    endif ()
    find_library(
            LIBXC xc
            PATHS ${Libxc_DIR}/lib
            NO_DEFAULT_PATH REQUIRED)

    message(STATUS "LibXC directory: " ${Libxc_DIR})
    message(STATUS "LibXC library:   " ${LIBXC})
    message(STATUS "LibXC headers:   " ${Libxc_DIR}/include)

    target_compile_definitions(lsmscore PUBLIC "USE_LIBXC")
    target_include_directories(lsmscore PUBLIC ${Libxc_DIR}/include)
    target_link_directories(lsmscore PUBLIC ${Libxc_DIR}/lib)
    target_link_libraries(lsmscore PUBLIC ${LIBXC})

    # DL
    target_link_libraries(lsmscore PUBLIC ${CMAKE_DL_LIBS})

    target_compile_definitions(lsmscore PUBLIC HAS_BACKTRACE)

    #
    # Enable Testing with CTest
    #
    option(BUILD_TESTING "Build the testing tree." ON)
    include(CTest)
    enable_testing()
    message(STATUS "Testing is " ${BUILD_TESTING})
    add_subdirectory(tests)

endif ()