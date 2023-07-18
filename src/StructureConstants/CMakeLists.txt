

add_library(structureconsts SHARED)

add_executable(sc_main main.cpp)

target_link_libraries(sc_main PUBLIC structureconsts)


target_sources(structureconsts PUBLIC
        Structure.cpp
        Structure.hpp
        NDArray.hpp
        #        Dcoeff.cpp
        #        Dcoeff.h
        )

target_include_directories(structureconsts PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


#
# Find Eigen3
#

find_package(Eigen3 REQUIRED)

target_link_libraries(structureconsts PUBLIC Eigen3::Eigen)



include(FetchContent)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

add_executable(test_ndarray test_ndarray.cpp)
target_include_directories(test_ndarray PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_ndarray PUBLIC structureconsts)

gtest_discover_tests(test_ndarray)