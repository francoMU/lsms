
#
#
#
add_executable(test_relaxation)

target_sources(test_relaxation
        PUBLIC
        relaxTest.cpp
        )

target_link_libraries(test_relaxation PUBLIC relaxation)

#
#
#
add_executable(test_poscar)

target_sources(test_poscar
        PUBLIC
        POSCARTest.cpp
        )
