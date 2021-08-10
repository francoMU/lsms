

if (BUILD_TESTING)

    find_package(MuSTTests)

    message(STATUS "Generating tests...")

    add_subdirectory(Core)

    add_subdirectory(Multipole)

    add_subdirectory(Relaxation)

endif()