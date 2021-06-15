

if (BUILD_TESTING)

    find_package(MuSTTests)

    message(STATUS "Generating tests...")

    add_subdirectory(Core)

endif()