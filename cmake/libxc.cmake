

if (NOT DEFINED Libxc_LIBRARIES)
    set(Libxc_LIBRARIES
            ${PROJECT_BINARY_DIR}/external/Libxc/lib/${CMAKE_STATIC_LIBRARY_PREFIX}Libxc${CMAKE_STATIC_LIBRARY_SUFFIX})

    set(Libxc_INCLUDE_DIR
            ${PROJECT_BINARY_DIR}/external/Libxc/include
            )
endif ()

set(Libxc_LIBRARIES ${Libxc_LIBRARIES}
        CACHE FILEPATH "Libxc library" FORCE)
set(Libxc_INCLUDE_DIR ${Libxc_INCLUDE_DIR}
        CACHE FILEPATH "Libxc include dirs" FORCE)


if (EXISTS ${Libxc_LIBRARIES} AND EXISTS ${Libxc_INCLUDE_DIR})
    set(Libxc_FOUND true)
    message(STATUS "Libxc was found")
    message(STATUS "Libxc library: " ${Libxc_LIBRARIES})
    message(STATUS "Libxc include: " ${Libxc_INCLUDE_DIR})
endif ()

if (NOT Libxc_FOUND)

    message(STATUS "LIBXC")

    find_program(AUTORECONF_EXECUTABLE
            NAMES autoreconf
            DOC "Autoreconf" REQUIRED)

    find_program(AUTOCONF_EXECUTABLE
            NAMES autoconf
            DOC "Autoconf" REQUIRED)

    find_program(AUTOMAKE_EXECUTABLE
            NAMES automake
            DOC "Automake" REQUIRED)

    find_program(MAKE_EXECUTABLE
            NAMES gmake make
            NAMES_PER_DIR
            DOC "GNU Make")

    #
    # dos2unix
    #

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Autotools
            REQUIRED_VARS AUTOCONF_EXECUTABLE AUTOMAKE_EXECUTABLE MAKE_EXECUTABLE)


    file(COPY ${PROJECT_SOURCE_DIR}/external/libxc-5.1.6
            DESTINATION ${PROJECT_BINARY_DIR}/external)

    set(_src ${PROJECT_BINARY_DIR}/external/libxc-5.1.6)
    get_filename_component(_src "${_src}" REALPATH)

    set(_install ${PROJECT_BINARY_DIR}/external/libxc)
    file(MAKE_DIRECTORY ${_install})
    get_filename_component(_install "${_install}" REALPATH)
    include(ExternalProject)

    ExternalProject_Add(libxc
            SOURCE_DIR ${_src}
            BUILD_IN_SOURCE true
            #CONFIGURE_COMMAND find . -type f -print0 | xargs -0 dos2unix
            CONFIGURE_COMMAND ${AUTORECONF_EXECUTABLE} -i
            COMMAND ./configure --prefix=${_install} CC=${CMAKE_C_COMPILER}
            BUILD_COMMAND ${MAKE_EXECUTABLE}
            INSTALL_COMMAND ${MAKE_EXECUTABLE} install
            )

endif()


if (NOT TARGET libxc::libxc)
    add_library(libxc::libxc INTERFACE IMPORTED GLOBAL)
    target_include_directories(libxc::libxc INTERFACE ${Libxc_INCLUDE_DIR})
    target_link_libraries(libxc::libxc INTERFACE ${Libxc_LIBRARIES})
endif ()

if (NOT libxc_FOUND)
    add_dependencies(libxc::libxc libxc)
endif ()