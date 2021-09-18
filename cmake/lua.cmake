

if (NOT DEFINED Lua_LIBRARIES)
    set(Lua_LIBRARIES
            ${PROJECT_BINARY_DIR}/external/lua/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lua${CMAKE_STATIC_LIBRARY_SUFFIX})

    set(Lua_INCLUDE_DIR
            ${PROJECT_BINARY_DIR}/external/lua/include
            )
endif ()

set(Lua_LIBRARIES ${Lua_LIBRARIES}
        CACHE FILEPATH "Lua library" FORCE)
set(Lua_INCLUDE_DIR ${Lua_INCLUDE_DIR}
        CACHE FILEPATH "Lua include dirs" FORCE)


if (EXISTS ${Lua_LIBRARIES} AND EXISTS ${Lua_INCLUDE_DIR})
    set(Lua_FOUND true)
    message(STATUS "Lua was found")
    message(STATUS "Lua library: " ${Lua_LIBRARIES})
    message(STATUS "Lua include: " ${Lua_INCLUDE_DIR})
endif ()

if (NOT Lua_FOUND)

    find_program(MAKE_EXECUTABLE NAMES gmake make REQUIRED)

    file(COPY ${PROJECT_SOURCE_DIR}/external/lua-5.2.4
            DESTINATION ${PROJECT_BINARY_DIR}/external)

    set(_src ${PROJECT_BINARY_DIR}/external/lua-5.2.4)
    get_filename_component(_src "${_src}" REALPATH)

    set(_install ${PROJECT_BINARY_DIR}/external/lua)
    file(MAKE_DIRECTORY ${_install})
    get_filename_component(_install "${_install}" REALPATH)

    set(_include ${PROJECT_BINARY_DIR}/external/lua/include)
    file(MAKE_DIRECTORY ${_include})
    get_filename_component(_include "${_include}" REALPATH)

    include(ExternalProject)

    ExternalProject_Add(Lua
            SOURCE_DIR ${_src}
            BUILD_IN_SOURCE false
            CONFIGURE_COMMAND sed -i "/^INSTALL_TOP/c INSTALL_TOP=${_install}" ${_src}/Makefile
            BUILD_COMMAND ${MAKE_EXECUTABLE} -C ${_src}
            INSTALL_COMMAND ${MAKE_EXECUTABLE} install -C ${_src}
            )


endif ()

if (NOT TARGET Lua::Lua)
    add_library(Lua::Lua INTERFACE IMPORTED GLOBAL)
    target_include_directories(Lua::Lua INTERFACE ${Lua_INCLUDE_DIR})
    target_link_libraries(Lua::Lua INTERFACE ${Lua_LIBRARIES})
endif ()

if (NOT Lua_FOUND)
    add_dependencies(Lua::Lua Lua)
endif ()