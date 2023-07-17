

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