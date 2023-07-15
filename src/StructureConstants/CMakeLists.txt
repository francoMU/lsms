

add_library(structureconsts SHARED)

target_sources(structureconsts PUBLIC
        Dcoeff.cpp
        Dcoeff.h
        )

target_include_directories(structureconsts PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})