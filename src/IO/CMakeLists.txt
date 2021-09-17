


target_sources(lsmscore
        PUBLIC
        StructureIO.hpp
        POSCARStructureIO.hpp
        POSCARStructureIO.cpp
        POSCARStructureType.hpp
        orderedPrint.hpp
        )

target_include_directories(lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )