


target_sources(lsmscore
        PUBLIC
        StructureIO.hpp
        POSCARStructureIO.hpp
        POSCARStructureIO.cpp
        POSCARStructureType.hpp
        orderedPrint.hpp
        orderedPrint.cpp
        mpi_chain.hpp
        chain_stream.hpp
        chain_stream.cpp
        io_utils.cpp
        io_utils.hpp
        )

target_include_directories(lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )