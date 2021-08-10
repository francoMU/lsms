
target_sources(lsmscore
        PUBLIC
        SystemParameters.cpp
        read_input.h
        read_input.cpp
        PotentialIO.cpp
        buildLIZandCommLists.cpp
        buildLIZandCommLists.hpp
        energyContourIntegration.cpp solveSingleScatterers.cpp
        calculateDensities.cpp calculateChemPot.cpp checkConsistency.cpp
        lsmsClass.cpp calculateEvec.cpp initializeAtom.cpp mixing.cpp
        ReplicaExchangeWL.cpp AlloyBankIO.cpp rotateToGlobal.cpp
        write_restart.cpp
        mixing_params.hpp
        mixing_params.cpp
        )

target_sources(lsmscore
        PUBLIC
        calculateChemPot.hpp
        calculateDensities.hpp
        calculateEvec.hpp
        checkConsistency.hpp
        EnergyContourIntegration.hpp
        EvecGenerator.h
        ExhaustiveIsing.h
        global.hpp
        Graph1d.hpp
        Graph1dMoments.hpp
        Graph1dMoments_new.hpp
        Graph1d_new.hpp
        Graph2d.hpp
        HDF5io.hpp
        initializeAtom.hpp
        lsmsClass.hpp
        LSMSMode.hpp
        mixing.hpp
        PotentialIO.hpp
        random_evec.h
        ReplicaExchangeWL.hpp
        SystemParameters.hpp
        WangLandau2d.h
        WangLandau.h
        WangLandau.hpp
        WangLandau_REWL.h
        write_restart.hpp
        )

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)


#
#
#
target_sources(lsms_main
        PUBLIC
        lsms.cpp)




