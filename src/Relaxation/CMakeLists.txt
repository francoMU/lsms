
target_sources(lsmscore
        PUBLIC
        relax_params.cpp
        relax_params.hpp
        )


#
#
#

target_sources(lsms_relaxation_main
        PUBLIC
        relax_params.cpp
        lsms_relaxation_main.cpp
        lsms_relaxation.cpp
        lsms_relaxation.hpp
        lsms_relaxation_function.cpp
        lsms_relaxation_function.hpp
        nl_optimization.cpp
        nl_optimization.hpp
        base_function.hpp
        utils.cpp
        utils.hpp
        )


#
#
#

add_library(relaxation SHARED)

target_sources(relaxation
        PUBLIC
        nl_optimization.cpp
        nl_optimization.hpp
        base_function.hpp
        rosenbrockFunction.cpp
        rosenbrockFunction.hpp
        )

target_include_directories(relaxation
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )




