

#
# Testing purposes
#

if (BUILD_MADELUNG_TEST)
    target_sources(
            multipole_madelung
            PUBLIC
            integer_factors.f90
            madelung_module.f90
    )

endif ()