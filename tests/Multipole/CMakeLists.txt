

#add_must_test(lsmscore multipole_madlung)


#
# Test for replaceing the madelung structure
#

add_executable(test_madelung_improved)

target_sources(
        test_madelung_improved
        PUBLIC
        madelung_improved.cpp
)

target_link_libraries(
        test_madelung_improved
        PUBLIC
        multipole_madelung)


#
#
#

#add_executable(test_multipole_madelung)


#target_sources(
#        test_multipole_madelung
#        PUBLIC
#        multipole_madlung.cpp
#)

#target_link_libraries(
#        test_multipole_madelung
#        PUBLIC
#        multipole_madelung
#)
