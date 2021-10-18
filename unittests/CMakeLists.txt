
###################################################################
# Google tests
add_executable(test_multipole_moments test_multipole_moments.cpp)
target_link_libraries(test_multipole_moments gtest_main)

include(GoogleTest)
gtest_discover_tests(test_multipole_moments)
###################################################################


###################################################################
# 1. Test Lua Libraries
add_executable(test_lua_library test_lua_library.cpp)
target_link_libraries(test_lua_library PUBLIC Lua::Lua)
add_test(
        NAME test_lua_library
        COMMAND $<TARGET_FILE:test_lua_library>
)
###################################################################
