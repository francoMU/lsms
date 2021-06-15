


target_sources(lsmscore
        PUBLIC
        LuaInterface.cpp
        LuaSupport.hpp
        LuaSupport.cpp
        )

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)