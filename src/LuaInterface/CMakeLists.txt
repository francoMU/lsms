
target_sources(
  lsmscore
  PUBLIC LuaInterface.cpp
         LuaInterface.hpp
         LuaSupport.hpp
         LuaSupport.cpp
         RadialGrid_lua.hpp
         RadialGrid_lua.cpp
         RadialPotential_lua.cpp
         RadialPotential_lua.hpp)

target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
